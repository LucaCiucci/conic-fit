
/// 2d column vector of doubles
pub type Vec2d = nalgebra::Vector2<f64>;

/// A conic section in 2d space in the form $a \\, x^2 + b \\, xy + c \\, y^2 + d \\, x + e \\, y - 1 = 0$
#[derive(Clone, Copy, Debug)]
pub struct Simplified2dConicSection {
    pub coefficients: nalgebra::SVector<f64, 5>,
}

impl Simplified2dConicSection {

    pub fn from_coefficients(coefficients: nalgebra::SVector<f64, 5>) -> Self {
        Self {
            coefficients
        }
    }

    pub fn from_slice(coefficients: &[f64; 5]) -> Self {
        Self {
            coefficients: nalgebra::SVector::from_iterator(coefficients.iter().cloned())
        }
    }

    /// coefficient for $x^2$
    pub fn a(&self) -> f64 {
        self.coefficients[0]
    }

    /// coefficient for $x^2$
    pub fn a_mut(&mut self) -> &mut f64 {
        &mut self.coefficients[0]
    }

    /// coefficient for $xy$
    pub fn b(&self) -> f64 {
        self.coefficients[1]
    }

    /// coefficient for $xy$
    pub fn b_mut(&mut self) -> &mut f64 {
        &mut self.coefficients[1]
    }

    /// coefficient for $y^2$
    pub fn c(&self) -> f64 {
        self.coefficients[2]
    }

    /// coefficient for $y^2$
    pub fn c_mut(&mut self) -> &mut f64 {
        &mut self.coefficients[2]
    }

    /// coefficient for $x$
    pub fn d(&self) -> f64 {
        self.coefficients[3]
    }

    /// coefficient for $x$
    pub fn d_mut(&mut self) -> &mut f64 {
        &mut self.coefficients[3]
    }

    /// coefficient for $y$
    pub fn e(&self) -> f64 {
        self.coefficients[4]
    }

    /// coefficient for $y$
    pub fn e_mut(&mut self) -> &mut f64 {
        &mut self.coefficients[4]
    }

    /// The number of coefficients in the conic section
    pub fn n_coefficients() -> usize {
        5
    }

    /// returns the value of the conic section at the point `p`:
    /// $$
    ///     a \\, x^2 + b \\, xy + c \\, y^2 + d \\, x + e \\, y =
    ///     \\begin{bmatrix}
    ///        a & b & c & d & e
    ///     \\end{bmatrix}
    ///     \\begin{bmatrix}
    ///        x^2 \\\\
    ///        xy \\\\
    ///        y^2 \\\\
    ///        x \\\\
    ///        y
    ///     \\end{bmatrix}
    /// $$
    pub fn eval(&self, p: nalgebra::Vector2<f64>) -> f64 {
        let data = Self::point_multiplier_vec(p);
        self.coefficients.dot(&data) - 1.0
    }

    /// takes a point and returns the column vector:
    /// $$
    /// \\begin{bmatrix}
    ///     x^2 \\\\
    ///     xy \\\\
    ///     y^2 \\\\
    ///     x  \\\\
    ///     y \\\\
    /// \\end{bmatrix}
    /// $$
    pub fn point_multiplier_vec(p: nalgebra::Vector2<f64>) -> nalgebra::SVector<f64, 5> {
        nalgebra::SVector::from_column_slice(&[
            sqr(p.x),
            p.x * p.y,
            sqr(p.y),
            p.x,
            p.y
        ])
    }

    /// An ellipse is defined as a conic section with a positive determinant.
    /// 
    /// The quadratic part of the conic section is:
    /// $$
    ///   a \\, x^2 + b \\, xy + c \\, y^2 = 
    ///   \\begin{bmatrix}
    ///     x & y
    ///   \\end{bmatrix}
    ///   \\begin{bmatrix}
    ///     a & b/2 \\\\
    ///     b/2 & c
    ///   \\end{bmatrix}
    ///   \\begin{bmatrix}
    ///     x \\\\
    ///     y
    ///   \\end{bmatrix}
    /// $$
    /// This product is positive if and only if the determinant of the quadratic part is positive.
    pub fn is_ellipse(&self) -> bool {
        let a = self.a();
        let b = self.b();
        let c = self.c();
        let determinant = a * c - sqr(b) / 4.0;
        determinant > 0.0
    }

    /// The center of a conic section is the intersection of the two lines that are perpendicular to the major and minor axes.
    /// 
    /// # Method
    /// 
    /// Given the "center" of the conic section, $C = (x_0, y_0)$, the conic equation can be written as:
    /// $$
    ///    a_R \\, (x - x_0)^2 + b_R \\, (x - x_0) \\, (y - y_0) + c_R \\, (y - y_0)^2 - 1 = 0
    /// $$
    /// where the "$R$" subscript indicates that the coefficients represents a conic whose axis are rotated a certain by $\\theta$.
    /// 
    /// As we can see, the jacobian of the conic equation with respect to $x$ and $y$ is 0 at the center of the conic section,
    /// therefore, we can find the center by solving the following system of equations:
    /// $$
    /// \\begin{cases}
    ///   \\partial_x C(x_0, y_0) = 0 \\\\
    ///   \\partial_y C(x_0, y_0) = 0
    /// \\end{cases}
    /// $$
    /// That is (see [`Simplified2dConicSection::eval`]):
    /// $$
    /// \\begin{cases}
    ///   2 \\, a \\, x_0 + b \\, y_0 + d = 0 \\\\
    ///   b \\, x_0 + 2 \\, c \\, y_0 + e = 0
    /// \\end{cases}
    /// $$
    /// And in matrix form:
    /// $$
    /// \\begin{bmatrix}
    ///   2 \\, a & b \\\\
    ///   b & 2 \\, c
    /// \\end{bmatrix}
    /// \\begin{bmatrix}
    ///   x_0 \\\\
    ///   y_0
    /// \\end{bmatrix} = \\begin{bmatrix}
    ///   -d \\\\
    ///   -e
    /// \\end{bmatrix}
    /// $$
    /// and solving for $x_0$ and $y_0$:
    /// $$
    /// \\begin{bmatrix}
    ///   x_0 \\\\
    ///   y_0
    /// \\end{bmatrix} = \\frac{1}{2 \\, a \\, c - b^2}
    /// \\begin{bmatrix}
    ///   c \\, (-d) - b \\, (-e) \\\\
    ///   a \\, (-e) - b \\, (-d)
    /// \\end{bmatrix}
    /// $$
    /// TODO this last equation and the implementation is written by copilot, it might be wrong, to be checked
    pub fn center(&self) -> Vec2d {
        let a = self.a();
        let b = self.b();
        let c = self.c();
        let d = self.d();
        let e = self.e();

        let matricial_form = nalgebra::Matrix2::new(
            2.0 * a, b,
            b, 2.0 * c
        );

        let rhs = nalgebra::Vector2::new(-d, -e);
        return rhs;

        let x = -(c * d - b * e) / (2.0 * a * c - sqr(b));
        let y = -(a * e - b * d) / (2.0 * a * c - sqr(b));

        Vec2d::new(x, y)
    }

    fn translation_gamma(&self, offset: &Vec2d) -> f64 {
        -(
            self.a() * sqr(offset.x) +
            self.b() * offset.x * offset.y +
            self.c() * sqr(offset.y) +
            -self.d() * offset.x +
            -self.e() * offset.y +
            - 1.0
        )
    }

    /// TODO calculations
    pub fn translate(&self, offset: Vec2d) -> Self {
        let gamma = self.translation_gamma(&offset);
        println!("translation gamma: {}", gamma);
        Self::from_slice(&[
            self.a() / gamma,
            self.b() / gamma,
            self.c() / gamma,
            (self.d() - 2.0 * self.a() * offset.x - self.b() * offset.y) / gamma,
            (self.e() - 2.0 * self.c() * offset.y - self.b() * offset.x) / gamma,
        ])
    }

    /// If the conic has its axis aligned in a primed coordinate system, in this c.s we can write
    /// the conic section equation in the form:
    /// $$
    ///   a' x'^2 + b' y'^2 = 1
    /// $$
    /// The primed system is defined:
    /// $$
    /// p = \\begin{bmatrix}
    ///   cos(\\theta) & -sin(\\theta) \\\\
    ///   sin(\\theta) & cos(\\theta)
    /// \\end{bmatrix} \\, p'
    /// \\Longleftrightarrow
    /// p' = \\begin{bmatrix}
    ///  cos(\\theta) & sin(\\theta) \\\\
    ///  -sin(\\theta) & cos(\\theta)
    /// \\end{bmatrix} \\, p
    /// $$
    /// So, the conic equation in the unprimed system is:
    /// $$
    /// a' (x \\, \\cos(\\theta) + y \\, \\sin(\\theta))^2 + b' (-x \\, \sin(\\theta) + y \\, \\cos(\\theta))^2 = 1
    /// $$
    /// That can be expanded:
    /// $$
    /// (a' \\, \\cos^2(\\theta) + b' \\, \\sin^2(\\theta)) \\, x^2 + 2 \\, (a' - b') \\, \\sin(\\theta) \\, \\cos(\\theta) \\, x \\, y + (a' \\, \\sin^2(\\theta) + b' \\, \\cos^2(\\theta)) \\, y^2 = 1
    /// $$
    /// That is:
    /// $$
    /// a \\, x^2 + b \\, xy + c \\, y^2 = 1
    /// $$
    /// where:
    /// $$
    /// \\begin{cases}
    ///   a = a' \\, \\cos^2(\\theta) + b' \\, \\sin^2(\\theta) \\\\
    ///   b = 2 (a' - b') \\, \\sin(\\theta) \\, \\cos(\\theta) \\\\
    ///   c = a' \\, \\sin^2(\\theta) + b' \\, \\cos^2(\\theta)
    /// \\end{cases}
    /// $$
    /// We now want to find the angle $\\theta$ is the rotation of the primed system with respect to the unprimed system. To do that, we can consider the quantity:
    /// $$
    /// a - c = a' (\\cos^2(\\theta) - \\sin^2(\\theta)) + b' (\\sin^2(\\theta) - \\cos^2(\\theta)) = (a' - b') \\, \\sin(2 \\, \\theta)
    /// $$
    /// white the therm $b$ can be written as:
    /// $$
    /// b = 2 (a' - b') \\, \\sin(\\theta) \\, \\cos(\\theta) = (a' - b') \\, \\sin(2 \\, \\theta)
    /// $$
    /// Thus, we can find the angle $\\theta$ as:
    /// $$
    /// \\theta = \\frac{1}{2} \\, \\text{atan2}(b, a - c)
    /// $$
    /// 
    /// We can now use this angle to find the length of the major and minor axis of the ellipse:
    /// $$
    /// \\begin{cases}
    ///   a = a' \\, \\cos^2(\\theta) + b' \\, \\sin^2(\\theta) \\\\
    ///   c = a' \\, \\sin^2(\\theta) + b' \\, \\cos^2(\\theta)
    /// \\end{cases}
    /// $$
    /// $$
    /// \begin{bmatrix}
    ///   a \\\\
    ///   c
    /// \end{bmatrix} = \begin{bmatrix}
    ///   \\cos^2(\\theta) & \\sin^2(\\theta) \\\\
    ///   \\sin^2(\\theta) & \\cos^2(\\theta)
    /// \end{bmatrix} \\, \begin{bmatrix}
    ///    a' \\\\
    ///    c'
    /// \end{bmatrix}
    /// $$
    /// and this can be solved to find the length of the major and minor axis.
    /// TODO soluzione ottimizzata e senza fare sin e cos di theta due volte
    /// perchè conosciamo già il seno e il coseno di theta
    pub fn axes(&self) -> Option<ConicAxes> {
        let a = self.a();
        let b = self.b();
        let c = self.c();

        let theta = 0.5 * b.atan2(a - c);
        
        let matrix = nalgebra::Matrix2::new(
            sqr(theta.cos()), sqr(theta.sin()),
            sqr(theta.sin()), sqr(theta.cos())
        );

        let axes = matrix * nalgebra::Vector2::new(a, c);
        let gamma = self.translation_gamma(&self.center());
        //println!("axes: {:?}, gamma: {}", axes, gamma);
        //let axes = axes / self.translation_gamma(&self.center());

        if axes[0] < 0.0 || axes[1] < 0.0 {
            return None;
        }

        if axes[1] < axes[0] {
            // swap and fix the angle
            return Some(ConicAxes {
                a_coeff: axes.y,
                b_coeff: axes.x,
                major_x_axis_angle: theta + std::f64::consts::FRAC_PI_2,
            });
        }

        return Some(ConicAxes {
            a_coeff: axes.x,
            b_coeff: axes.y,
            major_x_axis_angle: theta,
        });
    }

    /// returns the quadratic form part
    pub fn quadratic_form(&self) -> nalgebra::Matrix2<f64> {
        let a = self.a();
        let b = self.b();
        let c = self.c();
        let d = self.d();
        let e = self.e();

        nalgebra::Matrix2::new(
            a, b / 2.0,
            b / 2.0, c
        )
    }
}

#[derive(Debug, Clone, Copy)]
pub struct ConicAxes {
    pub a_coeff: f64,
    pub b_coeff: f64,
    pub major_x_axis_angle: f64,
}

impl ConicAxes {
    pub fn ellipse_axis_1_length(&self) -> f64 {
        (1.0 / self.a_coeff).sqrt()
    }

    pub fn ellipse_axis_2_length(&self) -> f64 {
        (1.0 / self.b_coeff).sqrt()
    }
}

pub trait Conic2d {
    fn to_simplified_conic(&self) -> Simplified2dConicSection;
    fn from_simplified_conic(conic: &Simplified2dConicSection) -> Self;
}

/// An ellipse in 2d space
pub struct Ellipse2d {
    /// The center of the ellipse
    pub center: Vec2d,

    /// The radii of the ellipse
    pub radii: Vec2d,

    /// The rotation of the ellipse in radians
    pub rotation: f64,
}


/// Fits a conic section to a set of points, using the method described below.
/// If the number of points is less than the number of coefficients in the conic section, then `None` is returned.
/// 
/// # Method
/// We can approximate the residual of each point $p_i$ from the conic section $C$ as:
/// $$
/// r_i = 1 - C(p_i)
/// $$
/// where $C(p_i)$ is the value of the conic section at the point $p_i$:
/// $$
/// C(p_i) = a \\, x_i^2 + b \\, x_i y_i + c \\, y_i^2 + d \\, x_i + e \\, y_i
/// $$
/// 
/// We can then minimize the sum of the squares of the residuals:
/// $$
///  \\chi^2 = \\sum_{i=1}^n r_i^2 = \\sum_{i=1}^n (1 - C(p_i))^2
/// $$
/// 
/// We can now define the matrix $A$:
/// $$
///   A = \\begin{bmatrix}
///     x_1^2 & x_1 y_1 & y_1^2 & x_1 & y_1 \\\\
///     x_2^2 & x_2 y_2 & y_2^2 & x_2 & y_2 \\\\
///     \\vdots & \\vdots & \\vdots & \\vdots & \\vdots \\\\
///     x_n^2 & x_n y_n & y_n^2 & x_n & y_n
///   \\end{bmatrix}
/// $$
/// and the vector $p$:
/// $$
///  p = \\begin{bmatrix}
///    a \\\\
///    b \\\\
///    c \\\\
///    d \\\\
///    e
/// \\end{bmatrix}
/// $$
/// 
/// We can then rewrite the residuals as:
/// $$
///  r = b - A p
/// $$
/// where $b$ is a vector of ones with the same number of rows as $A$ (i.e. $n$).
/// The $\\chi^2$ can then be written as:
/// $$
/// \\chi^2 = r^T r = (b - A p)^T (b - A p)
/// $$
/// 
/// Now, the minimum of $\\chi^2$ is a stationary point of $\\chi^2$ with respect to $p$, thus
/// we can take the derivative of $\\chi^2$ with respect to $p$:
/// $$
/// \\frac{\\partial \\chi^2}{\\partial p} = 2 (Ap - b)^T A = 2 (p^T A^T - b^T) A
/// $$
/// This has to be zero for a stationary point, thus:
/// $$
/// (p^T A^T - b^T) A = 0 \\Longleftrightarrow p^T A^T A = b^T A
/// $$
/// We now take the transpose of both sides in order to get a nicer matrix equation in $p$:
/// $$
/// A^T A p = A^T b
/// $$
/// This system can be numerically explicited and solved for $p$, but we don't want to deal with
/// unnecessarily large matrices, so we can write what these pieces are.
/// 
/// Let's define the matrix $M$:
/// $$
/// M = A^T A
/// $$
/// $$
/// M _ {ij} = (A^T) _ {ik} A _ {kj} = A _ {ki} A _ {kj}
/// $$
/// And this is a 5x5 matrix.
/// 
/// Next, we define the vector $v$:
/// $$
/// v = A^T b
/// $$
/// $$
/// v _ {i} = (A^T) _ {ik} b _ {k} = A _ {ki} b _ {k}
/// $$
/// And this is a 5x1 vector (5 components column vector).
/// 
/// The last step to solve the system is to solve the matrix equation $M p = v$ and this can be done
/// using the LU decomposition.
pub fn fit_simplified_2d_conic_section(pts: &[Vec2d]) -> Option<Simplified2dConicSection> {

    if pts.len() < Simplified2dConicSection::n_coefficients() {
        return None;
    }

    let (matrix_m, vector_v) = {
        let mut matrix_m = nalgebra::SMatrix::<f64, 5, 5>::zeros();
        let mut vector_v = nalgebra::SVector::<f64, 5>::zeros();
        for p in pts {
            let data = Simplified2dConicSection::point_multiplier_vec(*p);
            matrix_m += data * data.transpose();
            vector_v += data;
        }
        (matrix_m, vector_v)
    };

    // solve for Mp = v for the coefficients p
    let p_solution = matrix_m.lu().solve(&vector_v);

    if let Some(p) = p_solution {
        Some(Simplified2dConicSection::from_coefficients(p))
    } else {
        None
    }
}

/// Square of a number.
pub fn sqr<T: std::ops::Mul<Output = T> + Copy>(x: T) -> T {
    x * x
}