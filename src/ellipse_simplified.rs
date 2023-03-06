use std::error::Error;


pub type Vec2d = nalgebra::Vector2<f64>;

#[derive(Debug, Clone, Copy)]
pub struct SimplifiedEllipse {
    pub radii: (f64, f64),
    pub x_axis_rotation: f64,
}

impl Default for SimplifiedEllipse {
    fn default() -> Self {
        SimplifiedEllipse {
            radii: (0.0, 0.0),
            x_axis_rotation: 0.0,
        }
    }
}


/// Parameters for the fitting algorithm.
pub struct FitParams {
    /// The subdivisions used to find the best x-axis rotation.
    /// 
    /// # Example
    /// If [`FitParams::angle_subdivisions`] is [300, 100, 50], the algorithm will
    /// divide the 360° arc in 300 subdivisions, then it will find the best
    /// angle in the 300 subdivisions, then it will divide the 360° arc in 100
    /// subdivisions, then it will find the best angle in the 100 subdivisions,
    /// and so on.
    pub angle_subdivisions: Vec<u32>,
}

impl Default for FitParams {
    fn default() -> Self {
        FitParams {
            angle_subdivisions: vec![300, 100],
        }
    }
}

impl SimplifiedEllipse {

    /// Fit an ellipse which center is the origin
    pub fn fit(points: impl Iterator<Item = Vec2d> + Clone, params: &FitParams) -> Option<Self> {
        let mut ellipse = SimplifiedEllipse::default();

        ellipse.x_axis_rotation = find_best_angle(
            points.clone(),
            &params.angle_subdivisions
        )
            .ok()?;

        let (a, b) = fit_aligned_ellipse_axes(
            points
                .map(|p| rotate_point(p, -ellipse.x_axis_rotation))
        );

        ellipse.radii.0 = a;
        ellipse.radii.1 = b;

        Some(ellipse)
    }

    /// Refit the ellipse with the given points, keeping the same x-axis rotation.
    pub fn refit_radii(mut self, points: impl Iterator<Item = Vec2d>) -> Self {
        let (a, b) = fit_aligned_ellipse_axes(
            points
                .map(|p| rotate_point(p, -self.x_axis_rotation))
        );

        self.radii.0 = a;
        self.radii.1 = b;

        self
    }

    /// Returns an iterator over the points of the ellipse discretized in `n` points.
    pub fn discretize(&self, n: u32, repeat_last: bool) -> impl Iterator<Item=Vec2d> {
        let (a, b) = self.radii;
        let angle_step = 2.0 * std::f64::consts::PI / (n as f64);
        let x_axis_rotation = self.x_axis_rotation;

        let range = if repeat_last { 0..n + 1} else { 0..n };

        range.into_iter().map(move |i| {
            let angle = i as f64 * angle_step;

            let x = a * angle.cos();
            let y = b * angle.sin();

            let pt = Vec2d::new(x, y);

            rotate_point(pt, x_axis_rotation)
        })
    }

    /// Evaluate the ellipse equation at the given point.
    pub fn eval(&self, pt: Vec2d) -> f64 {
        let pt = rotate_point(pt, -self.x_axis_rotation);

        let a = 1.0 / (self.radii.0 * self.radii.0);
        let b = 1.0 / (self.radii.1 * self.radii.1);

        a * pt.x * pt.x + b * pt.y * pt.y - 1.0
    }

    /// Returns `true` if the given point is inside the ellipse.
    pub fn contains(&self, pt: Vec2d) -> bool {
        self.eval(pt) <= 0.0
    }
}

/// # Method
/// 
/// The method used is the same as [`crate::conics2d::fit_simplified_2d_conic_section`] but we
/// use the following conic section formula:
/// $$
/// C(p_i) = a \\, x_i^2 + b \\, y_i^2
/// $$
/// 
/// The matrix $A$ is then:
/// $$
///   A = \\begin{bmatrix}
///     x_1^2 & y_1^2 \\\\
///     x_2^2 & y_2^2 \\\\
///     \\vdots & \\vdots \\\\
///     x_n^2 & y_n^2
///   \\end{bmatrix}
/// $$
/// and the vector $p$:
/// $$
///  p = \\begin{bmatrix}
///    a \\\\
///    b
/// \\end{bmatrix}
/// $$
/// 
/// # Returns
/// the ellipse radii.
pub fn fit_aligned_ellipse_axes(aligned_pts: impl Iterator<Item=Vec2d>) -> (f64, f64) {
    let (matrix_m, vector_v) = {
        let mut matrix_m = nalgebra::SMatrix::<f64, 2, 2>::zeros();
        let mut vector_v = nalgebra::SVector::<f64, 2>::zeros();
        for p in aligned_pts {
            let data = nalgebra::SVector::<f64, 2>::new(p.x * p.x, p.y * p.y);
            matrix_m += data * data.transpose();
            vector_v += data;
        }
        (matrix_m, vector_v)
    };

    // nalgebra way: let p_solution = matrix_m.lu().solve(&vector_v);
    let p_solution = solve_2x2_linear_system(&matrix_m, &vector_v);

    if let Some(p_solution) = p_solution {
        (1.0 / p_solution[0].sqrt(), 1.0 / p_solution[1].sqrt())
    } else {
        (0.0, 0.0)
    }
}

/// Solve the following system:
/// $$
/// M \\, p = v
/// $$
/// using Cramer's rule.
/// 
/// This function is a simplified version of [`nalgebra::SMatrix::solve`].
/// It only works for 2x2 matrices but it is quite simple so it can be used as a
/// reference.
pub fn solve_2x2_linear_system(matrix_m: &nalgebra::Matrix2<f64>, vector_v: &Vec2d) -> Option<Vec2d> {
    let det = matrix_m[(0, 0)] * matrix_m[(1, 1)] - matrix_m[(0, 1)] * matrix_m[(1, 0)];

    // The matrix is not invertible.
    if det == 0.0 {
        return None;
    }

    let inv_det = 1.0 / det;

    let p_solution = Vec2d::new(
        inv_det * (matrix_m[(1, 1)] * vector_v.x - matrix_m[(0, 1)] * vector_v.y),
        inv_det * (matrix_m[(0, 0)] * vector_v.y - matrix_m[(1, 0)] * vector_v.x),
    );

    Some(p_solution)
}

pub fn find_best_angle(pts: impl Iterator<Item=Vec2d> + Clone, angle_subdivisions: &[u32]) -> Result<f64, Box<dyn Error>> {
    if angle_subdivisions.is_empty() {
        return Err("angle_subdivisions is empty".into());
    }

    if pts.clone().count() < 3 {
        return Err("pts must have at least 2 points".into());
    }

    let score_for_angle = |angle: f64| -> Result<f64, Box<dyn Error>> {
        let rotated_pts = pts.clone().map(|p| rotate_point(p, -angle));
        let (var_x, var_y) = compute_variance2(rotated_pts, Some(Vec2d::new(0.0, 0.0)));

        if var_y <= 0.0 {
            return Err("var_y is 0".into());
        }

        Ok(var_x / var_y)
    };

    let mut range = (-std::f64::consts::PI / 2.0, std::f64::consts::PI / 2.0);
    for &subdivisions in angle_subdivisions {
        let step = (range.1 - range.0) / subdivisions as f64;
        let mut best_score = 0.0;
        let mut best_angle = 0.0;
        for i in 0..subdivisions {
            let angle = range.0 + step * i as f64;
            let score = score_for_angle(angle)?;
            if score > best_score {
                best_score = score;
                best_angle = angle;
            }
        }
        range = (best_angle - step, best_angle + step);
    }

    Ok((range.0 + range.1) / 2.0)
}

/// Compute the variance of a vector.
pub fn compute_variance(values_iter: impl Iterator<Item=f64> + Clone, mean: Option<f64>) -> f64 {
    let mean = if let Some(mean) = mean {
        mean
    } else {
        values_iter.clone().sum::<f64>() / values_iter.clone().count() as f64
    };
    values_iter.clone().map(|v| (v - mean).powi(2)).sum::<f64>() / values_iter.count() as f64
}

/// Compute the variance of each component of a 2d vector.
pub fn compute_variance2(values_iter: impl Iterator<Item=Vec2d> + Clone, center: Option<Vec2d>) -> (f64, f64) {
    (
        compute_variance(
            values_iter.clone().map(|v| v.x),
            center.map(|c| c.x)
        ),
        compute_variance(
            values_iter.clone().map(|v| v.y),
            center.map(|c| c.y)
        )
    )
}

/// Rotate a point by an angle.
pub fn rotate_point(point: Vec2d, angle: f64) -> Vec2d {
    let (sin, cos) = angle.sin_cos();
    Vec2d::new(
        point.x * cos - point.y * sin,
        point.x * sin + point.y * cos,
    )
}

/// Square of a number.
pub fn sqr<T: std::ops::Mul<Output = T> + Copy>(x: T) -> T {
    x * x
}