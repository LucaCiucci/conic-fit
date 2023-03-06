use std::ops::Range;

use conic_fit::conics2d::{Ellipse2d, Vec2d, fit_simplified_2d_conic_section};
use plotters::{style::WHITE, prelude::{ChartBuilder, Circle}};

use plotters::drawing::IntoDrawingArea;
use plotters::style::Color;


pub fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Hello, world!");

    let target_ellipse = Ellipse2d {
        center: Vec2d::new(1.0, 0.0),
        radii: Vec2d::new(2.0, 1.0),
        rotation: 0.0,
    };

    let rotation_matrix = nalgebra::Matrix2::<f64>::new(
        target_ellipse.rotation.cos(),
        -target_ellipse.rotation.sin(),
        target_ellipse.rotation.sin(),
        target_ellipse.rotation.cos(),
    );

    let n = 10;

    let mut points = Vec::new();
    for i in 0..n {
        let t = i as f64 / (n as f64) * 2.0 * std::f64::consts::PI;
        let p = Vec2d::new(t.cos() * target_ellipse.radii.x, t.sin() * target_ellipse.radii.y);
        points.push(rotation_matrix * p + target_ellipse.center);
    }

    //println!("Points: {:?}", points);

    let conic = fit_simplified_2d_conic_section(&points).unwrap();

    println!("Conic: {:?}", conic);
    println!("Conic is ellipse: {:?}", conic.is_ellipse());
    println!("Conic center: {:?}", conic.center());
    println!("Conic exes: {:?}", conic.axes());
    println!("Conic exis 1: {:?}", conic.axes().unwrap().ellipse_axis_1_length());
    println!("Conic exis 2: {:?}", conic.axes().unwrap().ellipse_axis_2_length());
    println!("Conic centered: {:?}", conic.translate(-conic.center()));
    println!("Conic centered: {:?}", conic.translate(-conic.center()).axes());
    println!("Conic centered: {:?}", conic.translate(-conic.center()).axes().unwrap().ellipse_axis_1_length());
    println!("Conic centered: {:?}", conic.translate(-conic.center()).axes().unwrap().ellipse_axis_2_length());

    let root = plotters::backend::BitMapBackend::new("a.png", (1024 / 2, 768 / 2)).into_drawing_area();

    root.fill(&WHITE).unwrap();

    let mut chart = ChartBuilder::on(&root).build_cartesian_2d::<Range<f64>, Range<f64>>(-5.0..5.0, -5.0..5.0)?;

    chart.draw_series(
        points.iter().map(|pt| {
            Circle::new((pt.x, pt.y), 1, plotters::style::RED.filled())
        })
    )?;

    Ok(())
}