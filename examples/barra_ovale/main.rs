

const IN_FILE: &str = "barra-ovale.asc";

use std::{fs::File, io::{BufRead, BufReader}};

use conic_fit::ellipse_simplified::*;

use plotters::{drawing::IntoDrawingArea, series::LineSeries, prelude::PathElement};
use plotters::prelude::{ChartBuilder, Circle};
use plotters::style::{WHITE, Color};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Hello, world!");

    let points = load_pts();

    println!("{:} points loaded", points.len());

    let ell = SimplifiedEllipse::fit(
        points.iter().map(|p| *p),
        &Default::default()
    )
        .unwrap();

    let filtered_pts = points.iter().filter(|p| ell.contains(**p)).map(|p| *p).collect::<Vec<_>>();

    let ell2 = ell.clone().refit_radii(filtered_pts.iter().map(|p| *p));

    let root = plotters::backend::BitMapBackend::new("a.png", (1024, 768)).into_drawing_area();

    root.fill(&WHITE).unwrap();

    let r = 11.0;
    let range = -r..r;
    let mut chart = ChartBuilder::on(&root)
        .set_all_label_area_size(20)
        .build_cartesian_2d(range.clone(), range)?;

    chart.draw_series(
        points.iter().map(|pt| {
            Circle::new((pt.x, pt.y), 1, plotters::style::RED.filled())
        })
    )?;

    chart.draw_series(
        LineSeries::new(ell.discretize(1000, true).map(|pt| (pt.x, pt.y)), &plotters::style::BLACK)
    )?
        .label("fit tutti")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &plotters::style::BLACK));

    chart.draw_series(
        LineSeries::new(ell2.discretize(1000, true).map(|pt| (pt.x, pt.y)), &plotters::style::GREEN)
    )?
        .label("fit punti interni")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &plotters::style::GREEN));

    chart.configure_series_labels()
        .border_style(&plotters::style::BLACK)
        .position(plotters::prelude::SeriesLabelPosition::UpperRight)
        .label_font(("sans-serif", 20))
        .draw()?;

    Ok(())
}

fn load_pts() -> Vec<Vec2d> {
    let mut points: Vec<Vec2d> = Vec::new();

    let file = File::open(IN_FILE).expect("file not found");
    BufReader::new(file).lines().for_each(|l| {
        let line = l.unwrap();

        if line.starts_with("#") {
            return;
        }

        let mut parts = line.split_whitespace();

        let mut pt = Vec2d::new(0.0, 0.0);

        let part1 = parts.next().unwrap();
        if !part1.starts_with("X:") {
            return;
        }
        let x = part1[2..].parse::<f64>().unwrap();
        pt.x = x;

        let part2 = parts.next().unwrap();
        if !part2.starts_with("Y:") {
            return;
        }

        let y = part2[2..].parse::<f64>().unwrap();
        pt.y = y + 0.0 * x;

        let part3 = parts.next().unwrap();
        if !part3.starts_with("Z:") {
            return;
        }

        let _z = part3[2..].parse::<f64>().unwrap();

        points.push(rotate_point(pt, 0.1));
    });

    points
}