//! Exact elliptical-arc construction through the public crate API.

use hypercurve::{elliptical_arc_path, elliptical_span, EllipseMap2, Point2, Real};

fn r(x: i8) -> Real {
    Real::from(x)
}

fn ratio(num: i8, den: i8) -> Real {
    (Real::from(num) / Real::from(den)).unwrap()
}

/// Rotation by a 3-4-5 Pythagorean angle keeps `cos`/`sin` exact rationals, so
/// the whole ellipse stays exact (no `f64` boundary).
fn rotated_map(a: i8, b: i8) -> EllipseMap2 {
    let cos = ratio(3, 5); // cos θ
    let sin = ratio(4, 5); // sin θ
    let a = Real::from(a);
    let b = Real::from(b);
    // columns a·û = (a cos, a sin), b·v̂ = b·(-sin, cos)
    [&a * &cos, &(-&sin) * &b, &a * &sin, &b * &cos]
}

/// Axis-aligned exact unit-circle samples for a full CCW turn (four quarters).
fn full_turn_samples() -> Vec<(Real, Real)> {
    vec![
        (r(1), r(0)),
        (r(0), r(1)),
        (r(-1), r(0)),
        (r(0), r(-1)),
        (r(1), r(0)),
    ]
}

#[test]
fn full_rotated_ellipse_is_closed_with_four_spans() {
    let center = Point2::new(r(2), r(-1));
    let path = elliptical_arc_path(&center, &rotated_map(5, 2), &full_turn_samples()).unwrap();
    assert_eq!(path.curves().len(), 4);
    assert_eq!(path.start(), path.end());
}

#[test]
fn rotated_arc_start_endpoint_is_exact() {
    // First quarter of a rotated ellipse: u0 = (1,0) maps to a·û.
    // With a = 5 and (cos,sin) = (3/5, 4/5): a·û = (3, 4).
    let center = Point2::new(r(2), r(-1));
    let span =
        elliptical_span(&center, &rotated_map(5, 2), (&r(1), &r(0)), (&r(0), &r(1))).unwrap();
    let cps = span.control_points();
    assert_eq!(cps[0], &Point2::new(r(5), r(3))); // center (2,-1) + (3,4)
}

#[test]
fn partial_arc_is_open() {
    // A single quarter-span open arc (two samples -> one span).
    let center = Point2::new(r(0), r(0));
    let samples = vec![(r(1), r(0)), (r(0), r(1))];
    let path = elliptical_arc_path(&center, &rotated_map(4, 3), &samples).unwrap();
    assert_eq!(path.curves().len(), 1);
    assert_ne!(path.start(), path.end());
}
