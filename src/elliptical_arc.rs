//! Exact elliptical-arc decomposition into rational quadratic Bezier spans.
//!
//! An ellipse is the affine image of a circle:
//! `E(φ) = center + M·(cos φ, sin φ)`, where `M = [a·û | b·v̂]` is the ellipse's
//! 2×2 linear map (semi-axes `a`,`b` along unit axis directions `û`,`v̂`).
//!
//! Because an affine map preserves rational-Bezier **weights** and maps control
//! nets linearly, each elliptical-arc span reuses the exact unit-circle
//! rational-quadratic weight `w = sqrt((1 + u0·u1)/2) = cos(Δ/2)` inherited from
//! the pre-image circular arc, and only the three control points are mapped
//! through `M`. This is the same rational-quadratic construction
//! [`crate::CircularArc2::rational_bezier_decomposition`] uses for circles
//! (`M = diag(r, r)`), generalized to a non-uniform / rotated map.
//!
//! Inputs are exact [`Real`] values. Callers that begin from an `f64` rotation or
//! sweep (for example meshup's `hcurve` bridge) perform the single trigonometric
//! conversion to `Real` sample directions before calling in, so this module stays
//! exact. Axis-aligned ellipses whose span boundaries fall on the axes
//! (`(±1, 0)`, `(0, ±1)`) are fully exact.

use std::cmp::Ordering;

use crate::{
    Curve2, CurveError, CurveFamily2, CurveGeometry2, CurveOperation2, CurvePath2,
    CurvePolicy, ExactCurveError, ExactCurveResult, Point2, RationalQuadraticBezier2, Real,
    UncertaintyReason,
};

/// The exact 2×2 linear map of an ellipse, row-major:
/// `M·(x, y) = (m[0]·x + m[1]·y, m[2]·x + m[3]·y)`.
///
/// For semi-axes `a`,`b` and unit axis directions `û`,`v̂` the columns are `a·û`
/// and `b·v̂`, i.e. `m = [a·û.x, b·v̂.x, a·û.y, b·v̂.y]`.
pub type EllipseMap2 = [Real; 4];

/// Applies the ellipse map `M` to a direction `(ux, uy)`.
fn apply_map(map: &EllipseMap2, ux: &Real, uy: &Real) -> (Real, Real) {
    (
        &(&map[0] * ux) + &(&map[1] * uy),
        &(&map[2] * ux) + &(&map[3] * uy),
    )
}

/// Builds one exact elliptical-arc span between two unit-circle directions.
///
/// `u0` and `u1` are `(cos φ, sin φ)` on (or near) the unit circle; the span
/// covers their sub-sweep, which must be strictly less than a half turn so the
/// affine-invariant middle weight `sqrt((1 + u0·u1)/2)` is positive.
pub fn elliptical_span(
    center: &Point2,
    map: &EllipseMap2,
    u0: (&Real, &Real),
    u1: (&Real, &Real),
) -> ExactCurveResult<RationalQuadraticBezier2> {
    // Affine-invariant middle weight from the pre-image circular arc:
    //   w^2 = (1 + u0·u1) / 2 = cos^2(Δ/2).
    let dot = &(u0.0 * u1.0) + &(u0.1 * u1.1);
    let weight_squared = ((Real::one() + &dot) / Real::from(2_i8))
        .map_err(|cause| span_invalid(cause.into()))?;
    match crate::classify::compare_reals(
        &weight_squared,
        &Real::zero(),
        &CurvePolicy::certified(),
    ) {
        Some(Ordering::Greater) => {},
        Some(_) => return Err(span_invalid(CurveError::InvalidArcSweep)),
        None => return Err(span_blocked(UncertaintyReason::RealSign)),
    }
    let weight = weight_squared
        .clone()
        .sqrt()
        .map_err(|cause| span_invalid(cause.into()))?;

    // Mapped endpoints E = center + M·u.
    let (e0x, e0y) = apply_map(map, u0.0, u0.1);
    let (e1x, e1y) = apply_map(map, u1.0, u1.1);
    let start = Point2::new(center.x() + &e0x, center.y() + &e0y);
    let end = Point2::new(center.x() + &e1x, center.y() + &e1y);

    // In circle space the control point's offset from center is
    // `(u0 + u1) / (2 w^2)` (the tangent intersection). `M` is linear, so the
    // mapped control offset is `(M·u0 + M·u1) / (2 w^2)`; the weight is unchanged.
    let control_denominator = Real::from(2_i8) * weight_squared.clone();
    let control_x =
        ((&e0x + &e1x) / &control_denominator).map_err(|cause| span_invalid(cause.into()))?;
    let control_y =
        ((&e0y + &e1y) / &control_denominator).map_err(|cause| span_invalid(cause.into()))?;
    let control = Point2::new(center.x() + &control_x, center.y() + &control_y);

    RationalQuadraticBezier2::try_unit_end_weights(start, control, end, weight)
        .map_err(span_invalid)
}

/// Builds an exact elliptical-arc (or full-ellipse) [`CurvePath2`].
///
/// `samples` are ordered unit-circle directions `(cos φ, sin φ)` at the span
/// boundaries; consecutive samples must be strictly less than a half turn apart.
/// Two samples yield one span. A closed full ellipse passes four (or more)
/// samples whose first and last are the same `Real` pair, so the emitted path
/// closes exactly.
pub fn elliptical_arc_path(
    center: &Point2,
    map: &EllipseMap2,
    samples: &[(Real, Real)],
) -> ExactCurveResult<CurvePath2> {
    if samples.len() < 2 {
        return Err(span_invalid(CurveError::InvalidCurveParameter));
    }
    let mut curves = Vec::with_capacity(samples.len() - 1);
    for pair in samples.windows(2) {
        let span = elliptical_span(
            center,
            map,
            (&pair[0].0, &pair[0].1),
            (&pair[1].0, &pair[1].1),
        )?;
        curves.push(Curve2::new(CurveGeometry2::RationalQuadraticBezier(span)));
    }
    CurvePath2::try_new(curves)
}

fn span_invalid(cause: CurveError) -> ExactCurveError {
    ExactCurveError::invalid(
        CurveOperation2::BezierDecomposition,
        CurveFamily2::RationalQuadraticBezier,
        cause,
    )
}

fn span_blocked(reason: UncertaintyReason) -> ExactCurveError {
    ExactCurveError::blocked(
        CurveOperation2::BezierDecomposition,
        CurveFamily2::RationalQuadraticBezier,
        reason,
    )
}

#[cfg(test)]
mod tests {
    use super::{elliptical_arc_path, elliptical_span, EllipseMap2};
    use crate::{CurvePolicy, Point2, Real};

    fn r(x: i8) -> Real {
        Real::from(x)
    }

    /// Axis-aligned exact unit-circle samples for a full turn (CCW quarters).
    fn quarter_samples() -> Vec<(Real, Real)> {
        vec![
            (r(1), r(0)),
            (r(0), r(1)),
            (r(-1), r(0)),
            (r(0), r(-1)),
            (r(1), r(0)),
        ]
    }

    /// Axis-aligned map `M = diag(a, b)`.
    fn diag(a: i8, b: i8) -> EllipseMap2 {
        [r(a), r(0), r(0), r(b)]
    }

    #[test]
    fn axis_aligned_full_ellipse_builds_four_closed_spans() {
        let center = Point2::new(r(0), r(0));
        let path = elliptical_arc_path(&center, &diag(3, 2), &quarter_samples()).unwrap();
        assert_eq!(path.curves().len(), 4);
        // Closed: last end coincides with first start.
        assert_eq!(path.start(), path.end());
    }

    #[test]
    fn ellipse_endpoints_lie_on_the_ellipse() {
        // Quarter arc of a 3x2 ellipse: (3,0) -> (0,2).
        let center = Point2::new(r(1), r(1));
        let span =
            elliptical_span(&center, &diag(3, 2), (&r(1), &r(0)), (&r(0), &r(1))).unwrap();
        let cps = span.control_points();
        // Start = center + (3,0), end = center + (0,2).
        assert_eq!(cps[0], &Point2::new(r(4), r(1)));
        assert_eq!(cps[2], &Point2::new(r(1), r(3)));
    }

    #[test]
    fn circle_map_matches_circular_decomposition_midpoint() {
        // With M = diag(r, r) the ellipse is a circle; the quarter-span point at
        // t = 1/2 sits on the circle of radius r about the center.
        let center = Point2::new(r(0), r(0));
        let span =
            elliptical_span(&center, &diag(5, 5), (&r(1), &r(0)), (&r(0), &r(1))).unwrap();
        let policy = CurvePolicy::certified();
        let half = (Real::from(1_i8) / Real::from(2_i8)).unwrap();
        let mid = match span.point_at(half, &policy) {
            crate::Classification::Decided(p) => p,
            other => panic!("expected decided midpoint, got {other:?}"),
        };
        // On the radius-5 circle: x^2 + y^2 == 25.
        let r2 = &(mid.x() * mid.x()) + &(mid.y() * mid.y());
        assert_eq!(
            crate::classify::compare_reals(&r2, &r(25), &policy),
            Some(std::cmp::Ordering::Equal)
        );
    }
}
