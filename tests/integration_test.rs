use numlib::integrate::*;
use assert_approx_eq::assert_approx_eq;

#[test]
fn trapezoid_quad() {
    fn square(x: f32) -> f32 {
        return x*x
    }
    let a : f32 =  0.66679;
    assert_approx_eq!(a as f32, composite_trapezoid(-1.0, 1.0, 100, square), 2f32);
}
