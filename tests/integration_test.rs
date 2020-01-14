use numlib::integrate::*;
use assert_approx_eq::assert_approx_eq;

fn square(x : f64) -> f64 {
    return x*x;
}

#[test]
fn trapezoid_quad() {
    let a : f64 =  0.66679;
    assert_approx_eq!(a, composite_trapezoid(-1.0, 1.0, 100, square), 2 as f64);
}

#[test]
fn simpson_test() {
    let a : f64 = 0.66679;
    assert_approx_eq!(a, simpson_quadrature(-1.0, 1.0, 100, square), 2 as f64);
}

#[test] 
fn three_eighth_simpson_test() { 
    let a: f64 = 0.66679;
    assert_approx_eq!(a, three_eighths_simpson(-1.0, 1.0, 100, square), 2 as f64);
}



