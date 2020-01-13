use numlib::integrate::*;
use std::f64::consts::PI;

fn main() {
    println!("Demonstration of composite trapezoid rule and simpson rule over periodic intervals. We integrate sin(x) over 0 to 4pi. The solution should be 0.");
    println!("Trapezoid solution: {}", composite_trapezoid(0.0, PI*4.0, 100, sin));
    println!("Simpson's solution: {}", simpson_quadrature(0.0, PI*4.0, 100, sin));

}

fn sin(x : f64) -> f64 {
    return x.sin();
}
