use numlib::integrate::*;

fn main() {
    println!("Demonstration of the trapezoid quadrature and simpson quadrature of x^2 over [-1, 1]. The exact solution is 2/3.");
    println!("Trapezoid solution {}", composite_trapezoid(-1.0, 1.0, 100, square));
    println!("Simpson's solution {}", simpson_quadrature(-1.0, 1.0, 100, square)); 
    println!("3/8's Simpson's solution {}", three_eighths_simpson(-1.0, 1.0, 100, square));
}

fn square(x : f64) -> f64{
    return x*x;
}
