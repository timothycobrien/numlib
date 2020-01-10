// implementation of integrating techniques

/* Composite Trapezoidal Quadrature
 * Function: Computes an integral estimation using the composite trapezoidal rule
 * Parameters: 
 * f - the function,
 * n - the number of subintervals, 
 * a - the intial endpoint, 
 * b - the final endpoint
 * Returns:
 * The approximate solution to the integral
*/

pub fn composite_trapezoid(a: f32, b: f32, n: i32, f: fn(f32) -> f32) -> f32{
    let h : f32 = (b-a)/(n as f32);
    let mut sum : f32 = 0.0;
    for i in 1..n {
        sum += f(a+(i as f32)*h);
    }
    return h*(sum + f(a)/(2 as f32) + f(b)/(2 as f32));
}

