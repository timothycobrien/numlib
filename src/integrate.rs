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

pub fn composite_trapezoid(a: f32, b: f32, n: i32, f: fn(f32) -> f32){
    let h = (b-a)/n;
    let sum = 0;
    for i in 0..n {
        sum += f(a+i*h);
    }
    return sum + f(a)/2 + f(b)/2;
}
