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

pub fn composite_trapezoid(a: f64, b: f64, n: i64, f: fn(f64) -> f64) -> f64{
    let h : f64 = (b-a)/(n as f64);
    let mut sum : f64 = 0.0;
    for i in 1..n {
        sum += f(a+(i as f64)*h);
    }
    return h*(sum + f(a)/(2 as f64) + f(b)/(2 as f64));
}

/* Simpsons Rule 
 * Function: Computes an integral estimation using simpsons rule. Preferable
 * to the trapezoidal rule in terms of error except on periodically smooth functions
 * Parameters: 
 * f - the function,
 * n - the number of subintervals,
 * a - the initial endpoint, 
 * b - the final endpoint
 * Returns:
 * The approximate solution to the integral
 *
 */

pub fn simpson_quadrature(a: f64, b: f64, n:i64, f: fn(f64) -> f64) -> f64{
    let h : f64 = (b-a)/(n as f64);
    let mut sum :f64 = 0.0;
    for i in 1..(n-1){
        sum += 2.0*f(a+(i as f64)*h) + 4.0*f(a+(i as f64)*h/(2 as f64));
    }
    return h*(sum + f(a) + f(b))/(6.0);
}
