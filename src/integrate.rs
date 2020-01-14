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
    let h : f64 = (b-a)/(n as f64).abs();
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
    let h : f64 = (b-a)/(n as f64).abs();
    let mut sum :f64 = 0.0;
    for i in 1..n{
        if i % 2 == 1 {
            sum += 4.0 * f(a+(i as f64)*h);
        } else {
            sum += 2.0 * f(a+(i as f64)*h);
        }
    }
    return h*(sum + f(a) + f(b))/(3.0);
}

/* 3/8 Simpson's Rule
 * Function: Computes an integral estimation using the 3/8s simpsons rule. About twice as accurate
 * as simpson's original rule. 
 * Parameters:
 * f - the function, 
 * n - the number of subintervals,
 * a - the initial endpoint, 
 * b - the final endpoint
 * Returns: The approximate solution to the integral. 
 *
 */

pub fn three_eighths_simpson(a: f64, b: f64, n: i64, f: fn(f64) -> f64) -> f64{
    let h : f64 = (b-a)/(n as f64).abs();
    let mut sum : f64 = 0.0;
    for i in 1..n{
        if i % 3 == 0 {
            sum += 2*f(a+(i as f64)*h);
        } else{
            sum += 3*f(a+(i as f64)*h);
        }
    }
    return 3.0*h*sum/8.0;
}





}
