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

/* Adaptive Simpson's rule
 * Function: Computes an integral estimation by using the simspson rule and adapting to some user
 * specified tolerence.
 * Parameters:
 * f - the function
 * a - the inital endpoint
 * b - the final endpoint
 * eps - the tolerence
 * Returns:
 * The result of the integral
 */
pub fn adaptive_simpson(a: f64, b: f64, f: fn(f64) -> f64, eps: f64) -> f64 {
    let fa : f64 = f(a);
    let fb : f64 = f(b);
    let (m, fm, whole) = simpson_single_eval(f, a, b, fa, fb);
    return adaptive_simpson_quadrature(f, a, b, fa, fb, eps, m, fm, whole);
}

/* Adaptive_Simpson_Quadrature
 * Function: Performs the recursive nature of calculating adaptive simpsons rule to a tolerence
 * Parameters:
 * f - the function
 * a - the initial endpoint
 * b - the final endpoint
 * fa - the function evaluated at the initial endpoint
 * fb - the function evaluated at the final endpoint
 * eps - the tolerence
 * m - the midpoint
 * fm - the function evaluated at the midpoint
 * whole - the result of a single eval of simpsons quadrature
 * Returns: 
 * Either a recursive call to itself or the solution to the integral
 */
fn adaptive_simpson_quadrature(f: fn(f64) -> f64, a: f64, b: f64, fa: f64, fb: f64, eps: f64, m: f64, fm: f64, whole: f64) -> f64 {
    let (lm, flm, left) = simpson_single_eval(f, a, m, fa, fm);
    let (rm, frm, right) = simpson_single_eval(f, m, b, fm, fb);
    let delta : f64 = left+right-whole;
    if delta.abs() <= 15.0*eps {
        return left+right+delta/15.0;
    }
    return adaptive_simpson_quadrature(f, a, m, fa, fm, eps/2.0, lm, flm, left) + adaptive_simpson_quadrature(f, m, b, fm, fb, eps/2.0, rm, frm, right);
}

/* Simpson_Single_eval
 * Private Function: Computes simpsons rule at its midpoint and returns the value
 * Parameters:
 * f - the function
 * a - the inital endpoint
 * b - the final endpoint
 * fa - the function evaluated at the initial endpoint
 * fb - the function evaluated at the final endpoint
 * Returns: 
 * The result of a single eval of simpson's rule
 */
fn simpson_single_eval(f: fn(f64) -> f64, a: f64, b: f64, fa: f64, fb: f64) -> (f64, f64, f64) {
    let m : f64 = (a+b)/2.0;
    let fm : f64 = f(m);
    return (m, fm, (b-a).abs() / 6.0 * (fa + 4.0*fm  + fb));
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
            sum += 2.0*f(a+(i as f64)*h);
        } else{
            sum += 3.0*f(a+(i as f64)*h);
        }
    }
    return 3.0*h*sum/8.0;
}


