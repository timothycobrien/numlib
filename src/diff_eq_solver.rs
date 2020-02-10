// implentation of solving ODEs


/* Euler's (RK1)
 * Function: Performs a runge-kutta interation solve of first order for approximating diff eqs
 * Parameters:
 * deriv - ther derivative of the function
 * x0 - the inital function point
 * n - the number of steps
 * x - the point to approximate at
 * y - the initial y value
 * Returns: 
 * The approximate solution to the ode at x
 */

pub fn euler(mut x0: f64, x: f64, mut y: f64, n: i64, deriv: fn(f64, f64) -> f64) -> f64 {
    let h : f64 = (x-x0)/(n as f64).abs();
    for _ in 1..n {
        y += h*deriv(x0, y);
        x0 = x0 + h;
    }
    return y;
}

/* Runge-Kutta 2 
 * Function: Performs a runge-kutta interation solver of order 2 for approximating diff eqs
 * Parameters:
 * deriv - the derivative of the function
 * x0 - the inital function point
 * n - the number of steps
 * x - the point to approximate at
 * y - the initial y value
 * Returns: 
 * The approximate solution to the ODE at x
 */

pub fn runge_kutta_2(mut x0: f64, x: f64, mut y: f64, n: i64, deriv: fn(f64, f64) -> f64) -> f64{
    let h : f64 = (x-x0)/(n as f64).abs();
    let mut k1 : f64;
    let mut k2 : f64;
    for _ in 1..n {
        k1 = h*deriv(x0, y);
        k2 = h*deriv(x0 + 0.5*h, y+0.5*k1);
        y = y + k2;
        x0 = x0+h;
    }
    return y;
}

/* Runge-kutta 4
 * Function: Performs a runge-kutta iteration solver of order 4 for approximating diff eqs
 * Parameters:
 * deriv - the derivative of the function
 * x0 - the inital function point
 * n - the number of stpes
 * x - the point ot approximate at
 * y - the initial y value
 * Returns:
 * The approximate solution to the ODE at x
 */

pub fn runge_kutta_4(mut x0: f64, x: f64, mut y: f64, n: i64, deriv: fn(f64, f64) -> f64) -> f64 {
    let h: f64 = (x-x0)/(n as f64).abs();
    let mut k1 : f64;
    let mut k2 : f64;
    let mut k3 : f64;
    let mut k4 : f64;
    for _ in 1..n {
        k1 = h*deriv(x0, y);
        k2 = h*deriv(x0+0.5*h, y+0.5*k1);
        k3 = h*deriv(x0+0.5*h, y+0.5*k2);
        k4 = h*deriv(x0+h, y+k3);
        y = y + (1.0/6.0)*(k1+k4)+(1.0/3.0)*(k2+k3);
        x0 = x0 + h;
    }
    return y;
}

/* Adams-Brashforth
 * Function: Uses the two-step Adams-Brashforth method by using euler to calculate the first point
 * and propogating Adams-Brashforth from there
 * Paramters:
 * deriv - the derivative of the function
 * x0 - the initial function point
 * n - the number of step
 * x - the point to approximate at
 * y0 - the initial y value
 * Returns:
 * The approximate solution to the ODE at x
 */
pub fn adams_brashforth(x0: f64, x: f64, mut y0: f64, n: i64, deriv: fn(f64, f64)->f64) -> f64{
    let h: f64 = (x-x0)/(n as f64).abs();
    // Euler approximation for first step
    let mut y1: f64 = y0 + h*deriv(x0, y0);
    let mut x1: f64 = x0 + h;
    for _ in 1..(n-1){
        let tmp: f64 = y1;
        y1 += 1.5*h*deriv(x1, y1) - 0.5*h*deriv(x0, y0);
        y0 = tmp;
        x1 += h;
        x1 += h;
    }
    return y1;
}

/* Adams-Moulton
 * Function: Uses the two steps Adams-Moulton method by using euler to calculate the first point
 * point
 * Paramters:
 * deriv - the derivative of the function
 * x0 = the initial function point
 * n - the number of steps
 * x - the point to approximate at 
 * y0 - the inital y value
 * Returns:
 * The approximate solution to the ODE at x
 */
pub fn adams_moulton(mut x0: f64, x: f64, mut y0: f64, n: i64, deriv: fn(f64, f64) -> f64) -> f64 {
    let h: f64 = (x-x0)/(n as f64).abs();
    // Euler Approx
    let mut y1: f64 = y0 + h*deriv(x0, y0);
    let mut x1: f64 = x0+h;
    let mut y2: f64 = y1 + h*deriv(x1, y1);
    let mut x2: f64 = x1+h;
    for _ in 1..(n-2){
        let tmp: f64 = y2;
        y2 = y1 + h*((5.0/12.0)*deriv(x2, y2) + (2.0/3.0)*deriv(x1, y1) - (1.0/12.0)*deriv(x0, y0));
        y0 = y1;
        x0 = x1;
        y1 = tmp;
        x1 = x2;
        x2 += h;
    }
    return y2;
}

        








