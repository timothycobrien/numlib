// implentation of solving ODEs

/* Runge-Kutta 2 
 * Function: Performs a runge-kutta interation solver of order 2 for approximating diff eqs
 * Parameters:
 * deriv - the derivative of the function
 * x0 - the inital function point
 * n - the number of steps
 * x - the point to approximate at
 * Returns: 
 * The approximate solution to the ODE at x
 */

pub fn runge_kutta_2(x0: f64, x: f64, n: i64, deriv: fn(f64, f64) -> f64) -> f64{
    let h : f64 = (x-x0)/(n as f64).abs()
    let mut y = f64 = 0.0;
    for i in 1..n {
        let mut k1 = h*deriv(x0);
        let mut k2 = h*deriv(x0 + 0.5*h, y+0.5*k1);
        y = y + k2;
        x = x+h
    }
    return y;
}


