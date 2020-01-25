// implentation of solving ODEs

/* Runge-Kutta 2 
 * Function: Performs a runge-kutta interation solver of order 2 for approximating diff eqs
 * Parameters:
 * deriv - the derivative of the function
 * x0 - the inital function point
 * n - the number of steps
 * x - the point to approximate at
 * y - the inital y value
 * Returns: 
 * The approximate solution to the ODE at x
 */

pub fn runge_kutta_2(x0: f64, x: f64, y: f64, n: i64, deriv: fn(f64, f64) -> f64) -> f64{
    let h : f64 = (x-x0)/(n as f64).abs();
    let k1 : f64 = 0;
    let k2 : f64 = 0;
    for i in 1..n {
        k1 = h*deriv(x0, y);
        k2 = h*deriv(x0 + 0.5*h, y+0.5*k1);
        y = y + k2;
        x = x+h;
    }
    return y;
}

/* Runge-kutta 4
 * Function: Performs a runge-kutta interation solve of order 2 for approximating diff eqs
 * Parameters:
 * deriv - the derivative of the function
 * x0 - the inital function point
 * n - the number of stpes
 * x - the point ot approximate at
 * y - the initial y value
 * Returns:
 * The approximate solution to the ODE at x
 */

pub fn runge_kutta_4(x0: f64, x: f64, y: f64, n: i64, deriv: fn(f64, f64) -> f64) -> f64 {
    let h: f64 = (x-x0)/(n as f64).abs();
    let k1 : f64 = 0;
    let k2 : f64 = 0;
    let k3 : f64 = 0;
    let k4 : f64 = 0;
    for i in 1..n {
        k1 = h*deriv(x0, y);
        k2 = h*deriv(x0+0.5*h, y+0.5*k1);
        k3 = h*deriv(x0+0.5*h, y+0.5*k2);
        k4 = h*deriv(x0+h, y+k3);
        y = y + (1.0/6.0)*(k1+k4)+(1.0/3.0)*(k2+k3);
        x = x + h;
    }
    return y;
}


        
