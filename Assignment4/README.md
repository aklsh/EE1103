# Integration

Integrate the same dataset used for interpolation and compare the trapezoid rule against the more advanced Romberg Quadrature (C code is in the Implementation here: [Romberg's Method](https://en.wikipedia.org/wiki/Romberg%27s_method)).

##### Use this for practice - Numerically integrate these functions:

    f(x) = exp(-x) * pow(cos(x),2) between 0 and pi
    g(x) = cos(2*acos(x)) between -1 and 1

Use the three different methods (midpoint, trapezoid and Simpson's) and compare the convergence of your integral for different step sizes.
