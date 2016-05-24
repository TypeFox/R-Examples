
# cubic hermite polynomial interpolation
interp.h3 <- function( x, x0, x1, y0, y1, yp0, yp1 ) {
    delx <- x1 - x0
    s    <- ( x - x0 )/delx
    s2   <- s*s
    s3   <- s*s2
    dydx <- ( y1 - y0 )/delx
    
    c    <- 3*dydx - 2*yp0 - yp1
    d    <- yp0 + yp1 - 2*dydx
    
    y    <- y0  + (delx*d*s3 + delx*c*s2 + delx*yp0*s)
    #yp   <- yp0 + (3*d*s2 + 2*c*s) 
    
    y
}
