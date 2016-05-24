"parabolaSlopes" <-
function (x, y) 
{
    m <- length(x)
    m1 <- m - 1
    dx <- diff(x)
    dy <- diff(y)
    dydx = dy/dx
    if (m == 2) {
        yp <- rep(dydx, 2)
    } else {
        yp <- c((dydx[1]*(2*dx[1]+dx[2])-dydx[2]*dx[1])/(dx[1]+dx[2]),
        (dydx[-m1]*dx[-1] + dydx[-1]*dx[-m1])/(dx[-1] + dx[-m1]),
        (dydx[m1]*(2*dx[m1]+dx[m1-1])-dydx[m1-1]*dx[m1])/(dx[m1]+dx[m1-1]))
    }
    yp
}
