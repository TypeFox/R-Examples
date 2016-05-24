##
##  r k 4 . R  Classical Runge-Kutta
##


rk4 <- function(f, a, b, y0, n) {

    h <- (b-a)/n
    x <- seq(a+h, b, by = h)
    y <- numeric(n)

    k1 <- h * f(a, y0)
    k2 <- h * f(a + h / 2 , y0 + k1 / 2 ) 
    k3 <- h * f(a + h / 2 , y0 + k2 / 2 )
    k4 <- h * f(a + h ,     y0 + k3)
    y[1] <- y0 + k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6

    for (i in 1:(n-1)) {
       k1 <- h * f(x[i], y[i])
       k2 <- h * f(x[i] + h / 2, y[i] + k1 / 2) 
       k3 <- h * f(x[i] + h / 2 ,y[i] + k2 / 2 )
       k4 <- h * f(x[i] + h ,    y[i] + k3)
       y[i+1] = y[i] + k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6
    }

    x <- c( a, x)
    y <- c(y0, y)
    return(list(x = x, y = y))
}


rk4sys <- function(f, a, b, y0, n){


    m <- length(y0)
    h <- (b-a)/n
    x <- seq(a+h, b, by = h)
    y <- matrix(0, nrow = n, ncol = m)

    k1 <- h * f(a, y0)
    k2 <- h * f(a + h / 2 , y0 + k1 / 2 )
    k3 <- h * f(a + h / 2 , y0 + k2 / 2 )
    k4 <- h * f(a + h     , y0 + k3)
    y[1, ] <- y0 + k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6

    for (i in 1:(n-1)) {
        k1 <- h * f(x[i], y[i, ])
        k2 <- h * f(x[i] + h / 2, y[i, ] + k1 / 2 ) 
        k3 <- h * f(x[i] + h / 2, y[i, ] + k2 / 2)
        k4 <- h * f(x[i] + h    , y[i, ] + k3 )
        y[i+1, ] <- y[i, ] + k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6
    }

    x <- c(a, x)
    y <- rbind(y0, y)
    return(list(x = x, y = y))
}
