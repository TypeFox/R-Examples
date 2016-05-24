lotkaVolterra <- function(t, y, parameters){
    x <- y[1]
    y <- y[2]
    lambda  <- parameters[1]
    epsilon <- parameters[2]
    eta     <- parameters[3]
    delta   <- parameters[4]
    dy    <- numeric(2)
    dy[1] <- lambda*x - epsilon*x*y
    dy[2] <- eta*x*y - delta*y
    list(dy)
}