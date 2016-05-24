mod8p <-
function(theta, x) 
{
    A <- theta[1]
    B <- theta[2]
    C <- theta[3]
    D <- theta[4]
    E <- theta[5]
    F <- theta[6]
    G <- theta[7]
    H <- theta[8]
    f.x <- A^((x + B)^C) + D * exp(-E * (log(x) - log(F))^2) + 
        (G * (H^x))/(1 + G * (H^x))
    f.x2 <- f.x
    f.x2[1e-06 > f.x] <- 1e-06
    f.x2[0.999999 < f.x] <- 0.999999
    return(f.x2)
}

