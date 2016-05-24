soilwater5 <- 
function(x, theta_R, theta_S, alpha, n, m = 1 - 1/n, b0, b1, b2) 
{
    sat.index <- (1 + (alpha * x)^n)^(-m)
    out <- theta_R + (theta_S - theta_R) * sat.index + b0 + b1*x + b2*x^2
    return(out)
}
