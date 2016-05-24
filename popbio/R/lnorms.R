lnorms <- function(n, mean=2,var=1)
{
    nmeans  <- log(mean) - 0.5 * log(var/mean^2 + 1)
    nvars   <- log(var/mean^2 + 1)
    normals <- rnorm(n) * sqrt(nvars) + nmeans
    lns     <- exp(normals)
    lns   
} 

