`Kclust` <-
function (r, sigma2, rho) 
{
    (pi * r^2) + ((1 - exp(-(r^2)/(4 * sigma2)))/rho)
}

