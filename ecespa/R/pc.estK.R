`pc.estK` <-
function (Kobs, r, sigma2 = NULL, rho = NULL) 
{
    theta = c(sigma2, rho)
    if (is.null(theta)) 
        theta = c(1, 1)
    D.theta = function(theta, Kobs, r) {
        sum((Kobs^(0.25) - Kclust(r, theta[1], theta[2])^(0.25))^2)
    }
    nsfit = optim(par = theta, fn = D.theta, Kobs = Kobs, r = r)
    return(list(sigma2 = nsfit$par[1], rho = nsfit$par[2]))
}

