ipc.estK2 <-
function (mippp, lambda = NULL, correction = "iso", r = NULL, 
    sigma2 = NULL, rho = NULL, q = 1/4, p = 2,nlarge=NULL,...) 
{
    # Define theoretical K function for a PCP
    Kclust <- function(r, sigma2, rho) {
        (pi * r^2) + ((1 - exp(-(r^2)/(4 * sigma2)))/rho)
    }
    # Define minimum contrats function, i.e., D(theta) function of Diggle
    D.theta <- function(theta, Kobs, r) {
        sum((Kobs^q - Kclust(r, theta[1], theta[2])^q)^p)
    }
    lambdaname <- deparse(substitute(lambda))
    if (is.null(lambda)) {
        #lambda <- predict(ppm(mippp), type = "trend")
lambda<- rep( mippp$n/area.owin(mippp$window),mippp$n)
        lambdaname <- NULL
    }
    Kobs <- Kinhom(mippp, r = r, correction = correction, lambda = lambda, nlarge=nlarge)
    if (is.null(r)) 
        r <- Kobs$r
    Kobs <- Kobs[[3]]
    theta <- c(sigma2, rho)
    if (is.null(theta)) 
        theta <- c(1, 1)
    nsfit <- optim(par = theta, fn = D.theta, Kobs = Kobs, r = r,...) 
    Kfit <- Kclust(r, nsfit$par[1], nsfit$par[2])
    dataname <- deparse(substitute(mippp))
    dtheta <- sum((Kobs^q - Kfit^q)^p)
    result <- list(sigma2 = nsfit$par[1], rho = nsfit$par[2], 
        d.theta = dtheta, Kobs = Kobs, Kfit = Kfit, r = r, data = mippp, 
        lambda = lambda, dataname = dataname, lambdaname = lambdaname, 
        p = p, q = q)
    class(result) <- c("ecespa.minconfit", class(result))
    return(result)
}
