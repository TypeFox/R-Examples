### Function to calculate the theoretical raw moments (moments about 0)
### of a generalized inverse Gaussian distribution given its parameters.
gigRawMom <- function(order, Theta) {
    Theta <- as.numeric(Theta)
    parResult <- gigCheckPars(Theta)
    case <- parResult$case
    errMessage <- parResult$errMessage
    if (case == "error"){
        stop(errMessage)
    } else if (case =="normal"){
        lambda <- Theta[1]
        chi <- Theta[2]
        psi <- Theta[3]
        omega <- sqrt(chi*psi)
        eta <- sqrt(chi/psi)
        mom <- eta^order*besselRatio(omega, lambda, order)
    } else if (case == "gamma") {
        shape <- Theta[1]
        scale <- Theta[3]/2
        mom <- gammaRawMom(order, shape, scale)
    } else if (case == "invgamma"){
        shape <- -Theta[1]
        scale <- 2/Theta[2]
        mom <- gammaRawMom(-order, shape, scale)
    }
    return(mom)
}## End of gigRawMom()


gigMom <- function(order, Theta, about = 0) {
    if ((about != 0)&(!is.wholenumber(order))&(length(order)==1)){
        stop("Order must be a whole number except for moments about 0")
    }
    if ((order < 0)&(about != 0)) {
        stop("Order must be positive except for moments about 0")
    }
    Theta <- as.numeric(Theta)
    lambda <- Theta[1]
    chi <- Theta[2]
    psi <- Theta[3]
    omega <- sqrt(chi*psi)
    eta <- sqrt(chi/psi)
    if (about == 0){
        mom <- gigRawMom(order, Theta)
    } else {
        if (order == 0) {
            mom <- 1
        } else {
            rawMoments <- sapply(1:order, gigRawMom, Theta = Theta)
            mom <- momChangeAbout(order, rawMoments, 0, about)
        }
    }
    return(mom)
}## End of gigMom()


### Calculate moments of special cases
### chi = 0, gamma; psi = 0, inverse gamma
### Inverse gamma raw moments can be obtained from gamma raw moments
### An alternative would be to use mgamma and minvgamma from actuar
gammaRawMom <- function(order, shape = 1, rate = 1, scale = 1/rate) {
    if (order <= -shape) {
        mom <- Inf
    } else {
        mom <- (scale^order)*gamma(shape + order)/gamma(shape)
    }
    return(mom)
}## End of gammaRawMom()






