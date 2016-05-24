AT.inverse.gamma.response <- function(surv, alpha, beta){
    return(-0.5*(alpha/beta)+sqrt((alpha/(2*beta))^2 - 1/beta * log(surv)))
}
