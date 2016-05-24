MWI<-function (itemBank, item, x, it.given, model = NULL, lower = -4, 
    upper = 4, nqp = 33, type = "MLWI", priorDist = "norm", priorPar = c(0, 
        1), D = 1) 
{
    if (type != "MLWI" & type != "MPWI") 
        stop("'type' must be either 'MLWI' or 'MPWI'", call. = FALSE)
    if (is.null(model)) {
        L <- function(th, x, par) prod(Pi(th, par, D = D)$Pi^x * 
            (1 - Pi(th, par, D = D)$Pi)^(1 - x))
        X <- seq(from = lower, to = upper, length = nqp)
        lik <- sapply(X, L, x, it.given)
        Iprov <- function(t) Ii(t, itemBank, D = D)$Ii[item]
        info <- sapply(X, Iprov)
        crit.value <- lik * info
        if (type == "MPWI") {
            pd <- NULL
            for (k in 1:length(X)) pd[k] <- switch(priorDist, 
                norm = dnorm(X[k], priorPar[1], priorPar[2]), 
                unif = dunif(X[k], priorPar[1], priorPar[2]))
            crit.value <- crit.value * pd
        }
    }
    else {
        LL <- function(th, it.given, x, model, D = 1) {
if (dim(it.given)[1]==0) res<-1
else{
            prob <- Pi(th, it.given, model = model, D = D)$Pi
            res <- 1
            for (i in 1:length(x)) res <- res * prob[i, x[i] + 
                1]
}
            return(res)
        }
        X <- seq(from = lower, to = upper, length = nqp)
        lik <- sapply(X, LL, it.given, x, model, D = D)
        Iprov <- function(t) sum(Ii(t, itemBank[item, ], model = model, 
            D = D)$Ii)
        info <- sapply(X, Iprov)
        crit.value <- lik * info
        if (type == "MPWI") {
            pd <- NULL
            for (k in 1:length(X)) pd[k] <- switch(priorDist, 
                norm = dnorm(X[k], priorPar[1], priorPar[2]), 
                unif = dunif(X[k], priorPar[1], priorPar[2]))
            crit.value <- crit.value * pd
        }
    }
    RES <- integrate.catR(X, crit.value)
    return(RES)
}
