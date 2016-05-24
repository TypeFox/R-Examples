CIdelta.fun <-
function (VARbeta, lambdafit, covariates, clevel = 0.95) 
{
    covariates <<- covariates
    VARbeta <<- VARbeta
    n <- dim(covariates)[1]
    Varlambda <- rep(NA, n)
    calcVlambda.fun <- function(i) {
        Varlambda[i] <<- matrix(covariates[i, ], nrow = 1) %*% 
            VARbeta %*% matrix(covariates[i, ], ncol = 1)
    }
    aux <- sapply(c(1:n), FUN = calcVlambda.fun)
    Varlambda <- Varlambda * lambdafit^2
    UIdlambda <- lambdafit + qnorm(1 - (1 - clevel)/2) * Varlambda^0.5
    LIdlambda <- pmax(0,lambdafit - qnorm(1 - (1 - clevel)/2) * Varlambda^0.5)
    return(list(UIlambda = UIdlambda, LIlambda = LIdlambda, lambdafit = lambdafit))
}
