#' Return updated sufficient statistics S with new data matrix z
#'
#' For internal use only.
#'
#'@param z data matrix
#'
#'@param S previous sufficient statistics
#'
#'@param ltn random effects
#'
#'@param hyperprior Default is \code{NULL}
#'
#'@keywords internal
#'
#'@export


update_SSsn <- function(z, S, ltn, hyperprior=NULL){

    b0_xi <- S[["b_xi"]]
    b0_psi <- S[["b_psi"]]
    D0_xi <- S[["D_xi"]]
    D0_psi <- S[["D_psi"]]
    nu0 <- S[["nu"]]
    lambda0 <- S[["lambda"]]

    if(is.null(dim(z))){
        z <- matrix(z, ncol=1)
    }
    n <- ncol(z)

    X <- matrix(c(rep(1, n), ltn), ncol=2, byrow=FALSE)
    B <- solve(crossprod(X)+diag(c(1/D0_xi, 1/D0_psi)))
    b <- (z%*%X + cbind(b0_xi/D0_xi, b0_psi/D0_psi))%*%B


    b_xi <- b[,1]
    b_psi <- b[,2]

    nu1 <- nu0 + n #c


    eps2 <- tcrossprod(z[,1] - b_xi - ltn[1]*b_psi)

    if(n>1){
        for (i in 2:n){
            eps2 <- eps2 + tcrossprod(z[,i] - b_xi - ltn[i]*b_psi)
        }
    }

    #conjugate hyperprior on lambda: whishart distribution
    if(!is.null(hyperprior)){
        #g0 <- ncol(lambda0) + 5
        sigma_inv <- try(solve(hyperprior[["Sigma"]]), silent = TRUE)
        if(!inherits(sigma_inv, "try-error")){
         g0 <- nu0
         lambda0 <- wishrnd(n=nu0+g0, Sigma=solve(solve(lambda0) + sigma_inv))
        }
    }

    lambda1 <- lambda0 + (eps2 + tcrossprod(b_xi-b0_xi)/D0_xi + tcrossprod(b_psi-b0_psi)/D0_psi)


    S_up <- list()
    S_up[["b_xi"]] <- b_xi
    S_up[["b_psi"]] <- b_psi
    S_up[["B"]] <- B
    S_up[["nu"]] <- nu1
    S_up[["lambda"]] <- lambda1
    S_up[["D_xi"]] <- D0_xi
    S_up[["D_psi"]] <- D0_psi


    return(S_up)
}