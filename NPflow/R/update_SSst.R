#' Return updated sufficient statistics S for skew t-distribution
#' with data matrix z
#'
#' For internal use only.
#'
#'@param z data matrix
#'
#'@param S previous sufficient statistics
#'
#'@param ltn random effects
#'
#'@param scale
#'
#'@param df skew t degrees of freedom
#'
#'@param hyperprior Default is \code{NULL}
#'
#'@keywords internal
#'
#'@export


update_SSst <- function(z, S, ltn, scale, df, hyperprior=NULL){

    b0_xi <- S[["b_xi"]]
    b0_psi <- S[["b_psi"]]
    nu0 <- S[["nu"]]
    lambda0 <- S[["lambda"]]
    B0 <- S[["B"]]
    if(is.null(B0)){
        D0_xi <- S[["D_xi"]]
        D0_psi <- S[["D_psi"]]
        B0inv <- diag(c(1/D0_xi, 1/D0_psi))
    }else{
        B0inv <- solve(B0)
        D0_xi <- B0[1,1]
        D0_psi <- B0[2,2]
    }


    sc_sr <- sqrt(scale)

    if(is.null(dim(z))){
        z <- matrix(z, ncol=1)
    }
    n <- ncol(z)

    X <- matrix(c(sc_sr, sc_sr*ltn), ncol=2, byrow=FALSE)
    B <- try(solve(crossprod(X) + B0inv))
    if(class(B)=="try-error"){
     #   browser()
        stop("error in inverting B")
    }
    temp <- apply(X=z, MARGIN=1, FUN=function(x){x*sc_sr})
    if(n<2){
        temp <- matrix(temp, ncol=nrow(z))
    }
    b <- (crossprod(temp,X) + cbind(b0_xi/D0_xi, b0_psi/D0_psi))%*%B


    b_xi <- b[,1]
    b_psi <- b[,2]

    nu1 <- nu0 + n

    eps2 <- tcrossprod(z[,1] - b_xi - ltn[1]*b_psi)*scale[1]

    #TODO optimize the way to calculate eps2 (Rcpp ?)
    if(n>1){
        for (i in 2:n){
            eps2 <- eps2 + tcrossprod(z[,i] - b_xi - ltn[i]*b_psi)*scale[i]
        }
    }


    #conjugate hyperprior on lambda: whishart distribution
    if(!is.null(hyperprior)){
        #g0 <- ncol(lambda0) + 5
        g0 <- nu0
        lambda0 <- wishrnd(n=round(nu0+g0), Sigma=solve(solve(lambda0)+solve(hyperprior[["Sigma"]])))
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
    S_up[["df"]] <- df


    return(S_up)
}