#'Probability density function of multiple structured Normal inverse Wishart
#'
#'Probability density function of structured Normal inverse Wishart (sNiW)
#'for multiple inputs, on the log scale.
#'
#'@param U_xi a list of length n of observed mean vectors, each of dimension p
#'@param U_psi a list of length n of observed skew vectors of dimension p
#'@param U_Sigma a list of length n of observed covariance matrices, each of dimension p x p
#'@param U_xi0 a list of length K of mean vector parameters for sNiW, each of dimension p
#'@param U_psi0 a list of length K of mean vector parameters for sNiW, each of dimension p
#'@param U_B0 a list of length K of sturcturing matrix parameters for sNiW, each of dimension 2 x 2
#'@param U_Sigma0 a list of length K of covariance matrix parameters for sNiW, each of dimension p x p
#'@param U_df0 a list of length K of degrees of freedom parameters for sNiW, each of dimension p x p
#'@export
#'
#'@examples
#'hyperG0 <- list()
#'hyperG0$b_xi <- c(-1.6983129, -0.4819131)
#'hyperG0$b_psi <- c(-0.0641866, -0.7606068)
#'hyperG0$kappa <- 0.001
#'hyperG0$D_xi <- 16.951313
#'hyperG0$D_psi <- 1.255192
#'hyperG0$nu <- 27.67656
#'hyperG0$lambda <- matrix(c(2.3397761, -0.3975259,-0.3975259, 1.9601773), ncol=2)
#'
#'xi_list <- list()
#'psi_list <- list()
#'S_list <- list()
#'for(k in 1:1000){
#'  NNiW <- rNNiW(hyperG0, diagVar=FALSE)
#'  xi_list[[k]] <- NNiW[["xi"]]
#'  psi_list[[k]] <- NNiW[["psi"]]
#'  S_list[[k]] <- NNiW[["S"]]
#'}
#'mmsNiWlogpdf(U_xi=xi_list, U_psi=psi_list, U_Sigma=S_list,
#'             U_xi0=list(hyperG0$b_xi), U_psi0=list(hyperG0$b_psi) ,
#'             U_B0=list(diag(c(hyperG0$D_xi, hyperG0$D_psi))) ,
#'             U_Sigma0=list(hyperG0$lambda), U_df0=list(hyperG0$nu))
#'
#'
#'hyperG0 <- list()
#'hyperG0$b_xi <- c(-1.6983129)
#'hyperG0$b_psi <- c(-0.0641866)
#'hyperG0$kappa <- 0.001
#'hyperG0$D_xi <- 16.951313
#'hyperG0$D_psi <- 1.255192
#'hyperG0$nu <- 27.67656
#'hyperG0$lambda <- matrix(c(2.3397761), ncol=1)
#'#'xi_list <- list()
#'psi_list <- list()
#'S_list <- list()
#'for(k in 1:1000){
#'  NNiW <- rNNiW(hyperG0, diagVar=FALSE)
#'  xi_list[[k]] <- NNiW[["xi"]]
#'  psi_list[[k]] <- NNiW[["psi"]]
#'  S_list[[k]] <- NNiW[["S"]]
#'}
#'
#' mmsNiWlogpdf(U_xi=xi_list, U_psi=psi_list, U_Sigma=S_list,
#'             U_xi0=list(hyperG0$b_xi), U_psi0=list(hyperG0$b_psi) ,
#'             U_B0=list(diag(c(hyperG0$D_xi, hyperG0$D_psi))) ,
#'             U_Sigma0=list(hyperG0$lambda), U_df0=list(hyperG0$nu))
#'


mmsNiWlogpdf <- function(U_xi, U_psi, U_Sigma, U_xi0, U_psi0, U_B0, U_Sigma0, U_df0){

    loglik <- function(xi, psi, Sigma, xi0, psi0, B0, Lambda0, nu0){
        d <- length(xi)
        mu <- c(xi,psi)
        mu0 <- c(xi0, psi0)
        Sigmainv <- solve(Sigma)
        resl <- (-d/2*log(2*pi)
                 -(nu0 + d + 1)/2*log(det(Sigma))
                 -nu0*d/2*log(2)
                 +nu0/2*log(det(Lambda0))
                 -lgamma_mv(nu0/2, p=d)
                 -1/2*log(det(kronecker(B0, Sigma)))
                 -1/2*sum(diag(Lambda0%*%Sigmainv))
                 -1/2*t(mu-mu0)%*%kronecker(solve(B0), Sigmainv)%*%(mu-mu0)
        )
        #resl <- exp(resl)
    }
    ml <- function(xi, psi, Sigma, U_xi0, U_psi0, U_B0, U_Sigma0, U_df0){
        resml <- mapply(FUN = loglik,
                        xi0 = U_xi0, psi0 = U_psi0, B0 = U_B0,
                        Lambda0 = U_Sigma0, nu0 = U_df0,
                        MoreArgs=list(xi, psi, Sigma))
#         constnorm <- sum(resml)
#         if(constnorm>0){
#             resml <- resml/constnorm
#         }
        return(resml)
    }

    mml <- function(U_xi, U_psi, U_Sigma, U_xi0, U_psi0, U_B0, U_Sigma0, U_df0){
        resmml <- mapply(FUN = ml,
                        xi = U_xi, psi = U_psi, Sigma = U_Sigma,
                        MoreArgs=list(U_xi0, U_psi0, U_B0, U_Sigma0, U_df0))
    }

    res <- mml(U_xi, U_psi, U_Sigma, U_xi0, U_psi0, U_B0, U_Sigma0, U_df0)
    return(res)

}

msNiWlogpdf <- function(xi, psi, Sigma, U_xi0, U_psi0, U_B0, U_Sigma0, U_df0){

    loglik <- function(xi, psi, Sigma, xi0, psi0, B0, Lambda0, nu0){
        d <- length(xi)
        mu <- c(xi,psi)
        mu0 <- c(xi0, psi0)
        Sigmainv <- solve(Sigma)
        resl <- (-d/2*log(2*pi)
                 -(nu0 + d + 1)/2*log(det(Sigma))
                 -nu0*d/2*log(2)
                 +nu0/2*log(det(Lambda0))
                 -lgamma_mv(nu0/2, p=d)
                 -1/2*log(det(kronecker(B0, Sigma)))
                 -1/2*sum(diag(Lambda0%*%Sigmainv))
                 -1/2*t(mu-mu0)%*%kronecker(solve(B0), Sigmainv)%*%(mu-mu0)
        )
        #resl <- exp(resl)
    }


        res <- mapply(FUN = loglik,
                        xi0 = U_xi0, psi0 = U_psi0, B0 = U_B0,
                        Lambda0 = U_Sigma0, nu0 = U_df0,
                        MoreArgs=list(xi, psi, Sigma))

    return(res)

}

sNiWlogpdf <- function(xi, psi, Sigma, U_xi0, U_psi0, U_B0, U_Sigma0, U_df0){

    loglik <- function(xi, psi, Sigma, xi0, psi0, B0, Lambda0, nu0){
        d <- length(xi)
        mu <- c(xi,psi)
        mu0 <- c(xi0, psi0)
        Sigmainv <- solve(Sigma)
        resl <- (-d/2*log(2*pi)
                 -(nu0 + d + 1)/2*log(det(Sigma))
                 -nu0*d/2*log(2)
                 +nu0/2*log(det(Lambda0))
                 -lgamma_mv(nu0/2, p=d)
                 -1/2*log(det(kronecker(B0, Sigma)))
                 -1/2*sum(diag(Lambda0%*%Sigmainv))
                 -1/2*t(mu-mu0)%*%kronecker(solve(B0), Sigmainv)%*%(mu-mu0)
        )
        #resl <- exp(resl)
    }

    res <- loglik(xi, psi, Sigma, xi0 = U_xi0, psi0 = U_psi0, B0 = U_B0,
                  Lambda0 = U_Sigma0, nu0 = U_df0
    )

    return(res)

}

gamma_mv <- function(x,p){
    pi^(p*(p-1)/4)*prod(gamma(x+(1-1:p)/2))
}

#'Multivariate log gamma function
#'@param x strictly positive real number
#'@param p integer
#'@export
lgamma_mv <- function(x,p){
    (p*(p-1)/4)*log(pi) + sum(lgamma(x+(1-1:p)/2))
}

digamma_mv <- function(x,p){
    sum(digamma(x+(1-1:p)/2))
}