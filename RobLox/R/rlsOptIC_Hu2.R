###############################################################################
## chi function
###############################################################################
lsHu2.chi <- function(x, c0){
    beta.c0 <- 2*pnorm(c0) - 1 - 2*c0*dnorm(c0) + 2*c0^2*pnorm(-c0)
    return(pmin(x^2, c0^2) - beta.c0) 
}

###############################################################################
## Computation of bias
###############################################################################
.Hu2rlsGetbias <- function(x, k, c0){
    beta.c0 <- 2*pnorm(c0) - 1 - 2*c0*dnorm(c0) + 2*c0^2*pnorm(-c0)
    sqrt(pmin(k^2, x^2)/(2*pnorm(k)-1)^2 + (pmin(c0^2, x^2) - beta.c0)^2/
        (2*(2*pnorm(c0) - 1) - 4*c0*dnorm(c0))^2)
}


###############################################################################
## computation of asymptotic variance
###############################################################################
.Hu2rlsGetvar <- function(k, c0){
    beta.c0 <- 2*pnorm(c0) - 1 - 2*c0*dnorm(c0) + 2*c0^2*pnorm(-c0)
    Var.loc <- (2*pnorm(k) - 1 - 2*k*dnorm(k) + 2*k^2*(1-pnorm(k)))/(2*pnorm(k)-1)^2

    hilf <- 3*(2*pnorm(c0)-1) - 2*c0^3*dnorm(c0) - 6*c0*dnorm(c0) + 2*c0^4*pnorm(-c0)
    Var.sc <- (hilf - beta.c0^2)/(2*(2*pnorm(c0) - 1) - 4*c0*dnorm(c0))^2

    return(Var.loc + Var.sc)
}

###############################################################################
## computation of maximum asymptotic MSE
###############################################################################
.Hu2rlsGetmse <- function(kc0, r, MAX){
    k <- kc0[1]; c0 <- kc0[2]

    # constraints for k, c0
    if(k < 0 || c0 < 0) return(MAX)

    beta.c0 <- 2*pnorm(c0) - 1 - 2*c0*dnorm(c0) + 2*c0^2*pnorm(-c0)
    A.loc <- 1/(2*pnorm(k)-1)
    A.sc <- 1/(2*(2*pnorm(c0) - 1) - 4*c0*dnorm(c0))

    bias <- max(.Hu2rlsGetbias(x = 0, k = k, c0 = c0), 
                .Hu2rlsGetbias(x = k, k = k, c0 = c0), 
                .Hu2rlsGetbias(x = c0, k = k, c0 = c0), 
                .Hu2rlsGetbias(x = sqrt(beta.c0), k = k, c0 = c0), 
                .Hu2rlsGetbias(x = sqrt(max(0, beta.c0-0.5*A.loc^2/A.sc^2)), 
                                k = k, c0 = c0))

    Var <- .Hu2rlsGetvar(k = k, c0 = c0)

    return(Var + r^2*bias^2)
}


###############################################################################
## optimal IC
###############################################################################
rlsOptIC.Hu2 <- function(r, k.start = 1.5, c.start = 1.5, delta = 1e-6, MAX = 100){
    res <- optim(c(k.start, c.start), .Hu2rlsGetmse, method = "Nelder-Mead", 
                control = list(reltol=delta), r = r, MAX = MAX)

    k <- res$par[1]; c0 <- res$par[2]

    beta.c0 <- 2*pnorm(c0) - 1 - 2*c0*dnorm(c0) + 2*c0^2*pnorm(-c0)
    A.loc <- 1/(2*pnorm(k)-1)
    A.sc <- 1/(2*(2*pnorm(c0) - 1) - 4*c0*dnorm(c0))

    bias <- max(.Hu2rlsGetbias(x = 0, k = k, c0 = c0), 
                .Hu2rlsGetbias(x = k, k = k, c0 = c0), 
                .Hu2rlsGetbias(x = c0, k = k, c0 = c0), 
                .Hu2rlsGetbias(x = sqrt(beta.c0), k = k, c0 = c0), 
                .Hu2rlsGetbias(x = sqrt(max(0, beta.c0-0.5*A.loc^2/A.sc^2)), 
                                k = k, c0 = c0))

    fct1 <- function(x){ A.loc*sign(x)*pmin(abs(x), k) }
    body(fct1) <- substitute({ A.loc*sign(x)*pmin(abs(x), k) },
                        list(k = k, A.loc = A.loc))
    fct2 <- function(x){ A.sc*(pmin(x^2, c0^2) - beta.c0) }
    body(fct2) <- substitute({ A.sc*(pmin(x^2, c0^2) - beta.c0) },
                        list(c0 = c0, beta.c0 = beta.c0, A = A.sc))

    return(IC(name = "IC of Hu2 type", 
              Curve = EuclRandVarList(RealRandVariable(Map = list(fct1, fct2), Domain = Reals())),
              Risks = list(asMSE = res$value, asBias = bias, asCov = res$value - r^2*bias^2), 
              Infos = matrix(c("rlsOptIC.Hu2", "optimally robust IC for Hu2 estimators and 'asMSE'",
                               "rlsOptIC.Hu2", paste("where k =", round(k, 3), "and c =", round(c0, 3))), 
                           ncol=2, byrow = TRUE, dimnames=list(character(0), c("method", "message"))), 
              CallL2Fam = call("NormLocationScaleFamily")))
}
