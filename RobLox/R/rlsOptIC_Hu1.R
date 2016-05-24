###############################################################################
## optimal psi function
###############################################################################
lsHu1.chi <- function(x, k){ 
    beta.k <- 2*pnorm(k) - 1 - 2*k*dnorm(k) + 2*k^2*pnorm(-k)
    return(pmin(x^2, k^2) - beta.k) 
}

###############################################################################
## computation of bias
###############################################################################
.Hu1rlsGetbias <- function(x, k){
    beta.k <- 2*pnorm(k) - 1 - 2*k*dnorm(k) + 2*k^2*pnorm(-k)
    sqrt(pmin(k^2, x^2)/(2*pnorm(k)-1)^2 + (pmin(k^2, x^2) - beta.k)^2/
        (2*(2*pnorm(k) - 1) - 4*k*dnorm(k))^2)
}


###############################################################################
## computation of asymptotic variance
###############################################################################
.Hu1rlsGetvar <- function(k){
    beta.k <- 2*pnorm(k) - 1 - 2*k*dnorm(k) + 2*k^2*pnorm(-k)
    E.psi.4 <- 3*(2*pnorm(k)-1) - 2*(k^3+3*k)*dnorm(k) + 2*k^4*pnorm(-k)

    Var.loc <- beta.k/(2*pnorm(k)-1)^2
    Var.sc <- (E.psi.4 - beta.k^2)/(2*(2*pnorm(k) - 1) - 4*k*dnorm(k))^2

    return(Var.loc+Var.sc)
}


###############################################################################
## computation of maximum asymptotic MSE
###############################################################################
.Hu1rlsGetmse <- function(k, r){
    beta.k <- 2*pnorm(k) - 1 - 2*k*dnorm(k) + 2*k^2*pnorm(-k)
    A.loc <- 1/(2*pnorm(k)-1)
    A.sc <- 1/(2*(2*pnorm(k) - 1) - 4*k*dnorm(k))

    bias <- max(.Hu1rlsGetbias(x = 0, k = k), 
                .Hu1rlsGetbias(x = k, k = k), 
                .Hu1rlsGetbias(x = sqrt(beta.k), k = k), 
                .Hu1rlsGetbias(x = sqrt(max(0,beta.k-0.5*A.loc^2/A.sc^2)), k = k))

    return(.Hu1rlsGetvar(k)+r^2*bias^2)
}


###############################################################################
## optimal IC
###############################################################################
rlsOptIC.Hu1 <- function(r, kUp = 2.5, delta = 1e-6){
    res <- optimize(f = .Hu1rlsGetmse, lower = 1e-4, upper = kUp, 
                tol = delta, r = r)

    k <- res$minimum
    beta.k <- 2*pnorm(k) - 1 - 2*k*dnorm(k) + 2*k^2*pnorm(-k)
    A.loc <- 1/(2*pnorm(k)-1)
    A.sc <- 1/(2*(2*pnorm(k) - 1) - 4*k*dnorm(k))

    bias <- max(.Hu1rlsGetbias(x = 0, k = k), 
                .Hu1rlsGetbias(x = k, k = k), 
                .Hu1rlsGetbias(x = sqrt(beta.k), k = k), 
                .Hu1rlsGetbias(x = sqrt(max(0,beta.k-0.5*A.loc^2/A.sc^2)), k = k))

    fct1 <- function(x){ A.loc*sign(x)*pmin(abs(x), k) }
    body(fct1) <- substitute({ A.loc*sign(x)*pmin(abs(x), k) },
                        list(k = k, A.loc = A.loc))
    fct2 <- function(x){ A.sc*(pmin(x^2, k^2) - beta.k) }
    body(fct2) <- substitute({ A.sc*(pmin(x^2, k^2) - beta.k) },
                        list(k = k, beta.k = beta.k, A.sc = A.sc))

    return(IC(name = "IC of Hu1 type", 
              Curve = EuclRandVarList(RealRandVariable(Map = list(fct1, fct2), Domain = Reals())),
              Risks = list(asMSE = res$objective, asBias = bias, asCov = res$objective - r^2*bias^2), 
              Infos = matrix(c("rlsOptIC.Hu1", "optimally robust IC for Hu1 estimators and 'asMSE'",
                               "rlsOptIC.Hu1", paste("where k =", round(k, 3))), 
                           ncol=2, byrow = TRUE, dimnames=list(character(0), c("method", "message"))), 
              CallL2Fam = call("NormLocationScaleFamily")))
}
