###############################################################################
## computation of asymptotic variance
###############################################################################
.HuMadrlsGetvar <- function(k){
    Var.loc <- (2*pnorm(k) - 1 - 2*k*dnorm(k) + 2*k^2*pnorm(-k))/(2*pnorm(k)-1)^2

    a.mad <- qnorm(0.75)
    Var.sc <- 1/(4*a.mad*dnorm(a.mad))^2

    return(Var.loc+Var.sc)
}

###############################################################################
## computation of maximum asymptotic MSE
###############################################################################
.HuMadrlsGetmse <- function(k, r){
    a.mad <- qnorm(0.75)
    b.mad <- 1/(4*a.mad*dnorm(a.mad))
    bias <- sqrt(k^2/(2*pnorm(k)-1)^2 + b.mad^2)

    Var <- .HuMadrlsGetvar(k = k)

    return(Var + r^2*bias^2)
}


###############################################################################
## optimal IC
###############################################################################
rlsOptIC.HuMad <- function(r, kUp = 2.5, delta = 1e-6){
    res <- optimize(f = .HuMadrlsGetmse, lower = 1e-4, upper = kUp, 
                tol = delta, r = r)
    k <- res$minimum
    a.mad <- qnorm(0.75)
    b.mad <- 1/(4*a.mad*dnorm(a.mad))
    bias <- sqrt(k^2/(2*pnorm(k)-1)^2 + b.mad^2)

    A.loc <- 1/(2*pnorm(k)-1)
    fct1 <- function(x){ A.loc*sign(x)*pmin(abs(x), k) }
    body(fct1) <- substitute({ A.loc*sign(x)*pmin(abs(x), k) },
                        list(k = k, A.loc = A.loc))
    fct2 <- function(x){ b.mad*sign(abs(x) - a.mad) }
    body(fct2) <- substitute({ b.mad*sign(abs(x) - a.mad) },
                        list(a.mad = a.mad, b.mad = b.mad))

    return(IC(name = "IC of HuMad type", 
              Curve = EuclRandVarList(RealRandVariable(Map = list(fct1, fct2), Domain = Reals())),
              Risks = list(asMSE = res$objective, asBias = bias, asCov = res$objective - r^2*bias^2), 
              Infos = matrix(c("rlsOptIC.HuMad", "optimally robust IC for HuMad estimators and 'asMSE'",
                               "rlsOptIC.HuMad", paste("where k =", round(k, 3))), 
                           ncol=2, byrow = TRUE, dimnames=list(character(0), c("method", "message"))), 
              CallL2Fam = call("NormLocationScaleFamily")))
}
