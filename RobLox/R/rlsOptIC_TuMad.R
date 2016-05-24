###############################################################################
## computation of asymptotic variance
###############################################################################
.TuMadrlsGetvar <- function(a, r){
    A.loc <- 1/(2*a*dnorm(a)*(a^2 - 15) + (a^4 - 6*a^2 + 15)*(2*pnorm(a)-1))
    h1 <- 2*integrate(f = function(x, a0){ x^2*(a0^2-x^2)^4*dnorm(x) }, lower = 0, 
                upper = a, rel.tol = .Machine$double.eps^0.5, a0 = a)$value
    a.mad <- qnorm(0.75)
    b.mad <- 1/(4*a.mad*dnorm(a.mad))

    return(A.loc^2*h1 + b.mad^2)
}

###############################################################################
## computation of maximum asymptotic MSE
###############################################################################
.TuMadrlsGetmse <- function(a, r){
    A.loc <- 1/(2*a*dnorm(a)*(a^2 - 15) + (a^4 - 6*a^2 + 15)*(2*pnorm(a)-1))
    a.mad <- qnorm(0.75)
    b.mad <- 1/(4*a.mad*dnorm(a.mad))

    b.2 <- sqrt(A.loc^2*(16/25/sqrt(5)*a^5)^2 + b.mad^2)^2

    return(.TuMadrlsGetvar(a = a) + r^2*b.2)
}

###############################################################################
## optimal IC
###############################################################################
rlsOptIC.TuMad <- function(r, aUp = 10, delta = 1e-6){
    res <- optimize(f = .TuMadrlsGetmse, lower = 1e-4, upper = aUp, 
                tol = delta, r = r)

    a <- res$minimum
    A.loc <- 1/(2*a*dnorm(a)*(a^2 - 15) + (a^4 - 6*a^2 + 15)*(2*pnorm(a)-1))
    a.mad <- qnorm(0.75)
    b.mad <- 1/(4*a.mad*dnorm(a.mad))

    bias <- sqrt(A.loc^2*(16/25/sqrt(5)*a^5)^2 + b.mad^2)

    fct1 <- function(x){ A.loc*x*(a^2 - x^2)^2*(abs(x) < a) }
    body(fct1) <- substitute({ A.loc*x*(a^2 - x^2)^2*(abs(x) < a) },
                        list(a = a, A.loc = A.loc))
    fct2 <- function(x){ b.mad*sign(abs(x) - a.mad) }
    body(fct2) <- substitute({ b.mad*sign(abs(x) - a.mad) },
                        list(a.mad = a.mad, b.mad = b.mad))

    return(IC(name = "IC of TuMad type", 
              Curve = EuclRandVarList(RealRandVariable(Map = list(fct1, fct2), Domain = Reals())),
              Risks = list(asMSE = res$objective, asBias = bias, asCov = res$objective - r^2*bias^2), 
              Infos = matrix(c("rlsOptIC.TuMad", "optimally robust IC for TuMad estimators and 'asMSE'",
                               "rlsOptIC.TuMad", paste("where a =", round(a, 3))), 
                           ncol=2, byrow = TRUE, dimnames=list(character(0), c("method", "message"))), 
              CallL2Fam = call("NormLocationScaleFamily")))
}
