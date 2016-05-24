###############################################################################
## computation of asymptotic variance
###############################################################################
.AnMadrlsGetvar <- function(a){
    h1 <- 2*integrate(f = function(x, a0){ sin(x/a0)^2*dnorm(x) }, lower = 0, 
                    upper = a*pi, rel.tol = .Machine$double.eps^0.5, a0 = a)$value
    A.loc <- 1/(2*integrate(f = function(x, a0){ cos(x/a0)*dnorm(x)/a0 }, lower = 0, 
                    upper = a*pi, rel.tol = .Machine$double.eps^0.5, a0 = a)$value)
    a.mad <- qnorm(0.75)
    b.mad <- 1/(4*a.mad*dnorm(a.mad))

    return(A.loc^2*h1 + b.mad^2)
}

###############################################################################
## computation of maximum asymptotic MSE
###############################################################################
.AnMadrlsGetmse <- function(a, r){
    A.loc <- 1/(2*integrate(f = function(x, a0){ cos(x/a0)*dnorm(x)/a0 }, lower = 0, 
                    upper = a*pi, rel.tol = .Machine$double.eps^0.5, a0 = a)$value)
    a.mad <- qnorm(0.75)
    b.mad <- 1/(4*a.mad*dnorm(a.mad))

    return(.AnMadrlsGetvar(a = a) + r^2*(A.loc^2 + b.mad^2))
}

###############################################################################
## optimal IC
###############################################################################
rlsOptIC.AnMad <- function(r, aUp = 2.5, delta = 1e-6){
    res <- optimize(f = .AnMadrlsGetmse, lower = 1e-4, upper = aUp, 
                tol = delta, r = r)

    a <- res$minimum
    A.loc <- 1/(2*integrate(f = function(x, a0){ cos(x/a0)*dnorm(x)/a0 }, lower = 0, 
                    upper = a*pi, rel.tol = .Machine$double.eps^0.5, a0 = a)$value)
    a.mad <- qnorm(0.75)
    b.mad <- 1/(4*a.mad*dnorm(a.mad))
    bias <- sqrt(A.loc^2 + b.mad^2)

    fct1 <- function(x){ A.loc*sin(x/a)*(abs(x) < a*pi) }
    body(fct1) <- substitute({ A.loc*sin(x/a)*(abs(x) < a*pi) },
                        list(a = a, A.loc = A.loc))
    fct2 <- function(x){ b.mad*sign(abs(x) - a.mad) }
    body(fct2) <- substitute({ b.mad*sign(abs(x) - a.mad) },
                        list(a.mad = a.mad, b.mad = b.mad))

    return(IC(name = "IC of AnMad type", 
              Curve = EuclRandVarList(RealRandVariable(Map = list(fct1, fct2), Domain = Reals())),
              Risks = list(asMSE = res$objective, asBias = bias, asCov = res$objective - r^2*bias^2), 
              Infos = matrix(c("rlsOptIC.AnMad", "optimally robust IC for AnMad estimators and 'asMSE'",
                               "rlsOptIC.AnMad", paste("where a =", round(a, 3))), 
                           ncol=2, byrow = TRUE, dimnames=list(character(0), c("method", "message"))), 
              CallL2Fam = call("NormLocationScaleFamily")))
}
