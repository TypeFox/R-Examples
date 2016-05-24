###############################################################################
## computation of bias
###############################################################################
.Tu1rlsGetbias <- function(x, a){
    A.loc <- 1/(2*a*dnorm(a)*(a^2 - 15) + (a^4 - 6*a^2 + 15)*(2*pnorm(a)-1))
    A.sc <- 1/(2*A.loc*integrate(f = function(x, a0){ 2*x^2*(a0^2-x^2)*(a0^2-3*x^2)*dnorm(x) }, 
                lower = 0, upper = a, rel.tol = .Machine$double.eps^0.5, a0 = a)$value)

    return(sqrt(x^2*(pmax((a^2-x^2),0))^4*A.loc^2 
                + (x^2*pmax((a^2-x^2),0)^2*A.loc - 1)^2*A.sc^2))
}


###############################################################################
## computation of asymptotic variance
###############################################################################
.Tu1rlsGetvar <- function(a){
    A.loc <- 1/(2*a*dnorm(a)*(a^2 - 15) + (a^4 - 6*a^2 + 15)*(2*pnorm(a)-1))
    h1 <- 2*integrate(f = function(x, a0){ x^2*(a0^2-x^2)^4*dnorm(x) }, lower = 0, 
                upper = a, rel.tol = .Machine$double.eps^0.5, a0 = a)$value

    A.sc <- 1/(2*A.loc*integrate(f = function(x, a0){ 2*x^2*(a0^2-x^2)*(a0^2-3*x^2)*dnorm(x) }, 
                lower = 0, upper = a, rel.tol = .Machine$double.eps^0.5, a0 = a)$value)
    h2 <- 2*integrate(f = function(x, a0){ x^4*(a0^2-x^2)^4*dnorm(x) }, lower = 0, 
                upper = a, rel.tol = .Machine$double.eps^0.5, a0 = a)$value*A.loc^2 - 1

    return(h1*A.loc^2 + h2*A.sc^2)
}

###############################################################################
## computation of maximum asymptotic MSE
###############################################################################
.Tu1rlsGetmse <- function(a, r){
    x <- seq(from = 0, to = a, by = 0.01)
    bias <- sapply(x, .Tu1rlsGetbias, a = a)  
    index <- which.max(bias)  
    if(index==length(x))
        b <- optimize(f=.Tu1rlsGetbias, lower=x[index-1], upper=x[index], 
                    maximum=TRUE, tol=.Machine$double.eps^0.5, a=a)$objective
    else
        b <- optimize(f=.Tu1rlsGetbias, lower=x[index-1], upper=x[index+1], 
                    maximum=TRUE, tol=.Machine$double.eps^0.5, a=a)$objective

    return(.Tu1rlsGetvar(a = a) + r^2*b^2)
}

###############################################################################
## optimal IC
###############################################################################
rlsOptIC.Tu1 <- function(r, aUp = 10, delta = 1e-6){
    res <- optimize(f = .Tu1rlsGetmse, lower = 1e-4, upper = aUp, 
                tol = delta, r = r)

    a <- res$minimum
    A.loc <- 1/(2*a*dnorm(a)*(a^2 - 15) + (a^4 - 6*a^2 + 15)*(2*pnorm(a)-1))
    A.sc <- 1/(2*A.loc*integrate(f = function(x, a0){ 2*x^2*(a0^2-x^2)*(a0^2-3*x^2)*dnorm(x) }, 
                lower = 0, upper = a, rel.tol = .Machine$double.eps^0.5, a0 = a)$value)

    x <- seq(from = 0, to = a, by = 0.01)
    bias <- sapply(x, .Tu1rlsGetbias, a = a)
    index <- which.max(bias)
    if(index==length(x))
        b <- optimize(f=.Tu1rlsGetbias, lower=x[index-1], upper=x[index], 
                    maximum=TRUE, tol=.Machine$double.eps^0.5, a=a)$objective
    else
        b <- optimize(f=.Tu1rlsGetbias, lower=x[index-1], upper=x[index+1], 
                    maximum=TRUE, tol=.Machine$double.eps^0.5, a=a)$objective

    fct1 <- function(x){ A.loc*x*(a^2 - x^2)^2*(abs(x) < a) }
    body(fct1) <- substitute({ A.loc*x*(a^2 - x^2)^2*(abs(x) < a) },
                        list(a = a, A.loc = A.loc))
    fct2 <- function(x){ A.sc*(A.loc*x^2*(a^2 - x^2)^2*(abs(x) < a) - 1) }
    body(fct2) <- substitute({ A.sc*(A.loc*x^2*(a^2 - x^2)^2*(abs(x) < a) - 1) },
                        list(a = a, A.loc = A.loc, A.sc = A.sc))

    return(IC(name = "IC of Tu1 type", 
              Curve = EuclRandVarList(RealRandVariable(Map = list(fct1, fct2), Domain = Reals())),
              Risks = list(asMSE = res$objective, asBias = b, asCov = res$objective - r^2*b^2), 
              Infos = matrix(c("rlsOptIC.Tu1", "optimally robust IC for Tu1 estimators and 'asMSE'",
                               "rlsOptIC.Tu1", paste("where a =", round(a, 3))), 
                           ncol=2, byrow = TRUE, dimnames=list(character(0), c("method", "message"))), 
              CallL2Fam = call("NormLocationScaleFamily")))
}
