###############################################################################
## psi function
###############################################################################
lsAn1.psi <- function(x, a){ sin(x/a)*(abs(x) < a*pi) }

###############################################################################
## chi function
###############################################################################
lsAn1.chi <- function(x, a){ 
    A.loc <- 1/(2*integrate(f = function(x, a0){cos(x/a0)*dnorm(x)/a0}, lower = 0, 
                    upper = a*pi, rel.tol = .Machine$double.eps^0.5, a0 = a)$value)

    return(A.loc*x*sin(x/a)*(abs(x) < a*pi) - 1)
}

###############################################################################
## computation of bias
###############################################################################
.An1rlsGetbias <- function(x, a){
    A.loc <- 1/(2*integrate(f = function(x, a0){cos(x/a0)*dnorm(x)/a0}, lower = 0, 
                    upper = a*pi, rel.tol = .Machine$double.eps^0.5, a0 = a)$value)
    Int1 <- integrate(f = function(x, a0){x*sin(x/a0)*dnorm(x)}, lower = 0, 
                    upper = a*pi, rel.tol = .Machine$double.eps^0.5, a0 = a)$value
    Int2 <- integrate(f = function(x, a0){x^2*cos(x/a0)/a0*dnorm(x)}, lower = 0, 
                    upper = a*pi, rel.tol = .Machine$double.eps^0.5, a0 = a)$value
    A.sc <- 1/(2*A.loc*(Int1 + Int2))

    return(sqrt(A.loc^2*(sin(x/a)*(abs(x) < a*pi))^2 
                + A.sc^2*(A.loc*x*sin(x/a)*(abs(x) < a*pi) - 1)^2))
}

###############################################################################
## computation of asymptotic variance
###############################################################################
.An1rlsGetvar <- function(a){
    A.loc <- 1/(2*integrate(f = function(x, a0){cos(x/a0)*dnorm(x)/a0}, lower = 0, 
                    upper = a*pi, rel.tol = .Machine$double.eps^0.5, a0 = a)$value)
    h1 <- 2*integrate(f = function(x, a0){sin(x/a0)^2*dnorm(x)}, lower = 0, 
                    upper = a*pi, rel.tol = .Machine$double.eps^0.5, a0 = a)$value

    Int1 <- integrate(f = function(x, a0){x*sin(x/a0)*dnorm(x)}, lower = 0, 
                    upper = a*pi, rel.tol = .Machine$double.eps^0.5, a0 = a)$value
    Int2 <- integrate(f = function(x, a0){x^2*cos(x/a0)/a0*dnorm(x)}, lower = 0, 
                    upper = a*pi, rel.tol = .Machine$double.eps^0.5, a0 = a)$value
    A.sc <- 1/(2*A.loc*(Int1 + Int2))
    h2 <- 2*integrate(f = function(x, a0){x^2*sin(x/a0)^2*dnorm(x)}, lower = 0, 
                    upper = a*pi, rel.tol = .Machine$double.eps^0.5, a0 = a)$value

    return(h1*A.loc^2 + A.sc^2*(h2*A.loc^2 - 1))
}

###############################################################################
## computation of maximum asymptotic MSE
###############################################################################
.An1rlsGetmse <- function(a, r){
    Var <- .An1rlsGetvar(a = a)

    x <- seq(from=0, to=a*pi, by = 0.01)
    bias <- sapply(x, .An1rlsGetbias, a = a)  
    index <- which.max(bias)  

    b <- optimize(f=.An1rlsGetbias, lower=x[index-1], upper=x[index+1], 
                maximum=TRUE, tol=1e-8, a=a)$objective

    return(Var + r^2*b^2)
}

###############################################################################
## optimal IC
###############################################################################
rlsOptIC.An1 <- function(r, aUp = 2.5, delta = 1e-6){
    res <- optimize(f = .An1rlsGetmse, lower = 1e-4, upper = aUp, 
                tol = delta, r = r)

    a <- res$minimum
    A.loc <- 1/(2*integrate(f = function(x, a0){cos(x/a0)*dnorm(x)/a0}, lower = 0, 
                    upper = a*pi, rel.tol = .Machine$double.eps^0.5, a0 = a)$value)
    Int1 <- integrate(f = function(x, a0){x*sin(x/a0)*dnorm(x)}, lower = 0, 
                    upper = a*pi, rel.tol = .Machine$double.eps^0.5, a0 = a)$value
    Int2 <- integrate(f = function(x, a0){x^2*cos(x/a0)/a0*dnorm(x)}, lower = 0, 
                    upper = a*pi, rel.tol = .Machine$double.eps^0.5, a0 = a)$value
    A.sc <- 1/(2*A.loc*(Int1 + Int2))

    x <- seq(from=0, to=a*pi, by = 0.01)
    bias <- sapply(x, .An1rlsGetbias, a = a)  
    index <- which.max(bias)  

    b <- optimize(f=.An1rlsGetbias, lower=x[index-1], upper=x[index+1], 
                maximum=TRUE, tol=1e-8, a=a)$objective

    fct1 <- function(x){ A.loc*sin(x/a)*(abs(x) < a*pi) }
    body(fct1) <- substitute({ A.loc*sin(x/a)*(abs(x) < a*pi) },
                        list(a = a, A.loc = A.loc))
    fct2 <- function(x){ A.sc*(A.loc*x*sin(x/a)*(abs(x) < a*pi) - 1) }
    body(fct2) <- substitute({ A.sc*(A.loc*x*sin(x/a)*(abs(x) < a*pi) - 1) },
                        list(a = a, A.loc = A.loc, A.sc = A.sc))

    return(IC(name = "IC of An1 type", 
              Curve = EuclRandVarList(RealRandVariable(Map = list(fct1, fct2), Domain = Reals())),
              Risks = list(asMSE = res$objective, asBias = b, asCov = res$objective - r^2*b^2), 
              Infos = matrix(c("rlsOptIC.An1", "optimally robust IC for An1 estimators and 'asMSE'",
                               "rlsOptIC.An1", paste("where a =", round(a, 3))), 
                           ncol=2, byrow = TRUE, dimnames=list(character(0), c("method", "message"))), 
              CallL2Fam = call("NormLocationScaleFamily")))
}
