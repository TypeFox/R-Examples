###############################################################################
## computation of bias
###############################################################################
.An2rlsGetbias <- function(x, a, k){
    beta.k <- 2*pnorm(k) - 1 - 2*k*dnorm(k) + 2*k^2*pnorm(-k)

    A.loc <- 1/(2*integrate(f = function(x, a0){ cos(x/a0)*dnorm(x)/a0 }, 
                    lower = 0, upper = a*pi, rel.tol = .Machine$double.eps^0.5, 
                    a0 = a)$value)
    return(sqrt(lsAn1.psi(x = x, a = a)^2*A.loc^2 
            + (pmin(k^2, x^2) - beta.k)^2/(2*(2*pnorm(k) - 1) - 4*k*dnorm(k))^2))
}

###############################################################################
## computation of asymptotic variance
###############################################################################
.An2rlsGetvar <- function(a, k){
    h1 <- 2*integrate(f=function(x, a0){ sin(x/a0)^2*dnorm(x) }, lower = 0, 
                upper = a*pi, rel.tol = .Machine$double.eps^0.5, a0 = a)$value
    A.loc <- 1/(2*integrate(f = function(x, a0){ cos(x/a0)*dnorm(x)/a0 }, 
                    lower = 0, upper = a*pi, rel.tol = .Machine$double.eps^0.5, 
                    a0 = a)$value)

    beta.k <- 2*pnorm(k) - 1 - 2*k*dnorm(k) + 2*k^2*pnorm(-k)
    E.psi.4 <- 3*(2*pnorm(k)-1) - 2*(k^3+3*k)*dnorm(k) + 2*k^4*pnorm(-k)
    A.sc <- 1/(2*(2*pnorm(k) - 1) - 4*k*dnorm(k))

    return(h1*A.loc^2 + (E.psi.4 - beta.k^2)*A.sc^2)
}

###############################################################################
## computation of maximum asymptotic MSE
###############################################################################
.An2rlsGetmse <- function(ak, r, MAX){
    a <- ak[1]; k <- ak[2]

    # constraints
    if(a < 0 || k < 0) return(MAX)

    Var <- .An2rlsGetvar(a = a, k = k)

    beta.k <- 2*pnorm(k) - 1 - 2*k*dnorm(k) + 2*k^2*pnorm(-k)
    b <- max(.An2rlsGetbias(x = 0, a = a, k = k), 
             .An2rlsGetbias(x = a*pi/2, a = a, k = k), 
             .An2rlsGetbias(x = a*pi, a = a, k = k), 
             .An2rlsGetbias(x = k, a = a, k = k), 
             .An2rlsGetbias(x = sqrt(beta.k), a = a, k = k))

    return(Var + r^2*b^2)
}

###############################################################################
## optimal IC
###############################################################################
rlsOptIC.An2 <- function(r, a.start = 1.5, k.start = 1.5, delta = 1e-6, MAX = 100){
    res <- optim(c(a.start, k.start), .An2rlsGetmse, method = "Nelder-Mead", 
                control = list(reltol=delta), r = r, MAX = MAX)

    a <- res$par[1]; k <- res$par[2]

    A.loc <- 1/(2*integrate(f = function(x, a0){ cos(x/a0)*dnorm(x)/a0 }, 
                    lower = 0, upper = a*pi, rel.tol = .Machine$double.eps^0.5, 
                    a0 = a)$value)
    beta.k <- 2*pnorm(k) - 1 - 2*k*dnorm(k) + 2*k^2*pnorm(-k)
    bias <- max(.An2rlsGetbias(x = 0, a = a, k = k), 
                .An2rlsGetbias(x = a*pi/2, a = a, k = k), 
                .An2rlsGetbias(x = a*pi, a = a, k = k), 
                .An2rlsGetbias(x = k, a = a, k = k), 
                .An2rlsGetbias(x = sqrt(beta.k), a = a, k = k))

    fct1 <- function(x){ A.loc*sin(x/a)*(abs(x) < a*pi) }
    body(fct1) <- substitute({ A.loc*sin(x/a)*(abs(x) < a*pi) },
                        list(a = a, A.loc = A.loc))
    A.sc <- 1/(2*(2*pnorm(k) - 1) - 4*k*dnorm(k))
    fct2 <- function(x){ A.sc*(pmin(x^2, k^2) - beta.k) }
    body(fct2) <- substitute({ A.sc*(pmin(x^2, k^2) - beta.k) },
                        list(k = k, beta.k = beta.k, 
                             A.sc = A.sc))
    return(IC(name = "IC of An2 type", 
              Curve = EuclRandVarList(RealRandVariable(Map = list(fct1, fct2), Domain = Reals())),
              Risks = list(asMSE = res$value, asBias = bias, asCov = res$value - r^2*bias^2), 
              Infos = matrix(c("rlsOptIC.An2", "optimally robust IC for An2 estimators and 'asMSE'",
                               "rlsOptIC.An2", paste("where a =", round(a, 3), "and k =", round(k, 3))), 
                           ncol=2, byrow = TRUE, dimnames=list(character(0), c("method", "message"))), 
              CallL2Fam = call("NormLocationScaleFamily")))
}
