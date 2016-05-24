###############################################################################
## computation of bias
###############################################################################
.Ha4rlsGetbias <- function(x, k, a, b, c0){
    beta.k <- 2*pnorm(k) - 1 - 2*k*dnorm(k) + 2*k^2*pnorm(-k)
    A.loc <- 2*(pnorm(a) - 0.5 - a/(c0-b)*(pnorm(c0)-pnorm(b)))

    return(sqrt(.Ha3rlsGetpsi(x = x, a = a, b = b, c0 = c0)^2/A.loc^2 
        + (pmin(k^2, x^2) - beta.k)^2/(2*(2*pnorm(k) - 1) - 4*k*dnorm(k))^2))
}


###############################################################################
## computation of asymptotic variance
###############################################################################
.Ha4rlsGetvar <- function(a, b, c0, k){
    h1 <- 2*(-a*dnorm(a) + pnorm(a) - 0.5 + a^2*(pnorm(b)-pnorm(a)) 
             + (a/(c0-b))^2*(c0^2*(pnorm(c0)-pnorm(b)) 
             + 2*c0*(dnorm(c0)-dnorm(b)) - c0*dnorm(c0) + b*dnorm(b) 
             + pnorm(c0) - pnorm(b)))
    A.loc <- 1/(2*(pnorm(a) - 0.5 - a/(c0-b)*(pnorm(c0)-pnorm(b))))
    Var.loc <- h1*A.loc^2

    beta.k <- 2*pnorm(k) - 1 - 2*k*dnorm(k) + 2*k^2*pnorm(-k)
    E.psi.4 <- 3*(2*pnorm(k)-1) - 2*(k^3+3*k)*dnorm(k) + 2*k^4*pnorm(-k)
    Var.sc <- (E.psi.4 - beta.k^2)/(2*(2*pnorm(k) - 1) - 4*k*dnorm(k))^2

    return(Var.loc+Var.sc)
}


###############################################################################
## computation of maximum asymptotic MSE
###############################################################################
.Ha4rlsGetmse <- function(abc0k, r, MAX){
    a <- abc0k[1]; b <- abc0k[2]; c0 <- abc0k[3]; k <- abc0k[4]

    #constraints
    if(a < 0 || a > b || b > c0 || k < 0) return(MAX)

    Var <- .Ha4rlsGetvar(a = a, b = b, c0 = c0, k = k)

    beta.k <- 2*pnorm(k) - 1 - 2*k*dnorm(k) + 2*k^2*pnorm(-k)
    bias <- max(.Ha4rlsGetbias(x = 0, k = k, a = a, b = b, c0 = c0), 
                .Ha4rlsGetbias(x = a, k = k, a = a, b = b, c0 = c0), 
                .Ha4rlsGetbias(x = b, k = k, a = a, b = b, c0 = c0), 
                .Ha4rlsGetbias(x = c0, k = k, a = a, b = b, c0 = c0), 
                .Ha4rlsGetbias(x = sqrt(beta.k), k = k, a = a, b = b, c0 = c0), 
                .Ha4rlsGetbias(x = k, k = k, a = a, b = b, c0 = c0))

    return(Var + r^2*bias^2)
}

###############################################################################
## optimal IC
###############################################################################
rlsOptIC.Ha4 <- function(r, a.start = 0.25, b.start = 2.5, c.start = 5.0, 
                        k.start = 1.0, delta = 1e-6, MAX = 100){
    res <- optim(c(a.start, b.start, c.start, k.start), .Ha4rlsGetmse, 
                method = "Nelder-Mead", control = list(reltol=delta), 
                r = r, MAX = MAX)

    a <- res$par[1]; b <- res$par[2]; c0 <- res$par[3]; k <- res$par[4] 
    A.loc <- 1/(2*(pnorm(a) - 0.5 - a/(c0-b)*(pnorm(c0)-pnorm(b))))
    beta.k <- 2*pnorm(k) - 1 - 2*k*dnorm(k) + 2*k^2*pnorm(-k)
    A.sc <- 1/(2*(2*pnorm(k) - 1) - 4*k*dnorm(k))
    bias <- max(.Ha4rlsGetbias(x = 0, k = k, a = a, b = b, c0 = c0), 
                .Ha4rlsGetbias(x = a, k = k, a = a, b = b, c0 = c0), 
                .Ha4rlsGetbias(x = b, k = k, a = a, b = b, c0 = c0), 
                .Ha4rlsGetbias(x = c0, k = k, a = a, b = b, c0 = c0), 
                .Ha4rlsGetbias(x = sqrt(beta.k), k = k, a = a, b = b, c0 = c0), 
                .Ha4rlsGetbias(x = k, k = k, a = a, b = b, c0 = c0))

    fct1 <- function(x){ Ind1 <- (abs(x) < a); Ind2 <- (abs(x) < b); Ind3 <- (abs(x) < c0)
                         A.loc*(x*Ind1 + a*sign(x)*(Ind2-Ind1) + a*(c0-abs(x))/(c0-b)*sign(x)*(Ind3-Ind2))}
    body(fct1) <- substitute({ Ind1 <- (abs(x) < a); Ind2 <- (abs(x) < b); Ind3 <- (abs(x) < c0)
                         A.loc*(x*Ind1 + a*sign(x)*(Ind2-Ind1) + a*(c0-abs(x))/(c0-b)*sign(x)*(Ind3-Ind2)) },
                        list(A.loc = A.loc, a = a, b = b, c0 = c0))
    fct2 <- function(x){ A.sc*(pmin(x^2, k^2) - beta.k) }
    body(fct2) <- substitute({ A.sc*(pmin(x^2, k^2) - beta.k) },
                        list(k = k, beta.k = beta.k, A.sc = A.sc))

    return(IC(name = "IC of Ha4 type", 
              Curve = EuclRandVarList(RealRandVariable(Map = list(fct1, fct2), Domain = Reals())),
              Risks = list(asMSE = res$value, asBias = bias, asCov = res$value - r^2*bias^2), 
              Infos = matrix(c("rlsOptIC.Ha4", "optimally robust IC for Ha4 estimators and 'asMSE'",
                               "rlsOptIC.Ha4", paste("where a =", round(a, 3), ", b =", round(b, 3),
                                                     ", c =", round(c0, 3), "and k =", round(k, 3))), 
                           ncol=2, byrow = TRUE, dimnames=list(character(0), c("method", "message"))), 
              CallL2Fam = call("NormLocationScaleFamily")))
}
