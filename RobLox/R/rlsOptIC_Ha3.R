###############################################################################
## psi function
###############################################################################
.Ha3rlsGetpsi <- function(x, a, b, c0){
    Ind1 <- (abs(x) < a); Ind2 <- (abs(x) < b); Ind3 <- (abs(x) < c0)

    return(x*Ind1 + a*sign(x)*(Ind2-Ind1) + a*(c0-abs(x))/(c0-b)*sign(x)*(Ind3-Ind2))
}

###############################################################################
## chi function
###############################################################################
lsHa3.chi <- function(x, a, b, c0){
    A.loc <- 1/(2*(pnorm(a) - 0.5 - a/(c0-b)*(pnorm(c0)-pnorm(b))))

    return(A.loc*x*.Ha3rlsGetpsi(x = x, a = a, b = b, c0 = c0) - 1)
}

###############################################################################
## computation of bias
###############################################################################
.Ha3rlsGetbias <- function(x, a, b, c0){
    A.loc <- 1/(2*(pnorm(a) - 0.5 - a/(c0-b)*(pnorm(c0)-pnorm(b))))
    A.sc <- 1/(2*A.loc*(2*(-a*dnorm(a) + pnorm(a) - 0.5) - a*(dnorm(b)-dnorm(a)) 
                + a/(c0-b)*(c0*dnorm(c0)+(c0-2*b)*dnorm(b) + 2*(pnorm(b) - pnorm(c0)))))

    return(sqrt(.Ha3rlsGetpsi(x = x, a = a, b = b, c0 = c0)^2*A.loc^2 
                + (x*.Ha3rlsGetpsi(x = x, a = a, b = b, c0 = c0)*A.loc - 1)^2*A.sc^2))
}


###############################################################################
## computation of asymptotic variance
###############################################################################
.Ha3rlsGetvar <- function(a, b, c0){
    h1 <- 2*(-a*dnorm(a) + pnorm(a) - 0.5 + a^2*(pnorm(b)-pnorm(a)) 
             + (a/(c0-b))^2*(c0^2*(pnorm(c0)-pnorm(b))
             + 2*c0*(dnorm(c0)-dnorm(b)) - c0*dnorm(c0) + b*dnorm(b) 
             + pnorm(c0) - pnorm(b)))
    A.loc <- 1/(2*(pnorm(a) - 0.5 - a/(c0-b)*(pnorm(c0)-pnorm(b))))
    Var.loc <- h1*A.loc^2

    A.sc <- 1/(2*A.loc*(2*(-a*dnorm(a) + pnorm(a) - 0.5) - a*(dnorm(b)-dnorm(a)) 
                + a/(c0-b)*(c0*dnorm(c0)+(c0-2*b)*dnorm(b) 
                + 2*(pnorm(b) - pnorm(c0)))))
    h2 <- 2*(-3*a*dnorm(a) - a^2*b*dnorm(b) + 3*(pnorm(a)-0.5) 
             + a^2*(pnorm(b)-pnorm(a)) + (a/(c0-b))^2*(c0*dnorm(c0) 
             + (c0^2*b-2*c0*b^2-4*c0+b^3+3*b)*dnorm(b) 
             + (c0^2+3)*(pnorm(c0)-pnorm(b))))*A.loc^2 - 1
    Var.sc <- h2*A.sc^2

    return(Var.loc+Var.sc)
}

###############################################################################
## computation of maximum asymptotic MSE
###############################################################################
.Ha3rlsGetmse <- function(abc0, r, MAX){
    a <- abc0[1]; b <- abc0[2]; c0 <- abc0[3]

    #constraints
    if(a < 0 || a > b || b > c0) return(MAX)

    Var <- .Ha3rlsGetvar(a = a, b = b, c0 = c0)

    x <- seq(from=0, to=c0, by=0.01)
    bias <- sapply(x, .Ha3rlsGetbias, a=a, b=b, c0=c0)  
    index <- which.max(bias)  
    if(index==length(x))
        b1 <- optimize(f=.Ha3rlsGetbias, lower=x[index-1], upper=x[index], 
                    maximum=TRUE, tol=1e-8, a=a, b=b, c0=c0)$objective
    else
        if(index==1)
            b1 <- optimize(f=.Ha3rlsGetbias, lower=x[index], upper=x[index+1], 
                        maximum=TRUE, tol=1e-8, a=a, b=b, c0=c0)$objective
        else
            b1 <- optimize(f=.Ha3rlsGetbias, lower=x[index-1], upper=x[index+1], 
                        maximum=TRUE, tol=1e-8, a=a, b=b, c0=c0)$objective

    return(Var + r^2*b1^2)
}

###############################################################################
## optimal IC
###############################################################################
rlsOptIC.Ha3 <- function(r, a.start = 0.25, b.start = 2.5, c.start = 5.0, 
                        delta = 1e-6, MAX = 100){
    res <- optim(c(a.start, b.start, c.start), .Ha3rlsGetmse, method = "Nelder-Mead", 
                control = list(reltol=delta), r = r, MAX = MAX)

    a <- res$par[1]; b <- res$par[2]; c0 <- res$par[3]
    A.loc <- 1/(2*(pnorm(a) - 0.5 - a/(c0-b)*(pnorm(c0)-pnorm(b))))
    A.sc <- 1/(2*A.loc*(2*(-a*dnorm(a) + pnorm(a) - 0.5) - a*(dnorm(b)-dnorm(a)) 
                + a/(c0-b)*(c0*dnorm(c0)+(c0-2*b)*dnorm(b) + 2*(pnorm(b) - pnorm(c0)))))

    x <- seq(from=0, to=c0, by=0.01)
    bias <- sapply(x, .Ha3rlsGetbias, a=a, b=b, c0=c0)  
    index <- which.max(bias)  
    if(index==length(x))
        b1 <- optimize(f=.Ha3rlsGetbias, lower=x[index-1], upper=x[index], 
                    maximum=TRUE, tol=1e-8, a=a, b=b, c0=c0)$objective
    else
        if(index==1)
            b1 <- optimize(f=.Ha3rlsGetbias, lower=x[index], upper=x[index+1], 
                        maximum=TRUE, tol=1e-8, a=a, b=b, c0=c0)$objective
        else
            b1 <- optimize(f=.Ha3rlsGetbias, lower=x[index-1], upper=x[index+1], 
                        maximum=TRUE, tol=1e-8, a=a, b=b, c0=c0)$objective

    fct1 <- function(x){ Ind1 <- (abs(x) < a); Ind2 <- (abs(x) < b); Ind3 <- (abs(x) < c0)
                         A.loc*(x*Ind1 + a*sign(x)*(Ind2-Ind1) + a*(c0-abs(x))/(c0-b)*sign(x)*(Ind3-Ind2))}
    body(fct1) <- substitute({ Ind1 <- (abs(x) < a); Ind2 <- (abs(x) < b); Ind3 <- (abs(x) < c0)
                         A.loc*(x*Ind1 + a*sign(x)*(Ind2-Ind1) + a*(c0-abs(x))/(c0-b)*sign(x)*(Ind3-Ind2)) },
                        list(A.loc = A.loc, a = a, b = b, c0 = c0))
    fct2 <- function(x){ Ind1 <- (abs(x) < a); Ind2 <- (abs(x) < b); Ind3 <- (abs(x) < c0)
                         A.sc*(A.loc*(x*Ind1 + a*sign(x)*(Ind2-Ind1) + a*(c0-abs(x))/(c0-b)*sign(x)*(Ind3-Ind2)) - 1) }
    body(fct2) <- substitute({ Ind1 <- (abs(x) < a); Ind2 <- (abs(x) < b); Ind3 <- (abs(x) < c0)
                         A.sc*(A.loc*x*(x*Ind1 + a*sign(x)*(Ind2-Ind1) + a*(c0-abs(x))/(c0-b)*sign(x)*(Ind3-Ind2)) - 1) },
                        list(A.loc = A.loc, A.sc = A.sc, a = a, b = b, c0 = c0))

    return(IC(name = "IC of Ha3 type", 
              Curve = EuclRandVarList(RealRandVariable(Map = list(fct1, fct2), Domain = Reals())),
              Risks = list(asMSE = res$value, asBias = b1, asCov = res$value - r^2*b1^2), 
              Infos = matrix(c("rlsOptIC.Ha3", "optimally robust IC for Ha3 estimators and 'asMSE'",
                               "rlsOptIC.Ha3", paste("where a =", round(a, 3), ", b =", round(b, 3),
                                                     "and c =", round(c0, 3))), 
                           ncol=2, byrow = TRUE, dimnames=list(character(0), c("method", "message"))), 
              CallL2Fam = call("NormLocationScaleFamily")))
}
