###############################################################################
## psi function
###############################################################################
.MM2rlsGetpsi <- function(x, c0){
    Ind1 <- (abs(x) <= 0.8*c0); Ind2 <- (abs(x) <= 1.0*c0)
    return(sign(x)*(Ind1*abs(x)/c0 + (1-Ind1)*Ind2*(38.4 - 175*abs(x)/c0 
            + 300*(abs(x)/c0)^2 - 225*(abs(x)/c0)^3 + 62.5*(abs(x)/c0)^4) 
            + (1-Ind2)*0.9))
}

###############################################################################
## chi function
###############################################################################
.MM2rlsGetchi <- function(x, d){
    Ind1 <- (abs(x) <= d)
    return(Ind1*(3*(x/d)^2 - 3*(x/d)^4 + (x/d)^6) + (1-Ind1))
}

###############################################################################
## computation of bias
###############################################################################
.MM2rlsGetbias <- function(x, c0, d, A.loc, beta.d, A.sc){
    return(sqrt(A.loc^2*.MM2rlsGetpsi(x = x, c0 = c0)^2 
                + A.sc^2*(.MM2rlsGetchi(x = x, d = d) - beta.d)^2))
}


###############################################################################
## computation of asymptotic variance
###############################################################################
.MM2rlsGetvar <- function(c0, d){
    A.loc <- 1/(2*integrate(f = function(x, c0){ x*.MM2rlsGetpsi(x = x, c0 = c0)*dnorm(x) }, 
                    lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5, 
                    c0 = c0)$value)
    Var.loc <- A.loc^2*2*integrate(f = function(x, c0){ .MM2rlsGetpsi(x = x, c0 = c0)^2*dnorm(x) }, 
                    lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5, 
                    c0 = c0)$value

    beta.d <- 2*integrate(f = function(x, d){ .MM2rlsGetchi(x = x, d = d)*dnorm(x) }, 
                    lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5, 
                    d = d)$value
    A.sc <- 1/(2*integrate(f = function(x, d){ x^2*.MM2rlsGetchi(x = x, d = d)*dnorm(x) }, 
                    lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5, 
                    d = d)$value - beta.d)
    Var.sc <- A.sc^2*(2*integrate(f = function(x, d){ .MM2rlsGetchi(x = x, d = d)^2*dnorm(x) }, 
                    lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5, 
                    d = d)$value - beta.d^2)

    return(Var.loc + Var.sc)
}

###############################################################################
## computation of maximum asymptotic MSE
###############################################################################
.MM2rlsGetmse <- function(c0d, r, MAX){
    c0 <- c0d[1]; d <- c0d[2]

    # constraints
    if(c0 < 0 || d < 0) return(MAX)

    A.loc <- 1/(2*integrate(f = function(x, c0){ x*.MM2rlsGetpsi(x = x, c0 = c0)*dnorm(x) }, 
                    lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5, 
                    c0 = c0)$value)
    beta.d <- 2*integrate(f = function(x, d){ .MM2rlsGetchi(x = x, d = d)*dnorm(x) }, 
                    lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5, 
                    d = d)$value
    A.sc <- 1/(2*integrate(f = function(x, d){ x^2*.MM2rlsGetchi(x = x, d = d)*dnorm(x) }, 
                    lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5, 
                    d = d)$value - beta.d)

    x <- seq(from=0, to=max(c0,d), by = 0.01)
    bias <- sapply(x, .MM2rlsGetbias, c0 = c0, d = d, A.loc = A.loc, 
                    beta.d = beta.d, A.sc = A.sc)  
    index <- which.max(bias)  

    if(index==length(x))
        b1 <- optimize(f=.MM2rlsGetbias, lower=x[index-1], upper=x[index], maximum=TRUE, 
                    tol = .Machine$double.eps^0.5, c0=c0, d=d, A.loc = A.loc, 
                    beta.d = beta.d, A.sc = A.sc)$objective
    else{
        if(index==1)
            b1 <- optimize(f=.MM2rlsGetbias, lower=x[index], upper=x[index+1], maximum=TRUE, 
                        tol = .Machine$double.eps^0.5, c0=c0, d=d, A.loc = A.loc, 
                        beta.d = beta.d, A.sc = A.sc)$objective
        else 
            b1 <- optimize(f=.MM2rlsGetbias, lower=x[index-1], upper=x[index+1], maximum=TRUE, 
                        tol = .Machine$double.eps^0.5, c0=c0, d=d, A.loc = A.loc, 
                        beta.d = beta.d, A.sc = A.sc)$objective
    }

    return(.MM2rlsGetvar(c0 = c0, d = d) + r^2*b1^2)
}

###############################################################################
## optimal IC
###############################################################################
rlsOptIC.MM2 <- function(r, c.start = 1.5, d.start = 2.0, delta = 1e-6, MAX = 100){
    res <- optim(c(c.start, d.start), .MM2rlsGetmse, method = "Nelder-Mead", 
                control = list(reltol=delta), r = r, MAX = MAX)

    c0 <- res$par[1]; d <- res$par[2]

    A.loc <- 1/(2*integrate(f = function(x, c0){ x*.MM2rlsGetpsi(x = x, c0 = c0)*dnorm(x) }, 
                    lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5, 
                    c0 = c0)$value)
    beta.d <- 2*integrate(f = function(x, d){ .MM2rlsGetchi(x = x, d = d)*dnorm(x) }, 
                    lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5, 
                    d = d)$value
    A.sc <- 1/(2*integrate(f = function(x, d){ x^2*.MM2rlsGetchi(x = x, d = d)*dnorm(x) }, 
                    lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5, 
                    d = d)$value - beta.d)

    x <- seq(from=0, to=max(c0,d), by = 0.01)
    bias <- sapply(x, .MM2rlsGetbias, c0 = c0, d = d, A.loc = A.loc, 
                    beta.d = beta.d, A.sc = A.sc)
    index <- which.max(bias)

    if(index==length(x))
        b1 <- optimize(f=.MM2rlsGetbias, lower=x[index-1], upper=x[index], maximum=TRUE, 
                    tol = .Machine$double.eps^0.5, c0=c0, d=d, A.loc = A.loc, 
                    beta.d = beta.d, A.sc = A.sc)$objective
    else{
        if(index==1)
            b1 <- optimize(f=.MM2rlsGetbias, lower=x[index], upper=x[index+1], maximum=TRUE, 
                        tol = .Machine$double.eps^0.5, c0=c0, d=d, A.loc = A.loc, 
                        beta.d = beta.d, A.sc = A.sc)$objective
        else 
            b1 <- optimize(f=.MM2rlsGetbias, lower=x[index-1], upper=x[index+1], maximum=TRUE, 
                        tol = .Machine$double.eps^0.5, c0=c0, d=d, A.loc = A.loc, 
                        beta.d = beta.d, A.sc = A.sc)$objective
    }

    fct1 <- function(x){ Ind1 <- (abs(x) <= 0.8*c0); Ind2 <- (abs(x) <= 1.0*c0)
                         A.loc*sign(x)*(Ind1*abs(x)/c0 + (1-Ind1)*Ind2*(38.4 - 175*abs(x)/c0 
                                    + 300*(abs(x)/c0)^2 - 225*(abs(x)/c0)^3 
                                    + 62.5*(abs(x)/c0)^4) + (1-Ind2)*0.9) }
    body(fct1) <- substitute({ Ind1 <- (abs(x) <= 0.8*c0); Ind2 <- (abs(x) <= 1.0*c0)
                               A.loc*sign(x)*(Ind1*abs(x)/c0 + (1-Ind1)*Ind2*(38.4 - 175*abs(x)/c0 
                                          + 300*(abs(x)/c0)^2 - 225*(abs(x)/c0)^3 
                                          + 62.5*(abs(x)/c0)^4) + (1-Ind2)*0.9) },
                        list(c0 = c0, A.loc = A.loc))
    fct2 <- function(x){ Ind1 <- (abs(x) <= d)
                         A.sc*(Ind1*(3*(x/d)^2 - 3*(x/d)^4 + (x/d)^6) + (1-Ind1) - beta.d) }
    body(fct2) <- substitute({ Ind1 <- (abs(x) <= d)
                               A.sc*(Ind1*(3*(x/d)^2 - 3*(x/d)^4 + (x/d)^6) + (1-Ind1) - beta.d) },
                        list(d = d, beta.d = beta.d, A.sc = A.sc))
    return(IC(name = "IC of MM2 type", 
              Curve = EuclRandVarList(RealRandVariable(Map = list(fct1, fct2), Domain = Reals())),
              Risks = list(asMSE = res$value, asBias = b1, asCov = res$value - r^2*b1^2), 
              Infos = matrix(c("rlsOptIC.MM2", "optimally robust IC for MM2 estimators and 'asMSE'",
                               "rlsOptIC.MM2", paste("where c =", round(c0, 3), "and d =", round(d, 3))), 
                           ncol=2, byrow = TRUE, dimnames=list(character(0), c("method", "message"))), 
              CallL2Fam = call("NormLocationScaleFamily")))
}
