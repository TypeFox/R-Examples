###############################################################################
## Computation of C_1
###############################################################################
.MrlsGetC1 <- function(gg){
    integrandC1 <- function(x, gg){ dnorm(x)/(1+gg^2*x^2) }
    return(2*integrate(integrandC1, lower = 0, upper = Inf, 
                    rel.tol = .Machine$double.eps^0.5, gg = gg)$value)
}


###############################################################################
## Computation of C_3 
###############################################################################
.MrlsGetC3 <- function(gg){
    integrandC3 <- function(x, gg){ dnorm(x)*x^4/(1+gg^2*x^2) }
    Int <- integrate(integrandC3, lower = 0, upper = Inf, 
                    rel.tol = .Machine$double.eps^0.5, gg = gg)$value

    return(1/gg + 1 - 2*gg^2*Int)
}


###############################################################################
## weight function
###############################################################################
.MrlsGetw <- function(x, gg, b, a1, a3){
    gvct <- (a1*x + a3*x^3)/sqrt(1 + gg^2*x^2)
    bvct <- sqrt(b^2 - gg^2/(1 + gg^2*x^2))
    absgvct <- abs(gvct)    
    ind1 <- (absgvct < bvct)

    return(ind1 + (1-ind1)*bvct/absgvct)
}


###############################################################################
## psi function
###############################################################################
.MrlsGetpsi <- function(x, gg, b, a1, a3){
    gvct <- (a1*x + a3*x^3)/sqrt(1+gg^2*x^2)
    wvct <- .MrlsGetw(x = x, gg = gg, b = b, a1 = a1, a3 = a3)

    return(gvct*wvct/sqrt(1+gg^2*x^2) + gg^2*x/(1+gg^2*x^2))
}


###############################################################################
## computation of r
###############################################################################
.MrlsGetr <- function(b, r, gg, a1, a3){
    integrandr <- function(x, gg, b, a1, a3){
        gvct <- (a1*x + a3*x^3)/sqrt(1 + gg^2*x^2)
        bvct <- sqrt(b^2 - gg^2/(1 + gg^2*x^2))
        hvct <- gvct/bvct - 1
        hvct <- (hvct > 0)*hvct

        return(hvct*dnorm(x))
    }
    Int <- integrate(integrandr, lower = -Inf, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, gg = gg, b = b, 
                a1 = a1, a3 = a3)$value

    return(r-sqrt(2*Int))
}


###############################################################################
## computation of a1 und a3 
###############################################################################
.MrlsGeta1a3 <- function(gg, b, a1, a3){
    C1 <- .MrlsGetC1(gg = gg)
    C3 <- .MrlsGetC3(gg = gg)

    integrand1 <- function(x, gg, b, a1, a3){
        dnorm(x)*x^2*.MrlsGetw(x, gg, b, a1, a3)/(1+gg^2*x^2)
    }
    Int1 <- 2*integrate(integrand1, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, gg = gg, b = b, 
                a1 = a1, a3 = a3)$value

    integrand2 <- function(x, gg, b, a1, a3){
        dnorm(x)*x^4*.MrlsGetw(x, gg, b, a1, a3)/(1+gg^2*x^2)
    }
    Int2 <- - 2*integrate(integrand2, lower = 0, upper = Inf,
                rel.tol = .Machine$double.eps^0.5, gg = gg, b = b, 
                a1 = a1, a3 = a3)$value

    integrand3 <- function(x, gg, b, a1, a3){
        dnorm(x)*x^6*.MrlsGetw(x, gg, b, a1, a3)/(1+gg^2*x^2)
    }
    Int3 <- 2*integrate(integrand3, lower = 0, upper = Inf,
                rel.tol = .Machine$double.eps^0.5, gg = gg, b = b, 
                a1 = a1, a3 = a3)$value

    D1 <- Int1*Int3 - Int2^2
    a1 <- (Int3*C1 + Int2*C3)/D1
    a3 <- (Int2*C1 + Int1*C3)/D1

    return(list(a1 = a1, a3 = a3))
}


###############################################################################
## computation of b, a1 and a3
###############################################################################
.MrlsGetba1a3 <- function(r, gg, a1.start, a3.start, bUp, delta, itmax, check){
    a1 <- a1.start; a3 <- a3.start
    b <- try(uniroot(.MrlsGetr, lower = gg, upper = bUp, 
            tol = .Machine$double.eps^0.5, r = r, gg = gg, 
            a1 = a1, a3 = a3)$root, silent = TRUE)
    if(!is.numeric(b)) b <- gg

    iter <- 0
    repeat{
        iter <- iter + 1
        if(iter > itmax){
            cat("Algorithm did not converge!\n")
            cat("=> increase itmax or try different starting values")
            cat("for 'a1' and 'a3'\n")
            return(c(NA, NA, NA))
        }
        a1.old <- a1; a3.old <- a3; b.old <- b

        a1a3 <- .MrlsGeta1a3(gg = gg, b = b, a1 = a1, a3 = a3)
        a1 <- a1a3$a1; a3 <- a1a3$a3

        b <- try(uniroot(.MrlsGetr, lower = gg, upper = bUp, 
                tol = .Machine$double.eps^0.5, r =r, gg = gg, 
                a1 = a1, a3 = a3)$root, silent = TRUE)
        if(!is.numeric(b)) b <- gg

        if(max(abs(b-b.old), abs(a1-a1.old), abs(a3-a3.old)) < delta)
            break
    }

    if(check){ # check constraints
        integrand1 <- function(x, gg, b, a1, a3){
            x*.MrlsGetpsi(x, gg, b, a1, a3)*dnorm(x)
        }
        Int1 <- integrate(integrand1, lower = -Inf, upper = Inf,
                    rel.tol = .Machine$double.eps^0.5, gg = gg, b = b, 
                    a1 = a1, a3 = a3)$value

        integrand2 <- function(x, gg, b, a1, a3){
            x^3*.MrlsGetpsi(x, gg, b, a1, a3)*dnorm(x)
        }
        Int2 <- integrate(integrand2, lower = -Inf, upper = Inf,
                    rel.tol = .Machine$double.eps^0.5, gg = gg, b = b, 
                    a1 = a1, a3 = a3)$value
        cat("constraint a1:\t", Int1-1, "\n")
        cat("constraint a3:\t", Int2-1/gg-1, "\n")

        rvgl <- .MrlsGetr(b = b, r = r, gg = gg, a1 = a1, a3 = a3)
        cat("MSE equation:\t", rvgl, "\n")
    }

    return(list(a1 = a1, a3 = a3, b = b))
}


###############################################################################
## computation of asymptotic variance \Ew|\rho|^2
###############################################################################
.MrlsGetvar <- function(gg, b, a1, a3){
    C1 <- .MrlsGetC1(gg = gg)

    integrandVar <- function(x, gg, b, a1, a3){
        dnorm(x)*(a1*x + a3*x^3)^2*.MrlsGetw(x, gg, b, a1, a3)^2/(1+gg^2*x^2)
    }
    Int <- integrate(integrandVar, lower = 0, upper = Inf,
                rel.tol = .Machine$double.eps^0.5, gg = gg, b = b, 
                a1 = a1, a3 = a3)$value

    return(2*Int + gg^2*C1)
}

###############################################################################
## computation of maximum asymptotic MSE
###############################################################################
.MrlsGetmse <- function(gg, r, a1.start, a3.start, bUp, delta, itmax, check){
    ba1a3 <- .MrlsGetba1a3(r = r, gg = gg, a1.start = a1.start, a3.start = a3.start, 
                    bUp = bUp, delta = delta, itmax = itmax, check = check)
    b <- ba1a3$b; a1 <- ba1a3$a1; a3 <- ba1a3$a3

    Var <- .MrlsGetvar(gg = gg, b = b, a1 = a1, a3 = a3)
    mse <- Var + r^2*b^2

#    cat("current gamma:\t", gg, "\tcurrent MSE:\t", mse, "\n")
    return(mse)
}


###############################################################################
## optimal IC
###############################################################################
rlsOptIC.M <- function(r, ggLo = 0.5, ggUp = 1.5, a1.start = 0.75, a3.start = 0.25, 
                     bUp = 1000, delta = 1e-5, itmax = 100, check = FALSE){
    res <- optimize(f = .MrlsGetmse, lower = ggLo, upper = ggUp, maximum = FALSE, 
                tol = .Machine$double.eps^0.25, r = r, a1.start = a1.start, 
                a3.start = a3.start, bUp = bUp, delta = delta, itmax = itmax, 
                check = 0)
    ba1a3 <- .MrlsGetba1a3(r = r, gg = res$minimum, a1.start = a1.start, 
                    a3.start = a3.start, bUp = bUp, delta = delta, 
                    itmax = itmax, check = check)

    a1 <- ba1a3$a1; a3 <- ba1a3$a3; b <- ba1a3$b
    gg <- res$minimum

    w <- .MrlsGetw
    fct1 <- function(x){ 
        wfct <- w 
        return((a1*x + a3*x^3)/(1 + gg^2*x^2)*wfct(x, gg = gg, b = b, a1 = a1, a3 = a3) 
               + gg^2*x/(1 + gg^2*x^2))
    }
    body(fct1) <- substitute({wfct <- w 
                              return((a1*x + a3*x^3)/(1 + gg^2*x^2)*wfct(x, gg = gg, b = b, a1 = a1, a3 = a3) 
                                     + gg^2*x/(1 + gg^2*x^2))},
                        list(w = .MrlsGetw, a1 = a1, a3 = a3, b = b, gg = gg))
    fct2 <- function(x){ 
        wfct <- w 
        return(gg*(x*((a1*x + a3*x^3)/(1 + gg^2*x^2)*wfct(x, gg = gg, b = b, a1 = a1, a3 = a3) 
                      + gg^2*x/(1 + gg^2*x^2)) - 1))
    }
    body(fct2) <- substitute({wfct <- w 
                              return(gg*(x*((a1*x + a3*x^3)/(1 + gg^2*x^2)*wfct(x, gg = gg, b = b, a1 = a1, a3 = a3) 
                                            + gg^2*x/(1 + gg^2*x^2)) - 1))},
                        list(w = .MrlsGetw, a1 = a1, a3 = a3, b = b, gg = gg))
    return(IC(name = "IC of M type", 
              Curve = EuclRandVarList(RealRandVariable(Map = list(fct1, fct2), Domain = Reals())),
              Risks = list(asMSE = res$objective, asBias = b, asCov = res$objective - r^2*b^2), 
              Infos = matrix(c("rlsOptIC.M", "optimally robust IC for M estimators and 'asMSE'",
                               "rlsOptIC.M", paste("where a1 =", round(a1, 3), ", a3 =", round(a3, 3),
                                                   ", b =", round(b, 3), "and gamma =", round(gg, 3))), 
                           ncol=2, byrow = TRUE, dimnames=list(character(0), c("method", "message"))), 
              CallL2Fam = call("NormLocationScaleFamily")))
}
