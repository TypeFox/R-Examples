###############################################################################
## computation of a
###############################################################################
.BMrlsGeta <- function(a, bL, bS){
    tau1 <- bL/a; tau2 <- (1+bS)/bL; tau3 <- sqrt((1+bS)/a)

    if(tau1 <= tau3)
        return(2*(a*(pnorm(tau1)-0.5) - bL*dnorm(tau2) 
                + (1+bS)*pnorm(-tau2)) - 1)
    else
        return(2*(-a*tau3*dnorm(tau3) + a*(pnorm(tau3)-0.5) 
                + (1+bS)*pnorm(-tau3)) - 1)
}


###############################################################################
## computation of gg
###############################################################################
.BMrlsGetgg <- function(a, bL, bS){
    tau1 <- bL/a; tau2 <- (1+bS)/bL; tau3 <- sqrt((1+bS)/a)

    if(tau1 <= tau3)
        return(1/(2*(-bL*dnorm(tau1) + 3*a*(pnorm(tau1)-0.5) 
                  - 2*bL*dnorm(tau2) + (1+bS)*pnorm(-tau2)) - 1))
    else
        return(1/(2*(-3*a*tau3*dnorm(tau3) + 3*a*(pnorm(tau3)-0.5) 
                  + (1+bS)*pnorm(-tau3)) - 1))
}


###############################################################################
## computation of asymptotic variance E|rho|^2
###############################################################################
.BMrlsGetvar <- function(bL, bS, a, gg){
    tau1 <- bL/a; tau2 <- (1+bS)/bL; tau3 <- sqrt((1+bS)/a)

    if(tau1 <= tau3){
        res1 <- 2*(-a*bL*dnorm(tau1) + a^2*(pnorm(tau1)-0.5) + bL^2*(pnorm(tau2)-pnorm(tau1)) 
                    + bL*(1+bS)*dnorm(tau2) - (1+bS)^2*pnorm(-tau2))
        res2 <- 2*(-3*a*bL*dnorm(tau1) + 3*a^2*(pnorm(tau1)-0.5) - bL*(1+bS)*dnorm(tau2) 
                    + bL^2*(pnorm(tau2)-pnorm(tau1)) + (1+bS)^2*pnorm(-tau2))
    }else{
        res1 <- 2*(-a^2*tau3*dnorm(tau3) + a^2*(pnorm(tau3)-0.5) 
                    + (1+bS)^2*(dnorm(tau3)/tau3-pnorm(-tau3)))
        res2 <- 2*(-a^2*tau3^3*dnorm(tau3) + 3*a^2*(-tau3*dnorm(tau3) + pnorm(tau3)-0.5) 
                    + (1+bS)^2*pnorm(-tau3))
    }

    return(res1 + gg^2*(res2 - 1))
}


###############################################################################
## computation of minimal bL
###############################################################################
.BMrlsGetbLmin <- function(bL, bS){
    tau <- (1+bS)/bL
    return(2*bL*(1/sqrt(2*pi) - dnorm(tau) + tau*pnorm(-tau)) - 1)
}


###############################################################################
## computation of maximum asymptotic MSE
###############################################################################
.BMrlsGetmse <- function(bL.bS, r, MAX){
    bL <- bL.bS[1]; bS <- bL.bS[2]

    # constraint for bS
    if(bS < 1.0) return(MAX)

    bL.min <- uniroot(.BMrlsGetbLmin, lower = sqrt(pi/2), upper = 10, 
                    tol = .Machine$double.eps^0.5, bS = bS)$root

    # constraint for bL
    if(bL < bL.min) return(MAX)

    if(abs(bL-bL.min) < 1e-4){
        bL <- bL.min
        tau <- (1+bS)/bL
        gg <- 1/(4*bL*(1/sqrt(2*pi)-dnorm(tau) + 0.5*tau*pnorm(-tau)) - 1)
        bias.2 <- bL^2 + gg^2*bS^2

        Var <- 2*((1+gg^2)*bL^2*(pnorm(tau) - 0.5) + (1-gg^2)*bL*(1+bS)*dnorm(tau)
                    - (1-gg^2)*(1+bS)^2*pnorm(-tau)) - gg^2
    }else{
        a <- uniroot(.BMrlsGeta, lower = 0.01, upper = 1000, 
                    tol = .Machine$double.eps^0.5, bL = bL, bS = bS)$root
        gg <- .BMrlsGetgg(a = a, bL = bL, bS = bS)
        bias.2 <- bL^2 + gg^2*bS^2

        Var <- .BMrlsGetvar(bL = bL, bS = bS, a = a, gg = gg)
    }

    return(Var + r^2*bias.2)
}


###############################################################################
# optimal IC
###############################################################################
rlsOptIC.BM <- function(r, bL.start = 2, bS.start = 1.5, delta = 1e-6, MAX = 100){
    res <- optim(c(bL.start, bS.start), .BMrlsGetmse, method = "Nelder-Mead", 
                control = list(reltol=delta), r = r, MAX = MAX)
    bL <- res$par[1]; bS <- res$par[2]

    bL.min <- uniroot(.BMrlsGetbLmin, lower = sqrt(pi/2), upper = 10, 
                    tol = .Machine$double.eps^0.5, bS = bS)$root

    if(abs(bL-bL.min) < 1e-4){
        bL <- bL.min
        tau <- (1+bS)/bL
        gg <- 1/(4*bL*(1/sqrt(2*pi)-dnorm(tau) + 0.5*tau*pnorm(-tau)) - 1)
        a <- NULL
    }else{
        a <- uniroot(.BMrlsGeta, lower = 0.01, upper = 1000, 
                    tol = .Machine$double.eps^0.5, bL = bL, bS = bS)$root
        gg <- .BMrlsGetgg(a = a, bL = bL, bS = bS)
    }

    if(is.null(a)){
        fct1 <- function(x){ sign(x)*pmin(bL, (1+bS)/abs(x)) }
        body(fct1) <- substitute({ sign(x)*pmin(bL, (1+bS)/abs(x)) },
                            list(bL = bL, bS = bS))
        fct2 <- function(x){ gg*(abs(x)*pmin(bL, (1+bS)/abs(x)) - 1) }
        body(fct2) <- substitute({ gg*(abs(x)*pmin(bL, (1+bS)/abs(x)) - 1) },
                            list(bL = bL, bS = bS, gg = gg))
    }else{
        fct1 <- function(x){ sign(x)*pmin(a*abs(x), bL, (1+bS)/abs(x)) }
        body(fct1) <- substitute({ sign(x)*pmin(a*abs(x), bL, (1+bS)/abs(x)) },
                            list(bL = bL, bS = bS, a = a))
        fct2 <- function(x){ gg*(abs(x)*pmin(a*abs(x), bL, (1+bS)/abs(x)) - 1) }
        body(fct2) <- substitute({ gg*(abs(x)*pmin(a*abs(x), bL, (1+bS)/abs(x)) - 1) },
                            list(bL = bL, bS = bS, a = a, gg = gg))
    }

    b <- sqrt(bL^2 + gg^2*bS^2)
    return(IC(name = "IC of BM type", 
              Curve = EuclRandVarList(RealRandVariable(Map = list(fct1, fct2), Domain = Reals())),
              Risks = list(asMSE = res$value, asBias = b, asCov = res$value - r^2*b^2), 
              Infos = matrix(c("rlsOptIC.BM", "optimally robust IC for BM estimators and 'asMSE'",
                               "rlsOptIC.BM", paste(" where b.loc =", round(bL, 3), ", b.sc.0 =", round(bS, 3),
                                                    ", alpha =", round(a, 3), "and gamma =", round(gg, 3))), 
                           ncol=2, byrow = TRUE, dimnames=list(character(0), c("method", "message"))), 
              CallL2Fam = call("NormLocationScaleFamily")))
}
