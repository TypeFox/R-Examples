###############################################################################
# weight function
###############################################################################
.MsrgsGetw <- function(u, z, gg, b.rg, b.sc, a1.x, a3){
    q.vkt <- sqrt(z + gg^2*u^2)
    g.vkt <- (a1.x*u + a3*u^3)/q.vkt^2

    Ind1 <- (b.rg/sqrt(z) < (1+b.sc)/abs(u))
    b.vkt <- Ind1*b.rg/sqrt(z) + (1-Ind1)*(1+b.sc)/abs(u)
    
    g.abs.vkt <- abs(g.vkt)
    Ind3 <- (g.abs.vkt < b.vkt)
    
    return(Ind3 + (1-Ind3)*b.vkt/g.abs.vkt)
}


###############################################################################
# computation of b.rg
###############################################################################
.MsrgsGetrrg1 <- function(z.a1, gg, b.rg, b.sc, a3){
    z <- z.a1[1]; a1.x <- z.a1[2]
    
    integrandrrg1 <- function(u, z, gg, b.rg, b.sc, a1.x, a3){
        h.vkt <- abs(a1.x*u + a3*u^3)/sqrt(z) - (z + gg^2*u^2)*b.rg/z

        return((h.vkt > 0)*h.vkt*dnorm(u))
    }

    lower <- -(1+b.sc)*sqrt(z)/b.rg
    upper <- (1+b.sc)*sqrt(z)/b.rg
    return(integrate(integrandrrg1, lower = lower, upper = upper, 
                rel.tol = .Machine$double.eps^0.5, z = z, gg = gg, b.rg = b.rg, 
                b.sc = b.sc, a1.x = a1.x, a3 = a3)$value)
}
.MsrgsGetrrg <- function(b.rg, b.sc, Z.a1, prob, r, gg, a3){
    r1 <- apply(Z.a1, 1, .MsrgsGetrrg1, gg = gg, b.rg = b.rg, 
                b.sc = b.sc, a3 = a3)    
    
    return(r-sqrt(sum(prob*r1)/b.rg))
}


###############################################################################
# computation of b.sc
###############################################################################
.MsrgsGetrsc1 <- function(z.a1, gg, b.rg, b.sc, a3){
    z <- z.a1[1]; a1.x <- z.a1[2]
    
    integrandrsc1 <- function(u, z, gg, b.rg, b.sc, a1.x, a3){
        h.vkt <- abs(a1.x + a3*u^2) - (1+b.sc)*(z + gg^2*u^2)/u^2

        return((h.vkt > 0)*h.vkt*dnorm(u))
    }

    lower <- -Inf
    upper <- -(1+b.sc)*sqrt(z)/b.rg
    Int1 <- integrate(integrandrsc1, lower = lower, upper = upper, 
                rel.tol = .Machine$double.eps^0.5, z = z, gg = gg, 
                b.rg = b.rg, b.sc = b.sc, a1.x = a1.x, a3 = a3)$value

    lower <- (1+b.sc)*sqrt(z)/b.rg
    upper <- Inf
    Int2 <- integrate(integrandrsc1, lower = lower, upper = upper, 
                rel.tol = .Machine$double.eps^0.5, z = z, gg = gg, 
                b.rg = b.rg, b.sc = b.sc, a1.x = a1.x, a3 = a3)$value

    return(Int1+Int2)
}
.MsrgsGetrsc <- function(b.sc, b.rg, Z.a1, prob, r, gg, a3){
    r1 <- apply(Z.a1, 1, .MsrgsGetrsc1, gg = gg, b.rg = b.rg, 
                b.sc = b.sc, a3 = a3)    
    
    return(r-sqrt(sum(prob*r1)/b.sc)/gg)
}


###############################################################################
# computation of alpha_1(x), alpha_3
###############################################################################
.MsrgsGeth1 <- function(z.a1, gg, b.rg, b.sc, a3){
    z <- z.a1[1]; a1.x <- z.a1[2]

    integrandh1 <- function(u, z, gg, b.rg, b.sc, a1.x, a3){
        u^2/(z+gg^2*u^2)*.MsrgsGetw(u, z, gg, b.rg, b.sc, a1.x, a3)*dnorm(u)
    }
    return(2*integrate(integrandh1, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, z = z, gg = gg, 
                b.rg = b.rg, b.sc = b.sc, a1.x = a1.x, a3 = a3)$value)
}
.MsrgsGeth2 <- function(z.a1, gg, b.rg, b.sc, a3){
    z <- z.a1[1]; a1.x <- z.a1[2]

    integrandh2 <- function(u, z, gg, b.rg, b.sc, a1.x, a3){
        u^4/(z+gg^2*u^2)*.MsrgsGetw(u, z, gg, b.rg, b.sc, a1.x, a3)*dnorm(u)
    }
    return(2*integrate(integrandh2, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, z = z, gg = gg, 
                b.rg = b.rg, b.sc = b.sc, a1.x = a1.x, a3 = a3)$value)
}
.MsrgsGeth3 <- function(z.a1, gg, b.rg, b.sc, a3){
    z <- z.a1[1]; a1.x <- z.a1[2]

    integrandh3 <- function(u, z, gg, b.rg, b.sc, a1.x, a3){
        u^6/(z+gg^2*u^2)*.MsrgsGetw(u, z, gg, b.rg, b.sc, a1.x, a3)*dnorm(u)
    }
    return(2*integrate(integrandh3, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, z = z, gg = gg, 
                b.rg = b.rg, b.sc = b.sc, a1.x = a1.x, a3 = a3)$value)
}
.MsrgsGeta1a3 <- function(Z.a1, prob, gg, b.rg, b.sc, a3){
    h1.x <- apply(Z.a1, 1, .MsrgsGeth1, gg = gg, b.rg = b.rg, b.sc = b.sc, a3 = a3)
    h2.x <- apply(Z.a1, 1, .MsrgsGeth2, gg = gg, b.rg = b.rg, b.sc = b.sc, a3 = a3)
    h3.x <- apply(Z.a1, 1, .MsrgsGeth3, gg = gg, b.rg = b.rg, b.sc = b.sc, a3 = a3)

    a1.x <- (1 - a3*h2.x)/h1.x

    h3 <- sum(prob*h3.x)
    h4 <- sum(prob*a1.x*h2.x)

    return(list(a1.x = a1.x, a3 = (1+1/gg - h4)/h3))
}


###############################################################################
# check of constraints
###############################################################################
.MsrgsGetch1 <- function(z.a1, gg, b.rg, b.sc, a3){
    z <- z.a1[1]; a1.x <- z.a1[2]

    integrand <- function(u, z, gg, b.rg, b.sc, a1.x, a3){
        (a1.x*u^2 + a3*u^4)/(z+gg^2*u^2)*.MsrgsGetw(u, z, gg, b.rg, b.sc, a1.x, a3)*dnorm(u)
    }
    return(2*integrate(integrand, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, z = z, gg = gg, 
                b.rg = b.rg, b.sc = b.sc, a1.x = a1.x, a3 = a3)$value)
}
.MsrgsGetch2 <- function(z.a1, gg, b.rg, b.sc, a3){
    z <- z.a1[1]; a1.x <- z.a1[2]

    integrand <- function(u, z, gg, b.rg, b.sc, a1.x, a3){
        (a1.x*u^4 + a3*u^6)/(z+gg^2*u^2)*.MsrgsGetw(u, z, gg, b.rg, b.sc, a1.x, a3)*dnorm(u)
    }
    return(2*integrate(integrand, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, z = z, gg = gg, 
                b.rg = b.rg, b.sc = b.sc, a1.x = a1.x, a3 = a3)$value)
}


###############################################################################
# computation of b, alpha_1(x), alpha_3
###############################################################################
.MsrgsGetba1a3 <- function(r, Z, prob, gg, a1.x, a3, bUp, b.sc, delta, itmax){
    Z.a1 <- matrix(c(Z, a1.x), ncol=2)

    b.rg <- try(uniroot(.MsrgsGetrrg, lower = 1e-2, upper = bUp, 
                    tol = .Machine$double.eps^0.5, Z.a1 = Z.a1, prob = prob, 
                    r = r, gg = gg, b.sc = b.sc, a3 = a3)$root, silent = TRUE)
    if(!is.numeric(b.rg)) b.rg <- sqrt(pi/2)
    
    b.sc <- try(uniroot(.MsrgsGetrsc, lower = 1e-2, upper = bUp, 
                    tol = .Machine$double.eps^0.5, Z.a1 = Z.a1, prob = prob, 
                    r = r, gg = gg, b.rg = b.rg, a3 = a3)$root, silent = TRUE)
    if(!is.numeric(b.sc)) b.sc <- 1
    if(b.sc < 1) b.sc <- 1

    b.rg <- try(uniroot(.MsrgsGetrrg, lower = 1e-2, upper = bUp, 
                    tol = .Machine$double.eps^0.5, Z.a1 = Z.a1, prob = prob, 
                    r = r, gg = gg, b.sc = b.sc, a3 = a3)$root, silent = TRUE)
    if(!is.numeric(b.rg)) b.rg <- sqrt(pi/2)
    
    iter <- 0
    repeat{
        iter <- iter + 1
        if(iter > itmax){
            cat("Algorithm did not converge!\n")
            cat("=> increase itmax or try different starting values")
            cat(" for 'a1.x' and 'a3'\n")
            return(NA)
        }
        a1.x.old <- a1.x; a3.old <- a3; b.rg.old <- b.rg
        b.sc.old <- b.sc

        a1.a3 <- try(.MsrgsGeta1a3(Z.a1 = Z.a1, prob = prob, gg = gg, 
                            b.rg = b.rg, b.sc = b.sc, a3 = a3), silent = TRUE)

        if(!is.list(a1.a3))
            a1.a3 <- list(a1.x=a1.x.old+1e-4, a3=a3.old+1e-4)

        a1.x <- a1.a3$a1.x; a3 <- a1.a3$a3
        Z.a1 <- matrix(c(Z, a1.x), ncol=2)

        b.rg <- try(uniroot(.MsrgsGetrrg, lower = 1e-2, upper = bUp, 
                        tol = .Machine$double.eps^0.5, Z.a1 = Z.a1, prob = prob, 
                        r = r, gg = gg, b.sc = b.sc, a3 = a3)$root, silent = TRUE)
        if(!is.numeric(b.rg)) b.rg <- b.rg.old + 1e-4
        b.rg.old <- b.rg  
    
        b.sc <- try(uniroot(.MsrgsGetrsc, lower = 1e-2, upper = bUp, 
                        tol = .Machine$double.eps^0.5, Z.a1 = Z.a1, prob = prob, 
                        r = r, gg = gg, b.rg = b.rg, a3 = a3)$root, silent = TRUE)
        if(!is.numeric(b.sc)) b.sc <- 1
        if(b.sc < 1) b.sc <- 1
        b.sc.old <- b.sc      

        b.rg <- try(uniroot(.MsrgsGetrrg, lower = 1e-2, upper = bUp, 
                        tol = .Machine$double.eps^0.5, Z.a1 = Z.a1, prob = prob, 
                        r = r, gg = gg, b.sc = b.sc, a3 = a3)$root, silent = TRUE)
        if(!is.numeric(b.rg)) b.rg <- sqrt(pi/2)
        
        if(max(abs(a1.x.old-a1.x), abs(a3.old-a3), abs(b.rg.old-b.rg), abs(b.sc.old-b.sc)) < delta)
            break
    }

    return(list(b.rg=b.rg, b.sc=b.sc, a1.x=a1.x, a3=a3))
}


###############################################################################
# computation of asymptotic variance
###############################################################################
.MsrgsGetvar <- function(z.a1, gg, b.rg, b.sc, a3, inteps){
    z <- z.a1[1]; a1.x <- z.a1[2]

    integrand <- function(u, z, gg, b.rg, b.sc, a1.x, a3){
        (a1.x*u + a3*u^3)^2/(z+gg^2*u^2)*.MsrgsGetw(u, z, gg, b.rg, b.sc, a1.x, a3)^2*dnorm(u)
    }
    return(2*integrate(integrand, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, z = z, gg = gg, 
                b.rg = b.rg, b.sc = b.sc, a1.x = a1.x, a3 = a3)$value)
}


###############################################################################
# computation of maximum asymptotic MSE
###############################################################################
.MsrgsGetmse<- function(gg, Z, prob, r, a1.x, a3, b.sc, bUp, delta, itmax){
    b.a1.a3 <- .MsrgsGetba1a3(r = r, Z = Z, prob = prob, gg = gg, 
                    a1.x = a1.x, a3 = a3, bUp = bUp, b.sc = b.sc, 
                    delta = delta, itmax = itmax)
    b.rg <- b.a1.a3$b.rg; b.sc <- b.a1.a3$b.sc; a1.x <- b.a1.a3$a1.x
    a3 <- b.a1.a3$a3

    Z.a1 <- matrix(c(Z, a1.x), ncol=2)
    
    var.x <- apply(Z.a1, 1, .MsrgsGetvar, gg = gg, b.rg = b.rg, 
                    b.sc = b.sc, a3=a3)
    mse <- sum(prob*var.x) - gg^2 + r^2*(b.rg^2 + gg^2*b.sc^2)

    cat("current gamma:\t", gg, "current mse:\t", mse, "\n")    
    return(mse)
}


###############################################################################
# computation of optimally robust IC
###############################################################################
rgsOptIC.Ms <- function(r, K, a1.x.start, a3.start = 0.25, b.sc.start = 1.5, 
                         bUp = 1000, ggLo = 0.5, ggUp = 1.0, delta = 1e-6, 
                         itmax = 1000, check = FALSE){
    if(!is(K, "DiscreteDistribution"))
        stop("Ms estimators are only implemented for\n",
             "1-dimensional discrete regressor distributions!")
    Reg2Mom <- as.vector(.rgsDesignTest(K = K))
    if(is.logical(Reg2Mom))
        stop("second moment matrix of regressor distribution 'K'", 
             "is (numerically) not positive definite")

    supp <- support(K); prob = d(K)(supp)
    A <- 1/Reg2Mom
    Z <- (A*supp)^2

    if(missing(a1.x.start))
        a1.x <- Z - 0.25
    else
        a1.x <- a1.x.start

    res <- optimize(.MsrgsGetmse, lower = ggLo, upper = ggUp, 
                tol = .Machine$double.eps^0.3, Z = Z, prob = prob, 
                r = r, a1.x = a1.x, a3 = a3.start, b.sc = b.sc.start, 
                bUp = bUp, delta = delta, itmax = itmax)

    gg <- res$minimum
    b.a1.a3 <- .MsrgsGetba1a3(r = r, Z = Z, prob = prob, gg = gg, 
                    a1.x = a1.x, a3 = a3.start, bUp = bUp, b.sc = b.sc.start, 
                    delta = delta, itmax = itmax)
    b.rg <- b.a1.a3$b.rg; b.sc <- b.a1.a3$b.sc; a1.x <- b.a1.a3$a1.x
    a3 <- b.a1.a3$a3

    if(check){
        Z.a1 <- matrix(c(Z, a1.x), byrow=F, ncol=2)

        kont1 <- try(apply(Z.a1, 1, .MsrgsGetch1, gg = gg, b.rg = b.rg, 
                            b.sc = b.sc, a3 = a3), silent = TRUE)
        if(is.numeric(kont1))
            cat("constraint (Mc2):\t", max(abs(kont1 - 1)), "\n")
        else
            cat("could not determine constraint (Mc2):\n", kont1, "\n")

        kont21 <- try(apply(Z.a1, 1, .MsrgsGetch2, gg = gg, b.rg = b.rg, 
                            b.sc = b.sc, a3 = a3), silent = TRUE)
        if(is.numeric(kont21)){ 
            kont2 <- sum(prob*kont21)
            cat("constraint (Mc4):\t", kont2 - 1 - 1/gg, "\n")
        }else
            cat("could not determine constraint (Mc4):\n", kont21, "\n")

        r.rg.vgl <- .MsrgsGetrrg(b.rg = b.rg, Z.a1 = Z.a1, prob = prob, r = r, 
                            gg = gg, b.sc = b.sc, a3 = a3)
        r.sc.vgl <- .MsrgsGetrsc(b.sc = b.sc, Z.a1 = Z.a1, prob = prob, r = r, 
                            gg = gg, b.rg = b.rg, a3 = a3)
                                
        cat("MSE equation for b.rg:\t", r.rg.vgl ,"\n")
        cat("MSE equation for b.sc:\t", r.sc.vgl ,"\n")
    }
    
    intervall <- getdistrOption("DistrResolution") / 2
    supp.grid <- as.numeric(matrix(
                      rbind(supp - intervall, supp + intervall), nrow = 1))
    a1.grid <- c(as.numeric(matrix(rbind(0, a1.x), nrow = 1)), 0)
    a1fct <- function(x){ stepfun(x = supp.grid, y = a1.grid)(x) }

    w <- .MsrgsGetw
    fct1 <- function(x){ 
        z <- (A*x[1])^2
        a1fct1 <- a1fct
        wfct <- w 
        a1.x <- a1fct1(x[1])
        w.vct <- wfct(u = x[2], z = z, gg = gg, b.rg = b.rg, b.sc = b.sc, a1.x = a1.x, a3 = a3) 
        psi <- (a1.x*x[2] + a3*x[2]^3)/(z + gg^2*x[2]^2)*w.vct
        return(A*x[1]*psi)
    }
    body(fct1) <- substitute({ z <- (A*x[1])^2
                               a1fct1 <- a1fct
                               wfct <- w 
                               a1.x <- a1fct1(x[1])
                               w.vct <- wfct(u = x[2], z = z, gg = gg, b.rg = b.rg, b.sc = b.sc, a1.x = a1.x, a3 = a3) 
                               psi <- (a1.x*x[2] + a3*x[2]^3)/(z + gg^2*x[2]^2)*w.vct
                               return(A*x[1]*psi)},
                        list(w = w, b.rg = b.rg, b.sc = b.sc, a1fct = a1fct, a3 = a3, A = A, gg = gg))
    fct2 <- function(x){ 
        z <- (A*x[1])^2
        a1fct1 <- a1fct
        wfct <- w 
        a1.x <- a1fct1(x[1])
        w.vct <- wfct(u = x[2], z = z, gg = gg, b.rg = b.rg, b.sc = b.sc, a1.x = a1.x, a3 = a3) 
        psi <- (a1.x*x[2] + a3*x[2]^3)/(z + gg^2*x[2]^2)*w.vct
        return(gg*(x[2]*psi - 1))
    }
    body(fct2) <- substitute({ z <- (A*x[1])^2
                               a1fct1 <- a1fct
                               wfct <- w 
                               a1.x <- a1fct1(x[1])
                               w.vct <- wfct(u = x[2], z = z, gg = gg, b.rg = b.rg, b.sc = b.sc, a1.x = a1.x, a3 = a3) 
                               psi <- (a1.x*x[2] + a3*x[2]^3)/(z + gg^2*x[2]^2)*w.vct
                               return(gg*(x[2]*psi - 1)) },
                        list(w = w, b.rg = b.rg, b.sc = b.sc, a1fct = a1fct, a3 = a3, A = A, gg = gg))
    bias.2 <- b.rg^2 + gg^2*b.sc^2
    return(CondIC(name = "conditionally centered IC of Ms type", 
              Curve = EuclRandVarList(RealRandVariable(Map = list(fct1, fct2), 
                                                       Domain = EuclideanSpace(dimension = 2))),
              Risks = list(asMSE = res$objective, asBias = sqrt(bias.2), trAsCov = res$objective - r^2*bias.2), 
              Infos = matrix(c("rgsOptIC.Ms", "optimally robust IC for Ms estimators and 'asMSE'",
                               "rgsOptIC.Ms", paste("where a3 =", round(a3, 3), ", b.rg =", round(b.rg, 3), 
                                                    ", b.sc =", round(b.sc, 3), "and gamma =", round(gg, 3))), 
                           ncol=2, byrow = TRUE, dimnames=list(character(0), c("method", "message"))), 
              CallL2Fam = call("NormLinRegScaleFamily", theta = 0, RegDistr = K, 
                                                         Reg2Mom = as.matrix(Reg2Mom))))
}
