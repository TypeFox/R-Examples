###############################################################################
# weight function
###############################################################################
.McrgsGetw <- function(u, z, gg, b, a1.x, a3){
    q.vkt <- sqrt(z + gg^2*u^2)
    g.vkt <- (a1.x*u + a3*u^3)/q.vkt
    b.vkt <- sqrt(b^2 - gg^2*z/q.vkt^2)
    
    g.abs.vkt <- abs(g.vkt)
    ind1 <- (g.abs.vkt < b.vkt)   
    
    return(ind1 + (1-ind1)*b.vkt/g.abs.vkt)
}


###############################################################################
# computation of b
###############################################################################
.McrgsGetr1 <- function(z.a1, gg, b, a3){
    z <- z.a1[1]; a1.x <- z.a1[2]
    
    integrandr1 <- function(u, z, b, a1.x, a3){
        q.vkt <- sqrt(z + gg^2*u^2)
        g.vkt <- (a1.x*u + a3*u^3)/q.vkt
        b.vkt <- sqrt(b^2 - gg^2*z/q.vkt^2)
        h.vkt <- abs(g.vkt)/b.vkt - 1
        
        return((h.vkt > 0)*h.vkt*dnorm(u))
    }
    return(2*integrate(integrandr1, lower = 0, upper = Inf, 
                    rel.tol = .Machine$double.eps^0.5, z = z, b = b, 
                    a1.x = a1.x, a3 = a3)$value)
}
.McrgsGetr <- function(b, Z.a1, prob, r, gg, a3){
    r1 <- apply(Z.a1, 1, .McrgsGetr1, gg = gg, b = b, a3 = a3)
    
    return(r-sqrt(sum(prob*r1)))
}


###############################################################################
# computation of C1.x
###############################################################################
C1.x.fkt <- function(z, gg){
    integrandC1 <- function(u, z, gg){ z/(z+gg^2*u^2)*dnorm(u) }

    return(2*integrate(integrandC1, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, z = z, gg = gg)$value)
}

###############################################################################
# computation of C3
###############################################################################
C3.fkt <- function(z, gg){
    integrandC3 <- function(u, z, gg){ u^4/(z+gg^2*u^2)*dnorm(u) }

    return(2*integrate(integrandC3, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, z = z, gg = gg)$value)
}


###############################################################################
# computation of alpha.1(x), alpha.3
###############################################################################
.McrgsGeth1 <- function(z.a1, gg, b, a3){
    z <- z.a1[1]; a1.x <- z.a1[2]

    integrandh1 <- function(u, z, gg, b, a1.x, a3){
        u^2/(z+gg^2*u^2)*.McrgsGetw(u, z, gg, b, a1.x, a3)*dnorm(u)
    }
    return(2*integrate(integrandh1, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, z = z, gg = gg, 
                b = b, a1.x = a1.x, a3 = a3)$value)
}
.McrgsGeth2 <- function(z.a1, gg, b, a3){
    z <- z.a1[1]; a1.x <- z.a1[2]

    integrandh2 <- function(u, z, gg, b, a1.x, a3){
        u^4/(z+gg^2*u^2)*.McrgsGetw(u, z, gg, b, a1.x, a3)*dnorm(u)
    }
    return(2*integrate(integrandh2, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, z = z, gg = gg, 
                b = b, a1.x = a1.x, a3 = a3)$value)
}
.McrgsGeth3 <- function(z.a1, gg, b, a3){
    z <- z.a1[1]; a1.x <- z.a1[2]

    integrandh3 <- function(u, z, gg, b, a1.x, a3){
        u^6/(z+gg^2*u^2)*.McrgsGetw(u, z, gg, b, a1.x, a3)*dnorm(u)
    }
    return(2*integrate(integrandh3, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, z = z, gg = gg, 
                b = b, a1.x = a1.x, a3 = a3)$value)
}

.McrgsGeta1a3 <- function(Z.a1, prob, gg, b, a3, C1.x, C3){
    h1.x <- apply(Z.a1, 1, .McrgsGeth1, gg=gg, b=b, a3=a3)
    h2.x <- apply(Z.a1, 1, .McrgsGeth2, gg=gg, b=b, a3=a3)
    h3.x <- apply(Z.a1, 1, .McrgsGeth3, gg=gg, b=b, a3=a3)

    a1.x <- (C1.x - a3*h2.x)/h1.x
    h3 <- sum(prob*h3.x)
    h4 <- sum(prob*a1.x*h2.x)

    return(list(a1.x = a1.x, a3 = (C3 - h4)/h3))
}


###############################################################################
# check of constraints
###############################################################################
.McrgsGetch1 <- function(z.a1, gg, b, a3){
    z <- z.a1[1]; a1.x <- z.a1[2]

    integrand <- function(u, z, gg, b, a1.x, a3){
        (a1.x*u^2 + a3*u^4)/(z+gg^2*u^2)*.McrgsGetw(u, z, gg, b, a1.x, a3)*dnorm(u)
    }
    return(2*integrate(integrand, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, z = z, gg = gg, 
                b = b, a1.x = a1.x, a3 = a3)$value)
}
.McrgsGetch2 <- function(z.a1, gg, b, a3){
    z <- z.a1[1]; a1.x <- z.a1[2]

    integrand <- function(u, z, gg, b, a1.x, a3){
        (a1.x*u^4 + a3*u^6)/(z+gg^2*u^2)*.McrgsGetw(u, z, gg, b, a1.x, a3)*dnorm(u)
    }
    return(2*integrate(integrand, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, z = z, gg = gg, 
                b = b, a1.x = a1.x, a3 = a3)$value)
}


###############################################################################
# computation of b, alpha.1(x), alpha.3
###############################################################################
.McrgsGetba1a3 <- function(r, Z, prob, gg, a1.x, a3, bUp, delta, itmax, check){
    Z.a1 <- matrix(c(Z, a1.x), ncol=2)

    C1.x <- sapply(Z.a1[,1], C1.x.fkt, gg = gg)
    C3.x <- sapply(Z.a1[,1], C3.fkt, gg = gg)
    C3 <- 1 + 1/gg - gg^2*sum(prob*C3.x)

    b <- try(uniroot(.McrgsGetr, lower = gg, upper = bUp, 
                tol = .Machine$double.eps^0.5, Z.a1 = Z.a1, prob = prob,
                r = r, gg = gg, a3 = a3)$root, silent = TRUE)
    if(!is.numeric(b)) b <- gg

    iter <- 0
    repeat{
        iter <- iter + 1
        if(iter > itmax){
            cat("Algorithm did not converge!\n")
            cat("=> increase itmax or try different starting values")
            cat("for 'a1.x' and 'a3'\n")
            return(NA)
        }
        a1.x.old <- a1.x; a3.old <- a3; b.old <- b

        a1.a3 <- try(.McrgsGeta1a3(Z.a1 = Z.a1, prob = prob, gg = gg, 
                        b = b, a3 = a3, C1.x = C1.x, C3 = C3), silent = FALSE)

        if(!is.list(a1.a3)){ 
            a1.x <- a1.x.old + 1e-4
            a3 <- a3.old + 1e-4
            a1.a3 <- list(a1.x = a1.x, a3 = a3)
        }
        a1.x <- a1.a3$a1.x; a3 <- a1.a3$a3
        Z.a1 <- matrix(c(Z, a1.x), byrow=F, ncol=2)
        
        b <- try(uniroot(.McrgsGetr, lower = gg, upper = bUp, 
                    tol = .Machine$double.eps^0.5, Z.a1=Z.a1, r = r, 
                    prob = prob, gg = gg, a3 = a3)$root, silent = TRUE)
        if(!is.numeric(b)) b <- gg
        
        if(max(abs(a1.x.old-a1.x), abs(a3.old-a3), abs(b.old-b)) < delta)
            break
    }

    return(list(b = b, a1.x = a1.x, a3 = a3))
}

###############################################################################
# computation of asymptotic variance
###############################################################################
.McrgsGetvar <- function(z.a1, gg, b, a3, inteps){
    z <- z.a1[1]; a1.x <- z.a1[2]

    integrand <- function(u, z, gg, b, a1.x, a3){
        (a1.x*u + a3*u^3)^2/(z+gg^2*u^2)*.McrgsGetw(u, z, gg, b, a1.x, a3)^2*dnorm(u)
    }
    return(2*integrate(integrand, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, z = z, gg = gg, 
                b = b, a1.x = a1.x, a3 = a3)$value)
}


###############################################################################
# computation of maximum asymptotic MSE
###############################################################################
.McrgsGetmse<- function(gg, Z, prob, r, a1.x, a3, bUp, delta, itmax){
    b.a1.a3 <- .McrgsGetba1a3(r = r, Z = Z, prob = prob, gg = gg, a1.x = a1.x, 
                    a3 = a3, bUp=bUp, delta=delta, itmax=itmax)
    b <- b.a1.a3$b; a1.x <- b.a1.a3$a1.x; a3 <- b.a1.a3$a3
    
    Z.a1 <- matrix(c(Z, a1.x), ncol=2)
    var.x <- apply(Z.a1, 1, .McrgsGetvar, gg = gg, b = b, a3 = a3)
    C1.x <- sapply(Z, C1.x.fkt, gg = gg)
    mse <- sum(prob*var.x) + gg^2*sum(prob*C1.x) + r^2*b^2

    cat("current gamma:\t", gg, "current mse:\t", mse, "\n")
    return(mse)
}


###############################################################################
# computation of optimally robust IC
###############################################################################
rgsOptIC.Mc <- function(r, K, ggLo = 0.5, ggUp = 1.0, a1.x.start, a3.start = 0.25, 
                        bUp = 1000, delta = 1e-5, itmax = 1000, check = FALSE){
    if(!is(K, "DiscreteDistribution"))
        stop("Mc estimators are only implemented for\n",
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
    
    res <- optimize(f = .McrgsGetmse, lower = ggLo, upper = ggUp, 
                maximum = FALSE, tol = .Machine$double.eps^0.3, Z = Z, 
                prob = prob, r = r, a1.x = a1.x, a3 = a3.start, bUp = bUp, 
                delta = delta, itmax = itmax)
    
    gg <- res$minimum
    b.a1.a3 <- .McrgsGetba1a3(r = r, Z = Z, prob = prob, gg = gg, 
                    a1.x = a1.x, a3 = a3.start, bUp = bUp, 
                    delta = delta, itmax = itmax)
    b <- b.a1.a3$b; a1.x <- b.a1.a3$a1.x; a3 <- b.a1.a3$a3

    if(check){
        Z.a1 <- matrix(c(Z, a1.x), ncol=2)
        C1.x <- sapply(Z, C1.x.fkt, gg = gg)
        C3.x <- sapply(Z.a1[,1], C3.fkt, gg=gg)
        C3 <- 1 + 1/gg - gg^2*sum(prob*C3.x)
        kont1 <- apply(Z.a1, 1, .McrgsGetch1, gg=gg, b=b, a3=a3)

        kont21 <- apply(Z.a1, 1, .McrgsGetch2, gg=gg, b=b, a3=a3)
        kont2 <- sum(prob*kont21)

        rvgl <- .McrgsGetr(b = b, Z.a1 = Z.a1, prob = prob, r = r, gg = gg, a3 = a3)
        
        cat("constraint (Mc2):\t", max(abs(kont1-C1.x)), "\n")
        cat("constraint (M4):\t", kont2 - C3, "\n")
        cat("MSE equation:\t", rvgl ,"\n")
    }

    intervall <- getdistrOption("DistrResolution") / 2
    supp.grid <- as.numeric(matrix(
                      rbind(supp - intervall, supp + intervall), nrow = 1))
    a1.grid <- c(as.numeric(matrix(rbind(0, a1.x), nrow = 1)), 0)
    a1fct <- function(x){ stepfun(x = supp.grid, y = a1.grid)(x) }
    
    w <- .McrgsGetw
    fct1 <- function(x){ 
        z <- (A*x[1])^2
        a1fct1 <- a1fct
        wfct <- w 
        a1.x <- a1fct1(x[1])
        w.vct <- wfct(u = x[2], z = z, gg = gg, b = b, a1.x = a1.x, a3 = a3) 
        psi <- ((a1.x*x[2] + a3*x[2]^3)/(z + gg^2*x[2]^2)*w.vct
               + gg^2*x[2]/(z + gg^2*x[2]^2))
        return(A*x[1]*psi)
    }
    body(fct1) <- substitute({ z <- (A*x[1])^2
                               a1fct1 <- a1fct
                               wfct <- w 
                               a1.x <- a1fct1(x[1])
                               w.vct <- wfct(u = x[2], z = z, gg = gg, b = b, a1.x = a1.x, a3 = a3) 
                               psi <- ((a1.x*x[2] + a3*x[2]^3)/(z + gg^2*x[2]^2)*w.vct
                                       + gg^2*x[2]/(z + gg^2*x[2]^2))
                               return(A*x[1]*psi)},
                        list(w = w, b = b, a1fct = a1fct, a3 = a3, A = A, gg = gg))
    fct2 <- function(x){ 
        z <- (A*x[1])^2
        a1fct1 <- a1fct
        wfct <- w 
        a1.x <- a1fct1(x[1])
        w.vct <- wfct(u = x[2], z = z, gg = gg, b = b, a1.x = a1.x, a3 = a3) 
        psi <- ((a1.x*x[2] + a3*x[2]^3)/(z + gg^2*x[2]^2)*w.vct
               + gg^2*x[2]/(z + gg^2*x[2]^2))
        return(gg*(x[2]*psi - 1))
    }
    body(fct2) <- substitute({ z <- (A*x[1])^2
                               a1fct1 <- a1fct
                               wfct <- w 
                               a1.x <- a1fct1(x[1])
                               w.vct <- wfct(u = x[2], z = z, gg = gg, b = b, a1.x = a1.x, a3 = a3) 
                               psi <- ((a1.x*x[2] + a3*x[2]^3)/(z + gg^2*x[2]^2)*w.vct
                                       + gg^2*x[2]/(z + gg^2*x[2]^2))
                               return(gg*(x[2]*psi - 1)) },
                        list(w = w, b = b, a1fct = a1fct, a3 = a3, A = A, gg = gg))
    return(CondIC(name = "conditionally centered IC of Mc type", 
              Curve = EuclRandVarList(RealRandVariable(Map = list(fct1, fct2), 
                                                       Domain = EuclideanSpace(dimension = 2))),
              Risks = list(asMSE = res$objective, asBias = b, trAsCov = res$objective - r^2*b^2), 
              Infos = matrix(c("rgsOptIC.Mc", "optimally robust IC for Mc estimators and 'asMSE'",
                               "rgsOptIC.Mc", paste("where a3 =", round(a3, 3), ", b =", round(b, 3), 
                                                    "and gamma =", round(gg, 3))), 
                           ncol=2, byrow = TRUE, dimnames=list(character(0), c("method", "message"))), 
              CallL2Fam = call("NormLinRegScaleFamily", theta = 0, RegDistr = K, 
                                                         Reg2Mom = as.matrix(Reg2Mom))))
}
