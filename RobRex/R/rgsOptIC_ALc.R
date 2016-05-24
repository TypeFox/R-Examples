###############################################################################
# weight function
###############################################################################
.ALcrgsGetw <- function(u, b, A.rg.x, z.sc, A.sc){
    h.vkt <- sqrt((t(A.rg.x)%*%A.rg.x)*u^2 + A.sc^2*(u^2-z.sc)^2)
    ind1 <- (h.vkt < b)
    
    return(as.vector(ind1 + (1-ind1)*b/h.vkt))
}

###############################################################################
# computation of b
###############################################################################
.ALcrgsGetr1 <- function(z.sc.x, b, A.rg, A.sc){
    z.sc <- z.sc.x[1]; x <- z.sc.x[2:length(z.sc.x)]
    
    A.rg.x <- A.rg %*% x
    integrand <- function(u, b, A.rg.x, z.sc, A.sc){
        h1.vkt <- (t(A.rg.x)%*%A.rg.x)*u^2 + A.sc^2*(u^2-z.sc)^2
        h.vkt <- sqrt(h1.vkt)/b - 1
    
        return((h.vkt > 0)*h.vkt*dnorm(u))
    }
    return(2*integrate(integrand, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, b = b, A.rg.x = A.rg.x, 
                z.sc = z.sc, A.sc = A.sc)$value)
}
.ALcrgsGetr <- function(b, supp, prob, r, A.rg, z.sc, A.sc){
    z.sc.x <- cbind(z.sc, supp)
    r1 <- apply(z.sc.x, 1, .ALcrgsGetr1, b = b, A.rg = A.rg, A.sc = A.sc)
    
    return(r-sqrt(sum(prob*r1)))
}


###############################################################################
# computation of A.rg, A.sc and z.sc.x
###############################################################################
.ALcrgsGetArg <- function(z.sc.x, b, A.rg, A.sc){
    z.sc <- z.sc.x[1]; x <- z.sc.x[2:length(z.sc.x)]
    
    A.rg.x <- A.rg %*% x
    integrandArg <- function(u, b, A.rg.x, z.sc, A.sc){
        u^2*.ALcrgsGetw(u, b, A.rg.x, z.sc, A.sc)*dnorm(u)
    }
    return(2*integrate(integrandArg, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, b = b, A.rg.x = A.rg.x, 
                z.sc = z.sc, A.sc = A.sc)$value)
}
.ALcrgsGetzh1 <- function(z.sc.x, b, A.rg, A.sc){
    z.sc <- z.sc.x[1]; x <- z.sc.x[2:length(z.sc.x)]
    
    A.rg.x <- A.rg %*% x
    integrandzh1 <- function(u, b, A.rg.x, z.sc, A.sc){
        u^2*.ALcrgsGetw(u, b, A.rg.x, z.sc, A.sc)*dnorm(u)
    }
    return(integrate(integrandzh1, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, b = b, A.rg.x = A.rg.x, 
                z.sc = z.sc, A.sc = A.sc)$value)
}
.ALcrgsGetzh2 <- function(z.sc.x, b, A.rg, A.sc){
    z.sc <- z.sc.x[1]; x <- z.sc.x[2:length(z.sc.x)]
    
    A.rg.x <- A.rg %*% x
    integrandzh2 <- function(u, b, A.rg.x, z.sc, A.sc){
        .ALcrgsGetw(u, b, A.rg.x, z.sc, A.sc)*dnorm(u)
    }
    return(integrate(integrandzh2, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, b = b, A.rg.x = A.rg.x, 
                z.sc = z.sc, A.sc = A.sc)$value)
}
.ALcrgsGetAsc <- function(z.sc.x, b, A.rg, A.sc){
    z.sc <- z.sc.x[1]; x <- z.sc.x[2:length(z.sc.x)]
    
    A.rg.x <- A.rg %*% x
    integrandAsc <- function(u, b, A.rg.x, z.sc, A.sc){
        (u^4 - z.sc*u^2)*.ALcrgsGetw(u, b, A.rg.x, z.sc, A.sc)*dnorm(u)
    }
    return(2*integrate(integrandAsc, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, b = b, A.rg.x = A.rg.x, 
                z.sc = z.sc, A.sc = A.sc)$value)
}
.ALcrgsGetAz <- function(supp, prob, b, A.rg, z.sc, A.sc){
    z.sc.x <- cbind(z.sc, supp)
    
    A.rg1 <- apply(z.sc.x, 1, .ALcrgsGetArg, b = b, A.rg = A.rg, A.sc = A.sc)
    summe <- matrix(0, ncol=ncol(supp), nrow=ncol(supp))
    for(i in 1:nrow(supp)){
        summe <- summe + prob[i]*supp[i,]%*%t(supp[i,])*A.rg1[i]
    }
    A.rg <- solve(summe)
    
    A.sc1 <- apply(z.sc.x, 1, .ALcrgsGetAsc, b = b, A.rg = A.rg, A.sc = A.sc)
    A.sc <- 1/sum(prob*A.sc1)

    z.sc1 <- apply(z.sc.x, 1, .ALcrgsGetzh1, b = b, A.rg = A.rg, A.sc = A.sc)
    z.sc2 <- apply(z.sc.x, 1, .ALcrgsGetzh2, b = b, A.rg = A.rg, A.sc = A.sc)
    z.sc <- z.sc1/z.sc2

    return(list(A.rg = A.rg, z.sc = z.sc, A.sc = A.sc))
}

###############################################################################
# check of centering of eta.sc
###############################################################################
.ALcrgsGetcheck <- function(z.sc.x, b, A.rg, z.sc, A.sc){
    z.sc <- z.sc.x[1]; x <- z.sc.x[2:length(z.sc.x)]
    
    A.rg.x <- A.rg %*% x
    integrand <- function(u, b, A.rg.x, z.sc, A.sc){
        (u^2 - z.sc)*.ALcrgsGetw(u, b, A.rg.x, z.sc, A.sc)*dnorm(u)
    }
    return(2*integrate(integrand, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, b = b, A.rg.x = A.rg.x, 
                z.sc = z.sc, A.sc = A.sc)$value)
}


###############################################################################
# computation of optimally robust IC
###############################################################################
rgsOptIC.ALc <- function(r, K, theta, scale = 1, A.rg.start, a.sc.start, 
                         A.sc.start = 0.5, bUp = 1000, delta = 1e-6, itmax = 50, 
                         check = FALSE){
    if(!is(K, "DiscreteDistribution") & !is(K, "DiscreteMVDistribution"))
        stop("ALc estimators are only implemented for discrete distributions")
    Reg2Mom <- .rgsDesignTest(K = K)
    if(is.logical(Reg2Mom))
        stop("second moment matrix of regressor distribution 'K'", 
             "is (numerically) not positive definite")

    supp <- as.matrix(support(K)); prob <- d(K)(supp)
    k <- ncol(supp)
    if(!missing(theta))
        if(length(theta) != k)
            stop("'theta' has wrong dimension")

    A.sc <- A.sc.start 
    if(missing(a.sc.start))
        z.sc <- numeric(nrow(supp)) + 1
    else
        z.sc <- a.sc.start/A.sc + 1
    if(missing(A.rg.start))
        A.rg <- solve(Reg2Mom)
    else
        A.rg <- A.rg.start


    b <- uniroot(.ALcrgsGetr, lower = 1e-4, upper = bUp, 
            tol = .Machine$double.eps^0.5, supp = supp, prob = prob, 
            r = r, A.rg = A.rg, z.sc = z.sc, A.sc = A.sc)$root
    
    iter <- 0
    repeat{
        iter <- iter + 1
        if(iter > itmax){
            cat("Algorithm did not converge!\n")
            cat("=> increase itmax or try different starting values")
            cat("for 'A.rg', 'a.sc' and 'A.sc'\n")
            return(NA)
        }
        A.rg.old <- A.rg; z.sc.old <- z.sc; A.sc.old <- A.sc; b.old <- b
        
        A.a <- .ALcrgsGetAz(supp, prob, b = b, A.rg = A.rg, 
                        z.sc = z.sc, A.sc = A.sc)
        A.rg <- A.a$A.rg; z.sc <- A.a$z.sc; A.sc <- A.a$A.sc

        b <- uniroot(.ALcrgsGetr, lower = 1e-4, upper = bUp, 
                tol = .Machine$double.eps^0.5, supp = supp, prob = prob, 
                r = r, A.rg = A.rg, z.sc = z.sc, A.sc = A.sc)$root
        
        if(max(abs(A.rg.old-A.rg), abs(z.sc.old-z.sc), abs(A.sc.old-A.sc), abs(b.old-b))<delta)
            break
    }

    if(check){
        z.sc.x <- cbind(z.sc, supp)
        kont.rg1 <- apply(z.sc.x, 1, .ALcrgsGetArg, b = b, A.rg = A.rg, A.sc = A.sc)
        summe <- matrix(0, ncol=ncol(supp), nrow=ncol(supp))
        for(i in 1:nrow(supp)){
            summe <- summe + prob[i]*supp[i,]%*%t(supp[i,])*kont.rg1[i]
        }
        kont.rg <- A.rg%*%summe
    
        kont.sc21 <- apply(z.sc.x, 1, .ALcrgsGetAsc, b = b, A.rg = A.rg, A.sc = A.sc)
        kont.sc2 <- A.sc*sum(prob*kont.sc21)

        kont.sc11 <- apply(z.sc.x, 1, .ALcrgsGetcheck, b = b, A.rg = A.rg, A.sc = A.sc)
        kont.sc1 <- A.sc*kont.sc11
    
        rvgl <- .ALcrgsGetr(b = b, supp = supp, prob = prob, r = r, A.rg = A.rg, 
                        z.sc = z.sc, A.sc = A.sc)
        
        cat("Conditional centering of eta.sc:\t", max(abs(kont.sc1)), "\n")
        cat("Fisher consistency of eta.rg:\n") 
        print(kont.rg-diag(ncol(supp)))
        cat("Fisher consistency of eta.sc:\t", kont.sc2-1, "\n")
        cat("MSE equation:\t", rvgl ,"\n")
    }

    a.sc <- A.sc*(z.sc - 1)*scale
    A <- matrix(0, ncol = ncol(A.rg)+1, nrow = nrow(A.rg)+1)
    A[1:nrow(A.rg),1:ncol(A.rg)] <- scale^2*A.rg
    A[(nrow(A.rg)+1),(ncol(A.rg)+1)] <- scale^2*A.sc
    mse <- sum(diag(A.rg)) + A.sc
    
    b <- scale*b
    fct1 <- function(x){ numeric(k) }
    body(fct1) <- substitute({ numeric(k) }, list(k = k))
    if(is(K, "DiscreteMVDistribution")){
        fct2 <- function(x){ 
            if(liesInSupport(K, x[1:k])){
                ind <- colSums(apply(supp, 1, "==", x[1:k])) == k
                return(a.sc[ind])
            }else{
                return(NA)
            }
        }
        body(fct2) <- substitute({ if(liesInSupport(K, x[1:k])){
                                       ind <- colSums(apply(supp, 1, "==", x[1:k])) == k
                                       return(a.sc[ind])
                                   }else{
                                       return(NA)
                                   }}, list(k = k))
    }
    if(is(K, "DiscreteDistribution")){
        fct2 <- function(x){ 
            if(liesInSupport(K, x[1])){
                ind <- (round(x[1], 8) == round(supp, 8))
                return(a.sc[ind])
            }else{
                return(NA)
            }
        }
    }
    
    a.sc.fct <- EuclRandVarList(EuclRandVariable(Map = list(fct1), 
                                             Domain = EuclideanSpace(dimension = trunc(k) + 1),
                                             dimension = trunc(k)),
                                EuclRandVariable(Map = list(fct2), 
                                             Domain = EuclideanSpace(dimension = trunc(k) + 1),
                                             dimension = 1))

    if(missing(theta)) theta <- numeric(k)    

    return(generateIC(neighbor = Av1CondContNeighborhood(radius = r), 
                L2Fam = NormLinRegScaleFamily(theta = theta, scale = scale,
                            RegDistr = K, Reg2Mom = Reg2Mom), 
                res = list(A = A, a = a.sc.fct, b = b, 
                           d = NULL, risk = list(asMSE = mse, asBias = b, trAsCov = mse - r^2*b^2), 
                           info = c("rgsOptIC.ALc", "optimally robust IC for ALc estimators and 'asMSE'"))))
}
