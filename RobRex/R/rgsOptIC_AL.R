###############################################################################
# Is second moment matrix of regressor positive definite
###############################################################################
.rgsDesignTest <- function(K){
    Ex2 <- E(K, function(x){ x %*% t(x) })
    
    if(!all(eigen(Ex2)$values > 100*.Machine$double.eps)) FALSE else as.matrix(Ex2)
}

###############################################################################
# weight function
###############################################################################
.ALrgsGetw <- function(u, b, A.rg.x, z.sc, A.sc){
    h.vkt <- sqrt(as.vector(A.rg.x%*%A.rg.x)*u^2 + A.sc^2*(u^2-z.sc)^2)
    ind1 <- (h.vkt < b)
    
    return(as.vector(ind1 + (1-ind1)*b/h.vkt))
}

###############################################################################
# computation of b
###############################################################################
.ALrgsGetr1 <- function(x, b, A.rg, z.sc, A.sc){
    A.rg.x <- as.vector(A.rg %*% x)
    integrandr <- function(u, b, A.rg.x, z.sc, A.sc){
        h1.vkt <- as.vector(A.rg.x%*%A.rg.x)*u^2 + A.sc^2*(u^2-z.sc)^2
        h.vkt <- sqrt(h1.vkt)/b - 1
        Ind <- (h.vkt > 0)
    
        return(Ind*h.vkt*dnorm(u))
    }

    return(2*integrate(integrandr, lower = 0, upper = Inf, 
                    rel.tol = .Machine$double.eps^0.5, b = b, 
                    A.rg.x = A.rg.x, z.sc = z.sc, A.sc = A.sc)$value)
}


###############################################################################
# computation of b
###############################################################################
.ALrgsGetr <- function(b, K, r, A.rg, z.sc, A.sc){
    r1 <- E(K, .ALrgsGetr1, b = b, A.rg = A.rg, z.sc = z.sc, A.sc = A.sc)    

    return(r-sqrt(r1))
}


###############################################################################
# computation of A.rg, z.sc, A.sc
###############################################################################
.ALrgsGetArg <- function(x, b, A.rg, z.sc, A.sc){
    A.rg.x <- as.vector(A.rg %*% x)
    integrandArg <- function(u, b, A.rg.x, z.sc, A.sc){
        u^2*.ALrgsGetw(u, b, A.rg.x, z.sc, A.sc)*dnorm(u)
    }
    int <- 2*integrate(integrandArg, lower = 0, upper = Inf, 
                    rel.tol = .Machine$double.eps^0.5, b = b, A.rg.x = A.rg.x, 
                    z.sc = z.sc, A.sc = A.sc)$value
    if(length(x) == 1)
        return(x^2*int)
    else
        return(x%*%t(x)*int)
}
.ALrgsGetzsc1 <- function(x, b, A.rg, z.sc, A.sc){
    A.rg.x <- as.vector(A.rg %*% x)
    integrandz1 <- function(u, b, A.rg.x, z.sc, A.sc){
        u^2*.ALrgsGetw(u, b, A.rg.x, z.sc, A.sc)*dnorm(u)
    }
    return(integrate(integrandz1, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, b = b, A.rg.x = A.rg.x, 
                z.sc = z.sc, A.sc = A.sc)$value)
}
.ALrgsGetzsc2 <- function(x, b, A.rg, z.sc, A.sc){
    A.rg.x <- as.vector(A.rg %*% x)
    integrandz2 <- function(u, b, A.rg.x, z.sc, A.sc){
        .ALrgsGetw(u, b, A.rg.x, z.sc, A.sc)*dnorm(u)
    }
    return(integrate(integrandz2, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, b = b, A.rg.x = A.rg.x, 
                z.sc = z.sc, A.sc = A.sc)$value)
}
.ALrgsGetAsc <- function(x, b, A.rg, z.sc, A.sc){
    A.rg.x <- as.vector(A.rg %*% x)
    integrandAsc <- function(u, b, A.rg.x, z.sc, A.sc){
        (u^4 - z.sc*u^2)*.ALrgsGetw(u, b, A.rg.x, z.sc, A.sc)*dnorm(u)
    }
    return(2*integrate(integrandAsc, lower = 0, upper = Inf, 
                    rel.tol = .Machine$double.eps^0.5, b = b, A.rg.x = A.rg.x, 
                    z.sc = z.sc, A.sc = A.sc)$value)
}
.ALrgsGetAz <- function(K, b, A.rg, z.sc, A.sc){
    A.rg1 <- E(K, .ALrgsGetArg, b = b, A.rg = A.rg, 
                z.sc = z.sc, A.sc = A.sc)
    A.rg <- solve(A.rg1)

    A.sc1 <- E(K, .ALrgsGetAsc, b = b, A.rg = A.rg, 
                z.sc = z.sc, A.sc = A.sc)    
    A.sc <- 1/A.sc1

    z.sc1 <- E(K, .ALrgsGetzsc1, b = b, A.rg = A.rg, 
                z.sc = z.sc, A.sc = A.sc)
    z.sc2 <- E(K, .ALrgsGetzsc2, b = b, A.rg = A.rg, 
                z.sc = z.sc, A.sc = A.sc)
    z.sc <- z.sc1/z.sc2

    return(list(A.rg=A.rg, z.sc=z.sc, A.sc=A.sc))
}


###############################################################################
# check centering of eta.sc
###############################################################################
.ALrgsGetcheck <- function(x, b, A.rg, z.sc, A.sc){
    A.rg.x <- as.vector(A.rg %*% x)
    integrandcheck <- function(u, b, A.rg.x, z.sc, A.sc){
        (u^2 - z.sc)*.ALrgsGetw(u, b, A.rg.x, z.sc, A.sc)*dnorm(u)
    }
    return(2*integrate(integrandcheck, lower = 0, upper = Inf, 
                    rel.tol = .Machine$double.eps^0.5, b = b, A.rg.x = A.rg.x, 
                    z.sc = z.sc, A.sc = A.sc)$value)
}


###############################################################################
# computation of optimally robust IC
###############################################################################
rgsOptIC.AL <- function(r, K, theta, scale = 1, A.rg.start, a.sc.start = 0, 
                        A.sc.start = 0.5, bUp = 1000, delta=1e-06, itmax = 50, 
                        check = FALSE){
    k <- dimension(img(K))
    if(!missing(theta))
        if(length(theta) != k)
            stop("'theta' has wrong dimension")
    Reg2Mom <- .rgsDesignTest(K = K)
    if(is.logical(Reg2Mom))
        stop("second moment matrix of regressor distribution 'K'\n", 
             "is (numerically) not positive definite")

    A.sc <- A.sc.start; z.sc <- a.sc.start/A.sc + 1
    if(missing(A.rg.start))
        A.rg <- solve(Reg2Mom)
    else
        A.rg <- A.rg.start

    b <- uniroot(.ALrgsGetr, lower = 1e-4, upper = bUp, 
            tol = .Machine$double.eps^0.5, K = K,
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
        
        A.z <- .ALrgsGetAz(K = K, b = b, A.rg = A.rg, z.sc = z.sc, A.sc = A.sc)
        A.rg <- A.z$A.rg; z.sc <- A.z$z.sc; A.sc <- A.z$A.sc
        
        b <- uniroot(.ALrgsGetr, lower = 1e-4, upper = bUp, 
                tol = .Machine$double.eps^0.5, K = K, r = r, 
                A.rg = A.rg, z.sc = z.sc, A.sc = A.sc)$root
        if(max(abs(A.rg.old-A.rg), abs(z.sc.old-z.sc), abs(A.sc.old-A.sc), abs(b.old-b)) < delta)
            break
    }

    if(check){
        kont.rg1 <- E(K, .ALrgsGetArg, b = b, A.rg = A.rg, 
                        z.sc = z.sc, A.sc = A.sc)
        kont.rg <- A.rg%*%kont.rg1
    
        kont.sc21 <- E(K, .ALrgsGetAsc, b = b, A.rg = A.rg, 
                        z.sc = z.sc, A.sc = A.sc)    
        kont.sc2 <- A.sc*kont.sc21

        kont.sc11 <- E(K, .ALrgsGetcheck, b = b, A.rg = A.rg, 
                        z.sc = z.sc, A.sc = A.sc)    
        kont.sc1 <- A.sc*kont.sc11
    
        rvgl <- .ALrgsGetr(b = b, K = K, r = r, A.rg = A.rg, 
                        z.sc = z.sc, A.sc = A.sc)

        cat("centering of eta.sc:\t", kont.sc1, "\n")
        cat("Fisher consistency of eta.rg:\n")
        print(kont.rg-diag(dimension(img(K))))
        cat("Fisher consistency of eta.sc:\t", kont.sc2-1, "\n")
        cat("MSE equation:\t", rvgl ,"\n")
    }

    a <- c(numeric(nrow(A.rg)), A.sc*(z.sc - 1)*scale)
    A <- matrix(0, ncol = ncol(A.rg)+1, nrow = nrow(A.rg)+1)
    A[1:nrow(A.rg),1:ncol(A.rg)] <- scale^2*A.rg
    A[(nrow(A.rg)+1),(ncol(A.rg)+1)] <- scale^2*A.sc
    b <- scale*b
    if(missing(theta)) theta <- numeric(k)
    mse <- sum(diag(A.rg)) + A.sc
    
    return(generateIC(neighbor = ContNeighborhood(radius = r), 
                L2Fam = NormLinRegScaleFamily(theta = theta, scale = scale, 
                            RegDistr = K, Reg2Mom = Reg2Mom), 
                res = list(A = A, a = a, b = b, d = NULL, 
                           risk = list(asMSE = mse, asBias = b, trAsCov = mse - r^2*b^2), 
                           info = c("rgsOptIC.AL", "optimally robust IC for AL estimators and 'asMSE'"))))
}
