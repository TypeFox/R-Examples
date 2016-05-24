###############################################################################
## Duplication matrix
###############################################################################
.duplicationMatrix <- function(dimn){
    EM <- diag(dimn*(dimn + 1)/2)
    G <- matrix(0, nrow = dimn^2, ncol = dimn*(dimn + 1)/2)
    for(j in 1:dimn)
        for(i in j:dimn){
            G[((j-1)*dimn + i),] <- EM[((j-1)*(dimn-j/2)+i), ]
            G[((i-1)*dimn + j),] <- EM[((j-1)*(dimn-j/2)+i), ]
        }
    return(G)
}


###############################################################################
# Is second moment matrix of regressor positive definite
###############################################################################
.rgsRegressorCheck <- function(K){
    H <- E(K, function(x){ h <- x %*% t(x); h <- h[row(h) >= col(h)]; h %*% t(h)} )
    
    if(identical(all.equal(det(H), 0, tolerance = .Machine$double.eps^0.75), TRUE))
        TRUE
    else
        FALSE
}


###############################################################################
# weight function
###############################################################################
.MrgsGetw <- function(u, v, z, gg, b, a1, a3){
    q.vkt <- sqrt(z + gg^2*u^2)
    g.vkt <- ((a1 + v)*u + a3*u^3)/q.vkt
    b.vkt <- sqrt(b^2 - gg^2*z/q.vkt^2)
    
    g.abs.vkt <- abs(g.vkt)
    ind1 <- (g.abs.vkt < b.vkt)
    
    return(as.vector(ind1 + (1-ind1)*b.vkt/g.abs.vkt))
}


###############################################################################
# computation of b
###############################################################################
.MrgsGetr1 <- function(x, A, gg, b, a1, a3, B){
    z <- A %*% x
    z <- as.vector(t(z) %*% z)
    v <- as.vector(t(x) %*% B %*% x)
    
    integrandr1 <- function(u, v, z, b, a1, a3, gg){
        q.vkt <- sqrt(z + gg^2*u^2)
        g.vkt <- ((a1 + v)*u + a3*u^3)/q.vkt
        b.vkt <- sqrt(b^2 - gg^2*z/q.vkt^2)

        h.vkt <- abs(g.vkt)/b.vkt - 1
            
        return(as.vector((h.vkt > 0)*h.vkt*dnorm(u)))
    }

    return(2*integrate(integrandr1, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, v = v, z = z, b = b, 
                a1 = a1, a3 = a3, gg = gg)$value)
}
.MrgsGetr <- function(b, K, r, A, gg, a1, a3, B){
    r1 <- E(K, .MrgsGetr1, A = A, gg = gg, b = b, a1 = a1, a3 = a3, B = B)    

    return(r-sqrt(r1))
}


###############################################################################
# computation of C1
###############################################################################
.MrgsGetC1 <- function(x, A, gg){
    z <- A %*% x
    z <- as.vector(t(z) %*% z)
    integrandC1 <- function(u, z, gg){ z/(z+gg^2*u^2)*dnorm(u) }

    return(2*integrate(integrandC1, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, z = z, gg = gg)$value)
}


###############################################################################
# computation of C2
###############################################################################
.MrgsGetC2 <- function(x, A, gg){
    z <- A %*% x
    z <- as.vector(t(z) %*% z)
    integrandC2 <- function(u, z, gg){ u^2/(z+gg^2*u^2)*dnorm(u) }
    int <- 2*integrate(integrandC2, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, z = z, gg = gg)$value                
    if(length(x) == 1)
        return(x^2*int)
    else
        return(x %*% t(x)*int)
}


###############################################################################
# computation of C3
###############################################################################
.MrgsGetC3 <- function(x, A, gg){
    z <- A %*% x
    z <- as.vector(t(z) %*% z)
    integrandC3 <- function(u, z, gg){ u^4/(z+gg^2*u^2)*dnorm(u) }
    return(2*integrate(integrandC3, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, z = z, gg = gg)$value)
}


###############################################################################
# computation of alpha1, alpha3 and B
###############################################################################
.MrgsGeth1 <- function(x, A, gg, b, a1, a3, B){
    z <- A %*% x
    z <- as.vector(t(z) %*% z)
    v <- as.vector(t(x) %*% B %*% x)

    integrandh1 <- function(u, v, z, gg, b, a1, a3){
        u^2/(z+gg^2*u^2)*.MrgsGetw(u, v, z, gg, b, a1, a3)*dnorm(u)
    }
    return(2*integrate(integrandh1, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, v = v, z = z, gg = gg, 
                b = b, a1 = a1, a3 = a3)$value)
}
.MrgsGeth11 <- function(x, A, gg, b, a1, a3, B){
    z <- A %*% x
    z <- as.vector(t(z) %*% z)
    v <- as.vector(t(x) %*% B %*% x)

    integrandh1 <- function(u, v, z, gg, b, a1, a3){
        u^2/(z+gg^2*u^2)*.MrgsGetw(u, v, z, gg, b, a1, a3)*dnorm(u)
    }
    return(v*2*integrate(integrandh1, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, v = v, z = z, gg = gg, 
                b = b, a1 = a1, a3 = a3)$value)
}
.MrgsGeth2 <- function(x, A, gg, b, a1, a3, B){
    z <- A %*% x
    z <- as.vector(t(z) %*% z)
    v <- as.vector(t(x) %*% B %*% x)

    integrandh2 <- function(u, v, z, gg, b, a1, a3){
        u^4/(z+gg^2*u^2)*.MrgsGetw(u, v, z, gg, b, a1, a3)*dnorm(u)
    }
    return(2*integrate(integrandh2, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, v = v, z = z, gg = gg, 
                b = b, a1 = a1, a3 = a3)$value)
}
.MrgsGeth21 <- function(x, A, gg, b, a1, a3, B){
    z <- A %*% x
    z <- as.vector(t(z) %*% z)
    v <- as.vector(t(x) %*% B %*% x)

    integrandh2 <- function(u, v, z, gg, b, a1, a3){
        u^4/(z+gg^2*u^2)*.MrgsGetw(u, v, z, gg, b, a1, a3)*dnorm(u)
    }
    return(v*2*integrate(integrandh2, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, v = v, z = z, gg = gg, 
                b = b, a1 = a1, a3 = a3)$value)
}
.MrgsGeth3 <- function(x, A, gg, b, a1, a3, B){
    z <- A %*% x
    z <- as.vector(t(z) %*% z)
    v <- as.vector(t(x) %*% B %*% x)

    integrandh3 <- function(u, v, z, gg, b, a1, a3){
        u^6/(z+gg^2*u^2)*.MrgsGetw(u, v, z, gg, b, a1, a3)*dnorm(u)
    }
    return(2*integrate(integrandh3, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, v = v, z = z, gg = gg, 
                b = b, a1 = a1, a3 = a3)$value)
}
.MrgsGeth4 <- function(x, A, gg, b, a1, a3, B){
    z <- A %*% x
    z <- as.vector(t(z) %*% z)
    v <- as.vector(t(x) %*% B %*% x)

    integrandh3 <- function(u, v, z, gg, b, a1, a3){
        u^2/(z+gg^2*u^2)*.MrgsGetw(u, v, z, gg, b, a1, a3)*dnorm(u)
    }
    int <- 2*integrate(integrandh3, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, v = v, z = z, gg = gg, 
                b = b, a1 = a1, a3 = a3)$value
    if(length(x) == 1)
        return(x^2*int)
    else
        return(x %*% t(x)*int)
}
.MrgsGeth41 <- function(x, A, gg, b, a1, a3, B){
    z <- A %*% x
    z <- as.vector(t(z) %*% z)
    v <- as.vector(t(x) %*% B %*% x)

    integrandh3 <- function(u, v, z, gg, b, a1, a3){
        u^4/(z+gg^2*u^2)*.MrgsGetw(u, v, z, gg, b, a1, a3)*dnorm(u)
    }
    int <- 2*integrate(integrandh3, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, v = v, z = z, gg = gg, 
                b = b, a1 = a1, a3 = a3)$value
    if(length(x) == 1)
        return(x^2*int)
    else
        return(x %*% t(x)*int)
}
.MrgsGeth5 <- function(x, A, gg, b, a1, a3, B){
    z <- A %*% x
    z <- as.vector(t(z) %*% z)
    v <- as.vector(t(x) %*% B %*% x)

    integrandh1 <- function(u, v, z, gg, b, a1, a3){
        u^2/(z+gg^2*u^2)*.MrgsGetw(u, v, z, gg, b, a1, a3)*dnorm(u)
    }
    int <- 2*integrate(integrandh1, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, v = v, z = z, gg = gg, 
                b = b, a1 = a1, a3 = a3)$value
    if(length(x) == 1){
        return(x^4*int)
    }else{
        h <- as.vector(x %*% t(x))
        return(h %*% t(h)*int)
    }
}
.MrgsGeta1a3B <- function(K, A, gg, b, a1, a3, B, C1, C2, C3){
    h1 <- E(K, .MrgsGeth1, A = A, gg = gg, b = b, a1 = a1, a3 = a3, B = B)
    h11 <- E(K, .MrgsGeth11, A = A, gg = gg, b = b, a1 = a1, a3 = a3, B = B)
    h2 <- E(K, .MrgsGeth2, A = A, gg = gg, b = b, a1 = a1, a3 = a3, B = B)
    h21 <- E(K, .MrgsGeth21, A = A, gg = gg, b = b, a1 = a1, a3 = a3, B = B)
    h3 <- E(K, .MrgsGeth3, A = A, gg = gg, b = b, a1 = a1, a3 = a3, B = B)

    a1 <- (C1 - h11 - a3*h2)/h1
    a3 <- (C3 - a1*h2 - h21)/h3

    h4 <- E(K, .MrgsGeth4, A = A, gg = gg, b = b, a1 = a1, a3 = a3, B = B)
    h41 <- E(K, .MrgsGeth41, A = A, gg = gg, b = b, a1 = a1, a3 = a3, B = B)
    D <- C2 - (a1*h4 + a3*h41)
    vech.D <- D[row(D) >= col(D)]
    
    h5 <- E(K, .MrgsGeth5, A = A, gg = gg, b = b, a1 = a1, a3 = a3, B = B)
    
    k <- dimension(img(K))
    Gk <- .duplicationMatrix(dimn = k)
    Hk <- solve(t(Gk) %*% Gk)%*%t(Gk)
    h5 <- Hk %*% h5 %*% Gk

    vech.B <- solve(h5) %*% vech.D
    B <- matrix(0, nrow = k , ncol = k)
    B[row(B) >= col(B)] <- vech.B
    B[row(B) < col(B)] <- B[row(B) > col(B)]
    
    return(list(a1=a1, a3=a3, B=B))
}


###############################################################################
# check of the constraints
###############################################################################
.MrgsGetch1 <- function(x, A, gg, b, a1, a3, B){
    z <- A %*% x
    z <- as.vector(t(z) %*% z)
    v <- as.vector(t(x) %*% B %*% x)

    integrandch1 <- function(u, v, z, gg, b, a1.x, a3){
        ((a1 + v)*u^2 + a3*u^4)/(z+gg^2*u^2)*.MrgsGetw(u, v, z, gg, b, a1, a3)*dnorm(u)
    }
    return(2*integrate(integrandch1, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, v = v, z = z, gg = gg, 
                b = b, a1 = a1, a3 = a3)$value)
}
.MrgsGetch2 <- function(x, A, gg, b, a1, a3, B){
    z <- A %*% x
    z <- as.vector(t(z) %*% z)
    v <- as.vector(t(x) %*% B %*% x)

    integrandch1 <- function(u, v, z, gg, b, a1.x, a3){
        ((a1 + v)*u^2 + a3*u^4)/(z+gg^2*u^2)*.MrgsGetw(u, v, z, gg, b, a1, a3)*dnorm(u)
    }
    int <- 2*integrate(integrandch1, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, v = v, z = z, gg = gg, 
                b = b, a1 = a1, a3 = a3)$value
    if(length(x) == 1)
        return(x^2*int)
    else
        return(x %*% t(x)*int)
}
.MrgsGetch3 <- function(x, A, gg, b, a1, a3, B){
    z <- A %*% x
    z <- as.vector(t(z) %*% z)
    v <- as.vector(t(x) %*% B %*% x)

    integrand <- function(u, v, z, gg, b, a1.x, a3){
        ((a1 + v)*u^4 + a3*u^6)/(z+gg^2*u^2)*.MrgsGetw(u, v, z, gg, b, a1, a3)*dnorm(u)
    }
    return(2*integrate(integrand, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, v = v, z = z, gg = gg, 
                b = b, a1 = a1, a3 = a3)$value)
}


###############################################################################
# computation of b, alpha1, alpha3 and B
###############################################################################
.MrgsGetba1a3B <- function(r, K, A, gg, a1, a3, B, bUp, delta, itmax){
    C1 <- E(K, .MrgsGetC1, A = A, gg = gg)
    C2 <- solve(A) - gg^2*E(K, .MrgsGetC2, A = A, gg = gg)
    C3 <- 1 + 1/gg - gg^2*E(K, .MrgsGetC3, A = A, gg = gg)
    
    b <- try(uniroot(.MrgsGetr, lower = gg, upper = bUp, 
                tol = .Machine$double.eps^0.5, K = K, r = r, A = A, gg = gg,
                a1 = a1, a3 = a3, B = B)$root, silent = TRUE)
    if(!is.numeric(b)) b <- gg
    
    iter <- 0
    repeat{
        iter <- iter + 1
        if(iter > itmax){ 
            stop("Maximum number of iteration reached!\n",
                 "=> increase 'itmax' and try various starting values")
        }
        
        a1.old <- a1; a3.old <- a3; B.old <- B; b.old <- b

        a1.a3.B <- try(.MrgsGeta1a3B(K = K, A = A, gg = gg, b = b, 
                            a1 = a1, a3 = a3, B = B, C1 = C1, C2 = C2, 
                            C3 = C3), silent = TRUE)

        if(!is.list(a1.a3.B)){ 
            a1 <- a1.old + 5e-5; a3 <- a3.old + 5e-5; B <- B.old + 5e-5
            a1.a3.B <- list(a1 = a1, a3 = a3, B = B)
        }
        a1 <- a1.a3.B$a1; a3 <- a1.a3.B$a3; B <- a1.a3.B$B
    
        b <- try(uniroot(.MrgsGetr, lower = gg, upper = bUp, 
                    tol = .Machine$double.eps^0.5, K = K, r = r, A = A, 
                    gg = gg, a1 = a1, a3 = a3, B = B)$root, silent = TRUE)
        if(!is.numeric(b)) b <- gg

        kont1 <- try(E(K, .MrgsGetch1, A = A, gg = gg, b = b, a1 = a1, 
                            a3 = a3, B = B), silent = TRUE)
        if(!is.numeric(kont1)){
            cat("could not determine constraint (M2):\n", kont1, "\n")
            kont1 <- C1 # kont1 <- C1 + 1?
        }

        kont2 <- try(E(K, .MrgsGetch2, A = A, gg = gg, b = b, a1 = a1, 
                            a3 = a3, B = B), silent = TRUE)
        if(!is.numeric(kont2)){
            cat("could not determine constraint (M3):\n", kont2, "\n")
            kont2 <- C2
        }
        
        kont3 <- try(E(K, .MrgsGetch3, A = A, gg = gg, b = b, a1 = a1, 
                            a3 = a3, B = B), silent = TRUE)
        if(!is.numeric(kont3)){ 
            cat("could not determine constraint (M4):\n", kont3, "\n")
            kont3 <- C3
        }
        
#        cat("(M2):\t", kont1 - C1, "(M3):\t", max(abs(kont2 - C2)), "(M4):\t", kont3 - C3, "\n")
        if(max(abs(kont1-C1), abs(kont2-C2), abs(kont3-C3)) < delta) break
    }

    return(list(b=b, a1=a1, a3=a3, B=B))
}


###############################################################################
# computation of asymptotic variance
###############################################################################
.Mrgsgetvar <- function(x, A, gg, b, a1, a3, B){
    z <- A %*% x
    z <- as.vector(t(z) %*% z)
    v <- as.vector(t(x) %*% B %*% x)

    integrand <- function(u, v, z, gg, b, a1, a3){
        ((a1 + v)*u + a3*u^3)^2/(z+gg^2*u^2)*.MrgsGetw(u, v, z, gg, b, a1, a3)^2*dnorm(u)
    }
    return(2*integrate(integrand, lower = 0, upper = Inf, 
                rel.tol = .Machine$double.eps^0.5, v = v, z = z, gg = gg, 
                b = b, a1 = a1, a3 = a3)$value)
}

###############################################################################
# computation of maximum asymptotic mse
###############################################################################
.Mrgsgetmse<- function(vec.A.gg, K, r, a1, a3, B, bUp, delta, MAX, itmax){
    k <- dimension(img(K))
    m <- k*(k+1)/2
    vec.A <- vec.A.gg[1:m]
    gg <- vec.A.gg[m+1]
    cat("current vech(A):\t", vec.A, "current gamma:\t", gg, "\n")
    
    if(gg < 0.5){ 
        cat("gamma < 0.5 => current mse = MAX =", MAX, "\n")
        return(MAX)
    }

    if(m == 1)
        A <- matrix(vec.A)
    else{
        A <- matrix(0, nrow=k, ncol=k)
        A[row(A) >= col(A)] <- vec.A
        A[row(A) < col(A)] <- A[row(A) > col(A)]
    }

    b.a1.a3.B <- try(.MrgsGetba1a3B(r = r, K = K, A = A, 
                        gg = gg, a1 = a1, a3 = a3, B = B, bUp=bUp, 
                        delta = delta, itmax = itmax), silent = TRUE)

    if(!is.list(b.a1.a3.B))
        stop("Algorithm did not converge\n",
             "=> increase 'itmax' or try various starting values")
    
    b <- b.a1.a3.B$b; a1 <- b.a1.a3.B$a1; a3 <- b.a1.a3.B$a3
    B <- b.a1.a3.B$B
        
    var <- E(K, .Mrgsgetvar, A = A, gg = gg, b = b, a1 = a1, a3 = a3, B = B)
    C1 <- E(K, .MrgsGetC1, gg = gg, A = A)
    mse <- var + gg^2*C1 + r^2*b^2
    
    cat("current mse:\t", mse, "\n")
    return(mse)
}


###############################################################################
# computation of optimally robust IC
###############################################################################
rgsOptIC.M <- function(r, K, A.start, gg.start = 0.6, a1.start = -0.25, 
                       a3.start = 0.25, B.start, bUp = 1000, delta = 1e-5, 
                       MAX = 100, itmax = 1000, check = FALSE){
    Reg2Mom <- .rgsDesignTest(K = K)
    if(is.logical(Reg2Mom))
        stop("second moment matrix of regressor distribution 'K'", 
             "is (numerically) not positive definite")
    k <- dimension(img(K))
    if(k > 1)
        if(.rgsRegressorCheck(K))
            stop("Regressor is a.e. K concentrated on a conic")

    if(missing(A.start))
        A <- solve(Reg2Mom)
    else
        A <- A.start

    if(missing(B.start)) B.start <- A %*% A

    vec.A <- as.vector(A[row(A) >= col(A)])
    res <- optim(c(vec.A, gg.start), .Mrgsgetmse, method = "Nelder-Mead", 
                control = list(reltol = 10*delta), K = K, r = r, a1 = a1.start, 
                a3 = a3.start, B = B.start, bUp = bUp, delta = delta, MAX = MAX, 
                itmax = itmax)

    m <- k*(k+1)/2
    vec.A <- res$par[1:m]
    gg <- res$par[m+1]
    
    if(m == 1)
        A <- matrix(vec.A)
    else{
        A <- matrix(0, nrow=k, ncol=k)
        A[row(A) >= col(A)] <- vec.A
        A[row(A) < col(A)] <- A[row(A) > col(A)]
    }

    b.a1.a3.B <- .MrgsGetba1a3B(r = r, K = K, A = A, gg = gg, a1 = a1.start, 
                        a3 = a3.start, B = B.start, bUp = bUp, delta = delta, 
                        itmax = itmax)   
    b <- b.a1.a3.B$b; a1 <- b.a1.a3.B$a1; a3 <- b.a1.a3.B$a3
    B <- b.a1.a3.B$B
    
    if(check){
        C1 <- E(K, .MrgsGetC1, A = A, gg = gg)
        C2 <- solve(A) - gg^2*E(K, .MrgsGetC2, A = A, gg = gg)
        C3 <- 1 + 1/gg - gg^2*E(K, .MrgsGetC3, A = A, gg = gg)

        kont1 <- try(E(K, .MrgsGetch1, A = A, gg = gg, b = b, a1 = a1, 
                            a3 = a3, B = B), silent = TRUE)
        if(is.numeric(kont1))
            cat("constraint (M2):\t", kont1 - C1, "\n")
        else
            cat("could not determine constraint (M2):\n", kont1, "\n")

        kont2 <- try(E(K, .MrgsGetch2, A = A, gg = gg, b = b, a1 = a1, 
                            a3 = a3, B = B), silent = TRUE)
        if(is.numeric(kont2))
            cat("constraint (M3):\t", kont2 - C2, "\n")
        else
            cat("could not determine constraint (M3):\n", kont2, "\n")

        
        kont3 <- try(E(K, .MrgsGetch3, A = A, gg = gg, b = b, a1 = a1, 
                            a3 = a3, B = B), silent = TRUE)
        if(is.numeric(kont3)){ 
            cat("constraint (M4):\t", kont3 - C3, "\n")
        }else
            cat("could not determine constraint (M4):\n", kont3, "\n")

        rvgl <- try(.MrgsGetr(b = b, K = K, r = r, A = A, gg = gg, a1 = a1, 
                            a3 = a3, B = B), silent = TRUE)
        if(is.numeric(rvgl))
            cat("MSE equation:\t", rvgl ,"\n")
        else
            cat("could not determine MSE equation:\n", rvgl ,"\n")
    }

    w <- .MrgsGetw
    fct1 <- function(x){ 
        B.mat <- matrix(B, ncol = k)
        v <- as.vector(t(x[1:k]) %*% B.mat %*% x[1:k])
        A <- matrix(vec.A, ncol = k)
        z <- t(x[1:k]) %*% A 
        z <- as.vector(z %*% t(z))
        wfct <- w 
        w.vct <- wfct(u = x[k+1], v = v, z = z, gg = gg, b = b, a1 = a1, a3 = a3) 
        psi <- (((a1 + v)*x[k+1] + a3*x[k+1]^3)/(z + gg^2*x[(k+1)]^2)*w.vct
               + gg^2*x[k+1]/(z + gg^2*x[k+1]^2))
        return(A %*% x[1:k]*psi)
    }
    body(fct1) <- substitute({ B.mat <- matrix(B, ncol = k)
                               v <- as.vector(t(x[1:k]) %*% B.mat %*% x[1:k])
                               A <- matrix(vec.A, ncol = k)
                               z <- t(x[1:k]) %*% A
                               z <- as.vector(z %*% t(z))
                               wfct <- w 
                               w.vct <- wfct(u = x[k+1], v = v, z = z, gg = gg, b = b, a1 = a1, a3 = a3) 
                               psi <- (((a1 + v)*x[k+1] + a3*x[k+1]^3)/(z + gg^2*x[(k+1)]^2)*w.vct
                                       + gg^2*x[k+1]/(z + gg^2*x[k+1]^2))
                               return(A %*% x[1:k]*psi)},
                        list(w = w, b = b, a1 = a1, a3 = a3, B = B, vec.A = A, gg = gg, k = k))
    fct2 <- function(x){ 
        B.mat <- matrix(B, ncol = k)
        v <- as.vector(t(x[1:k]) %*% B.mat %*% x[1:k])
        A <- matrix(vec.A, ncol = k)
        z <- t(x[1:k]) %*% A 
        z <- as.vector(z %*% t(z))
        wfct <- w 
        w.vct <- wfct(u = x[k+1], v = v, z = z, gg = gg, b = b, a1 = a1, a3 = a3) 
        psi <- (((a1 + v)*x[k+1] + a3*x[k+1]^3)/(z + gg^2*x[(k+1)]^2)*w.vct
               + gg^2*x[k+1]/(z + gg^2*x[k+1]^2))
        return(gg*(x[k+1]*psi - 1))
    }
    body(fct2) <- substitute({ B.mat <- matrix(B, ncol = k)
                               v <- as.vector(t(x[1:k]) %*% B.mat %*% x[1:k])
                               A <- matrix(vec.A, ncol = k)
                               z <- t(x[1:k]) %*% A 
                               z <- as.vector(z %*% t(z))
                               wfct <- w 
                               w.vct <- wfct(u = x[k+1], v = v, z = z, gg = gg, b = b, a1 = a1, a3 = a3) 
                               psi <- (((a1 + v)*x[k+1] + a3*x[k+1]^3)/(z + gg^2*x[(k+1)]^2)*w.vct
                                       + gg^2*x[k+1]/(z + gg^2*x[k+1]^2))
                               return(gg*(x[k+1]*psi - 1)) },
                        list(w = w, b = b, a1 = a1, a3 = a3, B = B, vec.A = A, gg = gg, k = k))
    return(IC(name = "IC of M type", 
              Curve = EuclRandVarList(EuclRandVariable(Map = list(fct1), 
                                                       Domain = EuclideanSpace(dimension = (trunc(k)+1)),
                                                       dimension = trunc(k)),
                                      RealRandVariable(Map = list(fct2), 
                                                       Domain = EuclideanSpace(dimension = (trunc(k)+1)))),
              Risks = list(asMSE = res$value, asBias = b, trAsCov = res$value - r^2*b^2), 
              Infos = matrix(c("rgsOptIC.M", "optimally robust IC for M estimators and 'asMSE'",
                               "rgsOptIC.M", paste("where a1 =", round(a1, 3), ", a3 =", round(a3, 3),
                                                   ", b =", round(b, 3), "and gamma =", round(gg, 3))), 
                           ncol=2, byrow = TRUE, dimnames=list(character(0), c("method", "message"))), 
              CallL2Fam = call("NormLinRegScaleFamily", theta = numeric(k), 
                            RegDistr = K, Reg2Mom = Reg2Mom)))
}
