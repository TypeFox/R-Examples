###############################################################################
# computation of maximum asymptotic mse
###############################################################################
.MKrgsGetmse<- function(gg, A, K, r, a1, a3, B, bUp, delta, itmax){
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
    
    cat("current gamma:\t", gg, "current mse:\t", mse, "\n")
    return(mse)
}


###############################################################################
# computation of optimally robust IC
###############################################################################
rgsOptIC.MK <- function(r, K, ggLo = 0.5, ggUp = 1.0, a1.start = -0.25, 
                         a3.start = 0.25, B.start, bUp = 1000, delta = 1e-6, 
                         itmax = 1000, check = FALSE){
    Reg2Mom <- .rgsDesignTest(K = K)
    if(is.logical(Reg2Mom))
        stop("second moment matrix of regressor distribution 'K'", 
             "is (numerically) not positive definite")
    k <- dimension(img(K))
    if(k > 1)
        if(.rgsRegressorCheck(K))
            stop("Regressor is a.e. K concentrated on a conic")

    A <- solve(Reg2Mom)
    if(missing(B.start)) B.start <- A %*% A

    res <- optimize(.MKrgsGetmse, lower = ggLo, upper = ggUp, 
                tol = .Machine$double.eps^0.3, A = A, K = K, r = r, a1 = a1.start, 
                a3 = a3.start, B = B.start, bUp = bUp, delta = delta, itmax = itmax)

    gg <- res$minimum

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
        A.mat <- matrix(A, ncol = k)
        z <- t(x[1:k]) %*% A.mat
        z <- as.vector(z %*% t(z))
        wfct <- w 
        w.vct <- wfct(u = x[k+1], v = v, z = z, gg = gg, b = b, a1 = a1, a3 = a3) 
        psi <- (((a1 + v)*x[k+1] + a3*x[k+1]^3)/(z + gg^2*x[(k+1)]^2)*w.vct
               + gg^2*x[k+1]/(z + gg^2*x[k+1]^2))
        return(A.mat %*% x[1:k]*psi)
    }
    body(fct1) <- substitute({ B.mat <- matrix(B, ncol = k)
                               v <- as.vector(t(x[1:k]) %*% B.mat %*% x[1:k])
                               A.mat <- matrix(A, ncol = k)
                               z <- t(x[1:k]) %*% A.mat
                               z <- as.vector(z %*% t(z))
                               wfct <- w 
                               w.vct <- wfct(u = x[k+1], v = v, z = z, gg = gg, b = b, a1 = a1, a3 = a3) 
                               psi <- (((a1 + v)*x[k+1] + a3*x[k+1]^3)/(z + gg^2*x[(k+1)]^2)*w.vct
                                       + gg^2*x[k+1]/(z + gg^2*x[k+1]^2))
                               return(A.mat %*% x[1:k]*psi)},
                        list(w = w, b = b, a1 = a1, a3 = a3, B = B, A = A, gg = gg, k = k))
    fct2 <- function(x){ 
        B.mat <- matrix(B, ncol = k)
        v <- as.vector(t(x[1:k]) %*% B.mat %*% x[1:k])
        A.mat <- matrix(A, ncol = k)
        z <- t(x[1:k]) %*% A.mat
        z <- as.vector(z %*% t(z))
        wfct <- w 
        w.vct <- wfct(u = x[k+1], v = v, z = z, gg = gg, b = b, a1 = a1, a3 = a3) 
        psi <- (((a1 + v)*x[k+1] + a3*x[k+1]^3)/(z + gg^2*x[(k+1)]^2)*w.vct
               + gg^2*x[k+1]/(z + gg^2*x[k+1]^2))
        return(gg*(x[k+1]*psi - 1))
    }
    body(fct2) <- substitute({ B.mat <- matrix(B, ncol = k)
                               v <- as.vector(t(x[1:k]) %*% B.mat %*% x[1:k])
                               A.mat <- matrix(A, ncol = k)
                               z <- t(x[1:k]) %*% A.mat 
                               z <- as.vector(z %*% t(z))
                               wfct <- w 
                               w.vct <- wfct(u = x[k+1], v = v, z = z, gg = gg, b = b, a1 = a1, a3 = a3) 
                               psi <- (((a1 + v)*x[k+1] + a3*x[k+1]^3)/(z + gg^2*x[(k+1)]^2)*w.vct
                                       + gg^2*x[k+1]/(z + gg^2*x[k+1]^2))
                               return(gg*(x[k+1]*psi - 1)) },
                        list(w = w, b = b, a1 = a1, a3 = a3, B = B, A = A, gg = gg, k = k))
    return(IC(name = "IC of MK type", 
              Curve = EuclRandVarList(EuclRandVariable(Map = list(fct1), 
                                                       Domain = EuclideanSpace(dimension = (trunc(k)+1)),
                                                       dimension = trunc(k)),
                                      RealRandVariable(Map = list(fct2), 
                                                       Domain = EuclideanSpace(dimension = (trunc(k)+1)))),
              Risks = list(asMSE = res$objective, asBias = b, trAsCov = res$objective - r^2*b^2), 
              Infos = matrix(c("rgsOptIC.MK", "optimally robust IC for MK estimators and 'asMSE'",
                               "rgsOptIC.MK", paste("where a1 =", round(a1, 3), ", a3 =", round(a3, 3),
                                                   ", b =", round(b, 3), "and gamma =", round(gg, 3))), 
                           ncol=2, byrow = TRUE, dimnames=list(character(0), c("method", "message"))), 
              CallL2Fam = call("NormLinRegScaleFamily", theta = numeric(k), 
                            RegDistr = K, Reg2Mom = Reg2Mom)))
}
