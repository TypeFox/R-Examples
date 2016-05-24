## Test disabled
quit()

## test calculation of matrices for the classes
## rlmerPredD and rlmerPred_...
require(robustlmm)
.calcE.D.re <- robustlmm:::.calcE.D.re
.calcE.psi_bbt <- robustlmm:::.calcE.psi_bbt
.calcE.psi_bpsi_bt <- robustlmm:::.calcE.psi_bpsi_bt
setTheta <- robustlmm:::setTheta

calcMatrices <- function(object, numpoints=13) {
    X <- object@pp$X
    Zt <- object@pp$Zt
    rho.resp <- object@rho.e
    rho.re <- object@rho.b
    method <- object@method
    
    p <- object@pp$p
    q <- object@pp$q
    n <- object@pp$n
    
    ## D_e and D_b matrices, as well as
    ## Epsi2_e, Epsi2_b, Epsi_bbt and Epsi_bpsi_bt
    D_e <- Diagonal(x=rep(object@rho.e@EDpsi(), n))
    Epsi2_e <- object@rho.e@Epsi2()
    tmp <- numeric(0)
    t2 <- t3 <- list()
    for (l in seq_along(object@dim)) {
        ld <- object@dim[l]
        lq <- object@q[l]
        lr <- object@rho.b[[l]]
        tmp <- c(tmp, rep(.calcE.D.re(ld, lr), lq))
        t2 <- c(t2, rep(list(.calcE.psi_bbt(lr, ld)), lq/ld))
        t3 <- c(t3, rep(list(.calcE.psi_bpsi_bt(lr, ld)), lq/ld))
    }
    D_b <- Diagonal(x=tmp)
    Epsi_bbt <- bdiag(t2)
    Epsi_bpsi_bt <- bdiag(t3)
    Epsi2_b <- diag(Epsi_bpsi_bt)
    Lambda_b <- solve(D_b) * object@rho.e@EDpsi()
    
    ## Matrices we need for calculating the Jacobian
    DX <- D_e %*% X
    ZtD <- Zt %*% D_e
    XtDX <- crossprod(X, DX)
    ZtDX <- ZtD %*% X
    ZtDZ <- tcrossprod(ZtD, Zt)
    
    ## Initialize Jacobian Matrix (complete for given theta)
    J0 <- Matrix(0, p + q, p + q)
    J0[1:p, 1:p] <- XtDX
    
    ## CXt = C X\tr = solve(X\tr D.resp X, X\tr)
    CXt <- solve(XtDX, t(X))
    ## H = X CXt
    H <- X %*% CXt
    
    ## lfrac = la = \lambda_e / \lambda_b
    laD.re <- Lambda_b %*% D_b
    
    I <- Matrix(diag(n))
    ## P = I - D.resp H
    P <- I - D_e %*% H
    
    ZtPDZ <- tcrossprod(Zt %*% P, ZtD)
    CXtD <- CXt %*% D_e 
    
    ## Now we can assume that there are the matrices
    ## CXt, H, I, P, ZtPDZ, CXtD are in the environment
    
    ## get index of non-zero U
    idx <- !object@pp$zeroB
    
    if (any(idx)) {
        ## La = \Lambda_\theta
        Lat <- t(La <- object@pp$U_b)
        LtZt <- Lat %*% Zt
        
        ## Ms = M_\theta^* = solve(L\tr Z\tr D_resp P Z L + lD)
        ## Mst = Ms since D_resp is diagonal
        ## MsLtZt = solve(L\tr Z\tr P D_resp Z L + lD), L\tr Z\tr)
        ## LZMs = MsLtZt\tr
        Ms <- solve(Lat %*% ZtPDZ %*% La +
                    laD.re)
        MsLtZt <- Ms %*% LtZt
        
        ## Q = Q_\theta = CXt D.resp Z L Ms
        Qt <- tcrossprod(MsLtZt, CXtD)
        ## S = S_\theta = X Q
        St <- tcrossprod(Qt, X)
        ## T = Z L S\tr
        T <- crossprod(LtZt, St)
        
        ## K = K_\theta = Q\tr X\tr - Ms Z\tr
        K <- St - MsLtZt
        ## L = L_\theta = object@pp$lfrac * Ms
        L <- Ms %*% Lambda_b
        
        ## A = H - T - T\tr P + Z L Ms L\tr Z\tr
        A <- H - T - crossprod(T, P) + crossprod(MsLtZt, LtZt)
        ## B = lambda_e / lambda_b *(S - Z L Ms) 
        B <- t(K) %*% Lambda_b
        
        ## Complete Jacobian
        J <- J0
        J[p+(1:q), 1:p] <- t(J[1:p, p+(1:q)] <- crossprod(ZtDX, object@pp$U_b))
        t.idx <- (p+(1:q))[idx]
        J[t.idx, t.idx] <- crossprod(La[idx, idx], ZtDZ[idx, idx] %*% La[idx, idx]) + laD.re[idx, idx]
    } else {
        Ms <- solve(laD.re)
        A <- H
        ## Fixing Dimnames slot
        A@Dimnames <- list(A@Dimnames[[1]], NULL)
        B <- Matrix(0, object@pp$n, object@pp$q)
        K <- Matrix(0, object@pp$q, object@pp$n)
        L <- Ms %*% Lambda_b
        J <- J0
    }
    
    return(list(A = A, B = B, K = K, L = L, J = J,
                D_e = D_e, D_b = D_b, Lambda_b = Lambda_b,
                Epsi2_e = Epsi2_e, Epsi2_b = Epsi2_b,
                Epsi_bbt = Epsi_bbt, Epsi_bpsi_bt = Epsi_bpsi_bt,
                laD.re = laD.re))
}

cmp <- function(rfm) {
    las <- calcMatrices(rfm)
    res <- c(D_e=all.equal(las$D_e, rfm@pp$D_e),
             D_b=all.equal(las$D_b, rfm@pp$D_b),
             Lambda_b=all.equal(las$Lambda_b, rfm@pp$Lambda_b),
             A=all.equal(las$A, rfm@pp$A),
             B=all.equal(las$B, rfm@pp$B()),
             K=all.equal(las$K, rfm@pp$K()),
             L=all.equal(las$L, rfm@pp$L),
             J=all.equal(las$J, rfm@pp$J()),
             Epsi2_e=all.equal(las$Epsi2_e, rfm@pp$Epsi2_e),
             Epsi2_b=all.equal(las$Epsi2_b, rfm@pp$Epsi2_b),
             Epsi_bbt=all.equal(las$Epsi_bbt, rfm@pp$Epsi_bbt),
             Epsi_bpsi_bt=all.equal(las$Epsi_bpsi_bt, rfm@pp$Epsi_bpsi_bt))
    if (!all(sapply(res, isTRUE))) {
        for (i in seq_along(res))
            cat("Test i=", i, "(", names(res)[i], "), ", res[i], "\n")
        stop("Test failed")
    } else {
        cat("passed\n")
    }
}

test <- function(formula, data, ...) {
    rfm <- rlmer(formula, data, method="DASvar", ..., doFit=FALSE)
    rfm@pp$updateMatrices()
    theta0 <- theta(rfm)
    len <- length(theta0)
    ## test for original theta
    cat("Testing original theta... ")
    cmp(rfm)
    ## test for zero theta (all zero)
    cat("Testing all zero theta... ")
    setTheta(rfm, rep(0, len), fit.effects=FALSE)
    cmp(rfm)
    ## test one theta zero sequentially
    if (len > 1) {
        for (l in 1:len) {
            cat("Testing theta_", l, "=0... ", sep="")
            lth <- theta0
            lth[l] <- 0
            setTheta(rfm, lth, fit.effects=FALSE)
            cmp(rfm)
        }
    }
}

ttest <- function(formula, data) {
    test(formula, data)
    test(formula, data, rho.e = smoothPsi)
    test(formula, data, rho.b = smoothPsi)
    test(formula, data, rho.e = smoothPsi, rho.b = chgDefaults(smoothPsi, k=1, s=10))
}

cat("----- Dyestuff -----\n")
ttest(Yield ~ (1|Batch), Dyestuff)
cat("----- Sleepstudy -----\n")
ttest(Reaction ~ Days + (Days|Subject), sleepstudy)
cat("----- Sleepstudy2 -----\n")
sleepstudy2 <- within(sleepstudy, Group <- letters[1:4])
ttest(Reaction ~ Days + (Days|Subject) + (1|Group), sleepstudy2)

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
