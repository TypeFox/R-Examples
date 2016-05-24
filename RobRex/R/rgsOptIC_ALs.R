###############################################################################
# weight function for regression part
###############################################################################
.ALsrgsGetwrg <- function(u, b.rg, A.rg.x){
    h.vkt <- sqrt(as.vector(A.rg.x %*% A.rg.x))*abs(u)
    ind1 <- (h.vkt < b.rg)

    return(as.vector(ind1 + (1-ind1)*b.rg/h.vkt))
}

###############################################################################
# computation of b.rg
###############################################################################
.ALsrgsGetrrg1 <- function(x, b.rg, A.rg){
    A.rg.x <- as.vector(A.rg %*% x)
    integrandr <- function(u, b.rg, A.rg.x){
        h1.vkt <- as.vector(A.rg.x %*% A.rg.x)*u^2
        h.vkt <- sqrt(h1.vkt)/b.rg - 1
        Ind <- (h.vkt > 0)
    
        return(Ind*h.vkt*dnorm(u))
    }

    return(2*integrate(integrandr, lower = 0, upper = Inf, 
                    rel.tol = .Machine$double.eps^0.5, b.rg = b.rg,
                    A.rg.x = A.rg.x)$value)
}
.ALsrgsGetrrg <- function(b.rg, K, r, A.rg){
    r1 <- E(K, .ALsrgsGetrrg1, b.rg = b.rg, A.rg = A.rg)

    return(r-sqrt(r1))
}


###############################################################################
# computation of a.sc, c.sc
###############################################################################
.ALsGetasc <- function(a.sc, c.sc){
    c1 <- sqrt(pmax(a.sc - c.sc, 0))
    c2 <- sqrt(c.sc + a.sc)
    fa <- (c1*dnorm(c1) - c2*dnorm(c2) + (1-a.sc-c.sc)*pnorm(c2)
           - (1-a.sc+c.sc)*pnorm(c1) + 1.5*c.sc)

    return(fa)
}
.ALsGetcsc <- function(c.sc, r){
    c.sc <- c.sc^2/(1-c.sc^2)
    a <- uniroot(.ALsGetasc, lower = 0, upper = 1,
                 tol = .Machine$double.eps^0.5, c.sc=c.sc)$root

    c1 <- sqrt(pmax(a-c.sc, 0))
    c2 <- sqrt(c.sc+a)

    r1 <- 2*((1-a-c.sc)*(1-pnorm(c2)) + c2*dnorm(c2)
           + (c.sc-a+1)*(1/2-pnorm(c1)) + c1*dnorm(c1))/c.sc

    return(r-sqrt(r1))
}

###############################################################################
# computation of A.rg
###############################################################################
.ALsrgsGetArg1 <- function(x, b.rg, A.rg, z.sc, A.sc){
    A.rg.x <- as.vector(A.rg %*% x)
    integrandArg <- function(u, b.rg, A.rg.x){
        u^2*.ALsrgsGetwrg(u = u, b.rg = b.rg, A.rg.x = A.rg.x)*dnorm(u)
    }
    int <- 2*integrate(integrandArg, lower = 0, upper = Inf, 
                    rel.tol = .Machine$double.eps^0.5, b.rg = b.rg,
                    A.rg.x = A.rg.x)$value
    if(length(x) == 1)
        return(x^2*int)
    else
        return(x%*%t(x)*int)
}
.ALsrgsGetArg <- function(K, b.rg, A.rg){
    A.rg1 <- E(K, .ALsrgsGetArg1, b.rg = b.rg, A.rg = A.rg)

    return(solve(A.rg1))
}


###############################################################################
# check MSE equation for eta.sc
###############################################################################
.ALsGetcheck <- function(c.sc, a.sc, r){
    c1 <- sqrt(pmax(a.sc-c.sc, 0))
    c2 <- sqrt(c.sc+a.sc)

    r1 <- 2*((1-a.sc-c.sc)*(1-pnorm(c2)) + c2*dnorm(c2)
           + (c.sc-a.sc+1)*(1/2-pnorm(c1)) + c1*dnorm(c1))/c.sc

    return(r-sqrt(r1))
}

###############################################################################
# computation of optimally robust IC
###############################################################################
rgsOptIC.ALs <- function(r, K, A.rg.start, b.rg.Up = 1000, delta=1e-06,
                         itmax = 50, check = FALSE){
    Reg2Mom <- .rgsDesignTest(K = K)
    if(is.logical(Reg2Mom))
        stop("second moment matrix of regressor distribution 'K'\n", 
             "is (numerically) not positive definite")

    if(missing(A.rg.start))
        A.rg <- solve(Reg2Mom)
    else
        A.rg <- A.rg.start

    b.rg <- uniroot(.ALsrgsGetrrg, lower = 1e-4, upper = b.rg.Up,
            tol = .Machine$double.eps^0.5, K = K,
            r = r, A.rg = A.rg)$root
    
    iter <- 0
    repeat{
        iter <- iter + 1
        if(iter > itmax){
            cat("Algorithm did not converge!\n")
            cat("=> increase itmax or try different starting values ")
            cat("for 'A.rg'\n")
            return(NA)
        }
        A.rg.old <- A.rg; b.rg.old <- b.rg
        
        A.rg <- .ALsrgsGetArg(K = K, b.rg = b.rg, A.rg = A.rg)
        b.rg <- uniroot(.ALsrgsGetrrg, lower = 1e-4, upper = b.rg.Up,
                    tol = .Machine$double.eps^0.5, K = K, r = r,
                    A.rg = A.rg)$root

        if(max(abs(A.rg.old-A.rg), abs(b.rg.old-b.rg)) < delta)
            break
    }

    c.sc <- uniroot(.ALsGetcsc, lower = 0.01, upper = 0.991,
                tol = .Machine$double.eps^0.5, r = r)$root
    c.sc <- c.sc^2/(1-c.sc^2)

    a.sc <- uniroot(.ALsGetasc, lower = 0, upper = 1,
                tol = .Machine$double.eps^0.5, c.sc = c.sc)$root
    z.sc <- a.sc - 1

    c1 <- sqrt(pmax(a.sc-c.sc,0))
    c2 <- sqrt(c.sc + a.sc)

    aa <- dnorm(c2)*(-c2^3-3*c2+2*a.sc*c2)-dnorm(c1)*(-c1^3-3*c1+2*a.sc*c1)
    aa <- aa + (pnorm(c2)-pnorm(c1))*(a.sc^2+3-2*a.sc)
    aa <- aa + c.sc*(1-a.sc)*(1.5-pnorm(c2)-pnorm(c1))
    aa <- aa + c.sc*(c2*dnorm(c2)+c1*dnorm(c1))

    A.sc <- 1/(2*aa)

    if(check){
        kont.rg1 <- E(K, .ALsrgsGetArg1, b.rg = b.rg, A.rg = A.rg)
        kont.rg <- A.rg%*%kont.rg1
    
        kont.sc1 <- .ALsGetasc(a.sc = a.sc, c.sc = c.sc)

        rvgl.rg <- .ALsrgsGetrrg(b.rg = b.rg, K = K, r = r, A.rg = A.rg)
        rvgl.sc <- .ALsGetcheck(c.sc = c.sc, a.sc = a.sc, r = r)

        cat("Fisher consistency of eta.rg:\n")
        print(kont.rg-diag(dimension(img(K))))
        cat("MSE equation for eta.rg:\t", rvgl.rg ,"\n")
        cat("Centering of eta.sc:\t", kont.sc1, "\n")
        cat("MSE equation for eta.sc:\t", rvgl.sc, "\n")
    }

    k <- dimension(img(K))
    vec.A <- as.vector(A.rg)
    w <- .ALsrgsGetwrg
    fct1 <- function(x){
        A.rg <- matrix(vec.A, ncol = k)
        A.rg.x <- as.vector(A.rg %*% x)
        wfct <- w
        ww <- wfct(u = x[k+1], b.rg = b.rg, A.rg.x = A.rg.x)
        return(A.rg %*% x[1:k]*x[k+1]*ww)
    }
    body(fct1) <- substitute({ A.rg <- matrix(vec.A, ncol = k)
                               A.rg.x <- as.vector(A.rg %*% x[1:k])
                               wfct <- w
                               ww <- wfct(u = x[k+1], b.rg = b.rg, A.rg.x = A.rg.x)
                               return(A.rg %*% x[1:k]*x[k+1]*ww)},
                        list(w = w, b.rg = b.rg, vec.A = vec.A, k = k))
    fct2 <- function(x){
        return(A.sc*(x[k+1]^2 - a.sc)*pmin(1, c.sc/abs(x[k+1]^2 - a.sc)))
    }
    body(fct2) <- substitute({ return(A.sc*(x[k+1]^2 - a.sc)*pmin(1, c.sc/abs(x[k+1]^2 - a.sc))) },
                        list(A.sc = A.sc, a.sc = a.sc, c.sc = c.sc, k = k))
    asMSE <- sum(diag(A.rg)) + A.sc
    asBias <- sqrt(b.rg^2 + A.sc^2*c.sc^2)
    trAsCov <- asMSE - r^2*asBias^2

    return(IC(name = "IC of ALs type",
              Curve = EuclRandVarList(EuclRandVariable(Map = list(fct1),
                                                       Domain = EuclideanSpace(dimension = (trunc(k)+1)),
                                                       dimension = trunc(k)),
                                      RealRandVariable(Map = list(fct2),
                                                       Domain = EuclideanSpace(dimension = (trunc(k)+1)))),
              Risks = list(asMSE = asMSE, asBias = asBias, trAsCov = trAsCov),
              Infos = matrix(c("rgsOptIC.ALs", "optimally robust IC for ALs estimators and 'asMSE'",
                               "rgsOptIC.ALs", paste("where b.rg =", round(b.rg, 3), ", A.sc =", round(A.sc, 3), 
                                                     ", z.sc =", round(z.sc, 3), "and c.sc =", round(c.sc, 3))),
                           ncol=2, byrow = TRUE, dimnames=list(character(0), c("method", "message"))),
              CallL2Fam = call("NormLinRegScaleFamily", theta = numeric(k),
                            RegDistr = K, Reg2Mom = Reg2Mom)))
}
