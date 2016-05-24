###############################################################################
# computation of alpha(x)
###############################################################################
.BMrgsGetalpha <- function(alpha, b.rg.x, b.sc.0.x){
    tau1 <- b.rg.x/alpha
    tau2 <- (1+b.sc.0.x)/b.rg.x
    tau3 <- sqrt((1+b.sc.0.x)/alpha)
    
    if(tau1 <= tau3)
        return(2*(alpha*(pnorm(tau1)-0.5) - b.rg.x*dnorm(tau2) 
               + (1+b.sc.0.x)*pnorm(-tau2)) - 1)
    else
        return(2*(-alpha*tau3*dnorm(tau3) + alpha*(pnorm(tau3)-0.5) 
               + (1+b.sc.0.x)*pnorm(-tau3)) - 1)
}


###############################################################################
# computation of gamma.x
###############################################################################
.BMrgsGetgg <- function(alpha, b.rg.x, b.sc.0.x){
    tau1 <- b.rg.x/alpha
    tau2 <- (1+b.sc.0.x)/b.rg.x
    tau3 <- sqrt((1+b.sc.0.x)/alpha)
    
    if(tau1 <= tau3)
        return(1/(2*(-b.rg.x*dnorm(tau1) + 3*alpha*(pnorm(tau1)-0.5) 
                  - 2*b.rg.x*dnorm(tau2) + (1+b.sc.0.x)*pnorm(-tau2)) - 1))
    else
        return(1/(2*(-3*alpha*tau3*dnorm(tau3) + 3*alpha*(pnorm(tau3)-0.5) 
                     + (1+b.sc.0.x)*pnorm(-tau3)) - 1))
}


###############################################################################
# computation of asymptotic variance
###############################################################################
.BMrgsGetvar <- function(b.rg.x, b.sc.0.x, alpha.x, gg.x, z){
    tau1 <- b.rg.x/alpha.x
    tau2 <- (1+b.sc.0.x)/b.rg.x
    tau3 <- sqrt((1+b.sc.0.x)/alpha.x)
    
    if(tau1 <= tau3){
        erg1 <- 2*(-alpha.x*b.rg.x*dnorm(tau1) + alpha.x^2*(pnorm(tau1)-0.5) 
                   + b.rg.x^2*(pnorm(tau2)-pnorm(tau1)) 
                   + b.rg.x*(1+b.sc.0.x)*dnorm(tau2) 
                   - (1+b.sc.0.x)^2*pnorm(-tau2))
        erg2 <- 2*(-3*alpha.x*b.rg.x*dnorm(tau1) + 3*alpha.x^2*(pnorm(tau1)-0.5) 
                   - b.rg.x*(1+b.sc.0.x)*dnorm(tau2) 
                   + b.rg.x^2*(pnorm(tau2)-pnorm(tau1)) 
                   + (1+b.sc.0.x)^2*pnorm(-tau2))
    }else{
        erg1 <- 2*(-alpha.x^2*tau3*dnorm(tau3) + alpha.x^2*(pnorm(tau3)-0.5) 
                   + (1+b.sc.0.x)^2*(dnorm(tau3)/tau3-pnorm(-tau3)))
        erg2 <- 2*(-alpha.x^2*tau3^3*dnorm(tau3) 
                   + 3*alpha.x^2*(-tau3*dnorm(tau3)+pnorm(tau3)-0.5) 
                   + (1+b.sc.0.x)^2*pnorm(-tau3))
    }
        
    return(z*erg1 + gg.x^2*(erg2 - 1))
}

###############################################################################
# computation of minimum b.rg.x
###############################################################################
.BMrgsGetminbrgx <- function(b.rg.x, b.sc.0.x){
    delta <- (1+b.sc.0.x)/b.rg.x
    
    return(2*b.rg.x*(1/sqrt(2*pi) - dnorm(delta) + delta*pnorm(-delta)) - 1)
}


###############################################################################
# computation of maximum asymptotic MSE
###############################################################################
.BMrgsGetmse <- function(b, r, supp, prob, Z, delta, MAX){
    n <- length(supp)
    b.rg <- b[1]
    b.rg.x <- b.rg/sqrt(Z)
    b.sc.0.x <- b[2:(n+1)]
    
    # constraints
    if(any(b.sc.0.x < 1)) return(MAX)
    b.rg.min.x <- numeric(n)
    for(i in 1:n){
        b.rg.min.x[i] <- uniroot(.BMrgsGetminbrgx, lower = sqrt(pi/2), 
                            upper = 10, tol = .Machine$double.eps^0.5, 
                            b.sc.0.x = b.sc.0.x[i])$root
                            
        if(b.rg.x[i] < b.rg.min.x[i]) return(MAX)
    }

    alpha.x <- numeric(n)
    Var.x <- numeric(n)
    gg.x <- numeric(n)

    for(i in 1:n){
        if(abs(b.rg.x[i]-b.rg.min.x[i]) < 1e-4){
            b.rg.x[i] <- b.rg.min.x[i]
            delta <- (1+b.sc.0.x[i])/b.rg.x[i]
            gg.x[i] <- 1/(4*b.rg.x[i]*(1/sqrt(2*pi)-dnorm(delta) 
                            + 0.5*delta*pnorm(-delta)) - 1)
        
            erg1 <- 2*b.rg.x[i]^2*(pnorm(delta) - 0.5 + delta*dnorm(delta) 
                                   - delta^2*pnorm(-delta))
            erg2 <- 2*b.rg.x[i]^2*(-delta*dnorm(delta) + pnorm(delta) - 0.5 
                                   + delta^2*pnorm(-delta))
            Var.x[i] <- Z[i]*erg1 + gg.x[i]^2*(erg2 - 1)
        }else{
            alpha.x[i] <- uniroot(.BMrgsGetalpha, lower = 0.01, 
                                upper = 1000, tol = delta, b.rg.x = b.rg.x[i], 
                                b.sc.0.x = b.sc.0.x[i])$root
            gg.x[i] <- .BMrgsGetgg(alpha = alpha.x[i], b.rg.x = b.rg.x[i], 
                                b.sc.0.x = b.sc.0.x[i])
            Var.x[i] <- .BMrgsGetvar(b.rg.x = b.rg.x[i], b.sc.0.x = b.sc.0.x[i], 
                                alpha.x = alpha.x[i], gg.x = gg.x[i], z = Z[i])
        }
    }
    b.sc <- max(b.sc.0.x*gg.x)

    return(sum(prob*Var.x) + r^2*(b.rg^2 + b.sc^2))
}


###############################################################################
# computation of optimally robust IC
###############################################################################
rgsOptIC.BM <- function(r, K, b.rg.start = 2.5, b.sc.0.x.start, 
                         delta = 1e-6, MAX = 100, itmax = 1000){
    if(!is(K, "DiscreteDistribution"))
        stop("BM estimators are only implemented for\n",
             "1-dimensional discrete regressor distributions!")
    Reg2Mom <- as.vector(.rgsDesignTest(K = K))
    if(is.logical(Reg2Mom))
        stop("second moment matrix of regressor distribution 'K'", 
             "is (numerically) not positive definite")

    supp <- support(K); prob = d(K)(supp)
    A <- 1/Reg2Mom
    Z <- (A*supp)^2

    n <- length(supp)
    if(missing(b.sc.0.x.start))
        b.sc.0.x.start <- numeric(n) + 2.5

    res <- optim(c(b.rg.start, b.sc.0.x.start), .BMrgsGetmse, 
                method = "Nelder-Mead", 
                control = list(reltol = .Machine$double.eps^0.5, maxit = itmax), 
                r = r, supp = supp, prob = prob, Z = Z, delta = delta, MAX = MAX)

    b.rg <- res$par[1]
    b.rg.x <- b.rg/sqrt(Z)
    b.sc.0.x <- res$par[2:(n+1)]
    
    b.rg.min.x <- numeric(n)
    for(i in 1:n){
        b.rg.min.x[i] <- uniroot(.BMrgsGetminbrgx, lower=sqrt(pi/2), upper=10, tol=1e-8, 
                            maxiter=50, b.sc.0.x=b.sc.0.x[i])$root
    }

    alpha.x <- numeric(n)
    Var.x <- numeric(n)
    gg.x <- numeric(n)
    for(i in 1:n){
        if(abs(b.rg.x[i]-b.rg.min.x[i]) < 1e-4){
            b.rg.x[i] <- b.rg.min.x[i]
            delta <- (1+b.sc.0.x[i])/b.rg.x[i]
            gg.x[i] <- 1/(4*b.rg.x[i]*(1/sqrt(2*pi)-dnorm(delta) 
                            + 0.5*delta*pnorm(-delta)) - 1)
        
            erg1 <- 2*b.rg.x[i]^2*(pnorm(delta) - 0.5 + delta*dnorm(delta) 
                                   - delta^2*pnorm(-delta))
            erg2 <- 2*b.rg.x[i]^2*(-delta*dnorm(delta) + pnorm(delta) - 0.5 
                                   + delta^2*pnorm(-delta))
            Var.x[i] <- Z[i]*erg1 + gg.x[i]^2*(erg2 - 1)
        }else{
            alpha.x[i] <- uniroot(.BMrgsGetalpha, lower = 0.01, 
                                upper = 1000, tol = delta, b.rg.x = b.rg.x[i], 
                                b.sc.0.x = b.sc.0.x[i])$root
            gg.x[i] <- .BMrgsGetgg(alpha = alpha.x[i], b.rg.x = b.rg.x[i], 
                                b.sc.0.x = b.sc.0.x[i])
            Var.x[i] <- .BMrgsGetvar(b.rg.x = b.rg.x[i], b.sc.0.x = b.sc.0.x[i], 
                                alpha.x = alpha.x[i], gg.x = gg.x[i], z = Z[i])
        }
    }
    b.sc <- max(b.sc.0.x*gg.x)
    Var <- sum(prob*Var.x)
    bias <- sqrt(b.rg^2 + b.sc^2)
    mse <- Var + r^2*bias^2
    
    intervall <- getdistrOption("DistrResolution") / 2
    supp.grid <- as.numeric(matrix(
                      rbind(supp - intervall, supp + intervall), nrow = 1))
    a.grid <- c(as.numeric(matrix(rbind(0, alpha.x), nrow = 1)), 0)
    afct <- function(x){ stepfun(x = supp.grid, y = a.grid)(x) }
    bsc.grid <- c(as.numeric(matrix(rbind(0, b.sc.0.x), nrow = 1)), 0)
    bfct <- function(x){ stepfun(x = supp.grid, y = bsc.grid)(x) }
    gg.grid <- c(as.numeric(matrix(rbind(0, gg.x), nrow = 1)), 0)
    gfct <- function(x){ stepfun(x = supp.grid, y = gg.grid)(x) }

    fct1 <- function(x){ 
        z <- (A*x[1])^2
        afct1 <- afct
        alpha.x <- afct1(x[1])
        b.rg.x <- b.rg/sqrt(z)
        bfct1 <- bfct
        b.sc.0.x <- bfct1(x[1])
        psi <- sign(x[2])*min(alpha.x*abs(x[2]), b.rg.x, (1+b.sc.0.x)/abs(x[2]))
        return(A*x[1]*psi)
    }
    body(fct1) <- substitute({ z <- (A*x[1])^2
                               afct1 <- afct
                               alpha.x <- afct1(x[1])
                               b.rg.x <- b.rg/sqrt(z)
                               bfct1 <- bfct
                               b.sc.0.x <- bfct1(x[1])
                               psi <- sign(x[2])*min(alpha.x*abs(x[2]), b.rg.x, (1+b.sc.0.x)/abs(x[2]))
                               return(A*x[1]*psi)},
                        list(b.rg = b.rg, bfct = bfct, afct = afct, A = A))
    fct2 <- function(x){ 
        z <- (A*x[1])^2
        afct1 <- afct
        alpha.x <- afct1(x[1])
        b.rg.x <- b.rg/sqrt(z)
        bfct1 <- bfct
        b.sc.0.x <- bfct1(x[1])
        gfct1 <- gfct
        gg.x <- gfct1(x[1])
        psi <- sign(x[2])*min(alpha.x*abs(x[2]), b.rg.x, (1+b.sc.0.x)/abs(x[2]))
        return(gg.x*(x[2]*psi - 1))
    }
    body(fct2) <- substitute({ z <- (A*x[1])^2
                               afct1 <- afct
                               alpha.x <- afct1(x[1])
                               b.rg.x <- b.rg/sqrt(z)
                               bfct1 <- bfct
                               b.sc.0.x <- bfct1(x[1])
                               gfct1 <- gfct
                               gg.x <- gfct1(x[1])
                               psi <- sign(x[2])*min(alpha.x*abs(x[2]), b.rg.x, (1+b.sc.0.x)/abs(x[2]))
                               return(gg.x*(x[2]*psi - 1)) },
                        list(b.rg = b.rg, bfct = bfct, afct = afct, A = A, gg = gfct))
    return(CondIC(name = "conditionally centered IC of BM type", 
              Curve = EuclRandVarList(RealRandVariable(Map = list(fct1, fct2), 
                                                       Domain = EuclideanSpace(dimension = 2))),
              Risks = list(asMSE = mse, asBias = bias, trAsCov = Var), 
              Infos = matrix(c("rgsOptIC.BM", "optimally robust IC for BM estimators and 'asMSE'",
                               "rgsOptIC.BM", paste("b.rg =", round(b.rg, 3), " and b.sc =", round(b.sc, 3))), 
                           ncol=2, byrow = TRUE, dimnames=list(character(0), c("method", "message"))), 
              CallL2Fam = call("NormLinRegScaleFamily", theta = 0, RegDistr = K, 
                                                         Reg2Mom = as.matrix(Reg2Mom))))
}
