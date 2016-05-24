###############################################################################
# AR(1) with shift -- *=c, t=alpha=1
###############################################################################

require(distrEx)

###############################################################################
# Computation of optimal clipping bound
###############################################################################
.AR1Sgetc0 <- function(c0, M, H, rad, FFT, eps){
    integrandc0 <- function(x, M, c0, FFT, eps){
        c.H <- c0/pmax(1e-10,abs(x - M))
        c.H1 <- sqrt(c.H/(1+c.H)) 
        if(FFT)
            a.H <- sapply(c.H1, .AR1geta, eps = eps)
        else
            a.H <- .AR1geta.spline(c.H1)
        c1 <- pmax(a.H-c.H,0)
        c2 <- a.H + c.H
        r.H <- abs(x - M)*pmax(((1 + c1 - (a.H-c.H))*exp(-c1) + exp(-c2) + (a.H-c.H) - 1), 0)
    }
    
    if(FFT)
        r1 <- E(H, integrandc0, c0 = c0, M = M, FFT = FFT, eps = eps)
    else
        r1 <- mean(integrandc0(x = H, M = M, c0 = c0, FFT = FFT, eps = eps))

    r <- sqrt(r1/c0)
    
    return(rad-r)
}

##################################################################
# Computation of M
##################################################################
.AR1SgetM1 <- function(x, c0, M, FFT, eps){
    c.H <- c0/pmax(1e-10,abs(x-M))
    c.H1 <- sqrt(c.H/(1+c.H)) 
    if(FFT)
        a.H <- sapply(c.H1, .AR1geta, eps = eps)
    else
        a.H <- .AR1geta.spline(c.H1)
    c1 <- pmax(a.H-c.H,0)
    c2 <- a.H + c.H
    aa.H <- (c1*(2 + c1 - (a.H-c.H)) + 2 - (a.H-c.H))*exp(-c1) - (2+c2)*exp(-c2) - c.H

    return(aa.H)
}

.AR1SgetM2 <- function(x, c0, M, FFT, eps){
    c.H <- c0/pmax(1e-10,abs(x-M))
    c.H1 <- sqrt(c.H/(1+c.H)) 
    if(FFT)
        a.H <- sapply(c.H1, .AR1geta, eps = eps)
    else
        a.H <- .AR1geta.spline(c.H1)
    c1 <- pmax(a.H-c.H, 0)
    c2 <- a.H + c.H
    aa.H <- (c1*(2 + c1 - (a.H-c.H)) + 2 - (a.H-c.H))*exp(-c1) - (2+c2)*exp(-c2) - c.H

    return(x*aa.H)
}

###############################################################################
# Computation of A
###############################################################################
.AR1SgetA <- function(x, c0, M, FFT, eps){
    c.H <- c0/pmax(1e-10,abs(x-M))
    c.H1 <- sqrt(c.H/(1+c.H))
    if(FFT)
        a.H <- sapply(c.H1, .AR1geta, eps = eps)
    else
        a.H <- .AR1geta.spline(c.H1)
    c1 <- pmax(a.H-c.H,0)
    c2 <- a.H + c.H
    aa.H <- (c1*(2 + c1 - (a.H-c.H)) + 2 - (a.H-c.H))*exp(-c1) - (2+c2)*exp(-c2) - c.H

    return(x*(x-M)*aa.H)
}


###############################################################################
# Computation of optimal IC
###############################################################################
AR1SGumbel <- function(rad, phi, shift = 0, M.start, FFT = TRUE, 
                      delta = .Machine$double.eps^0.5, simn = 100000, 
                      offs = 10000, x.start = 0, eps = .Machine$double.eps^0.5, 
                      upper = 1000){  
    if(length(rad) != 1)
        stop("'rad' has to be of length 1")
    if(rad <= 0)
        stop("'rad' has to be in (0,Inf)")
    if(length(phi) != 1)
        stop("'phi' has to be of length 1")
    if((phi <= -1)||(phi >= 1))
        stop("'phi' has to be in (-1,1)")
    Innov <- Gumbel(loc = digamma(1))
    if(missing(M.start)) 
        M <- 0
    else
        M <- M.start

    if(FFT){
        convpow <- ceiling(log(delta)/log(abs(phi)))
        X <- Innov
        for(i in 1:convpow){
            X.alt <- X
            X <- shift + X.alt + (-phi)^i*Innov
        }
        rm(X.alt)
        H <- -X
        rm(X)
    }else{
        x <- x.start
        X <- numeric(simn)
        for(i in 1:(simn + offs)){
            x <- shift - phi*x + r(Innov)(1)
            if(i > offs) X[i-offs] <- x
        }
        H <- -X
        rm(X)
        if(!exists(".AR1geta.spline")){
            c.vct <- seq(from = 1e-6, to = 0.99, length = 10000)
            a.vct <- sapply(c.vct, .AR1geta, eps = eps)
            assign(".AR1geta.spline", splinefun(c.vct, a.vct), env = .GlobalEnv)
        }
    }

    c0 <- 0
    repeat{
        M.alt <- M
        c0.alt <- c0
        
        c0 <- uniroot(.AR1Sgetc0, lower = 1e-4, upper = upper, tol = eps, 
                      M = M, H = H, rad = rad, FFT = FFT, eps = eps)$root

        if(FFT){
            M1 <- E(H, .AR1SgetM1, c0 = c0, M = M, FFT = FFT, eps = eps)
            M2 <- E(H, .AR1SgetM2, c0 = c0, M = M, FFT = FFT, eps = eps)
        }else{
            M1 <- mean(.AR1SgetM1(x = H, c0 = c0, M = M, FFT = FFT, eps = eps))
            M2 <- mean(.AR1SgetM2(x = H, c0 = c0, M = M, FFT = FFT, eps = eps))
        }
        M <- M2/M1

        prec <- max(abs(M.alt-M), abs(c0.alt-c0))
        cat("current precision:\t", prec, "\n")
        if(prec < eps) break
    }

    if(FFT)
        A.H <- 1/E(H, .AR1SgetA, c0 = c0, M = M, FFT = FFT, eps = eps)
    else{
        A.H <- 1/mean(.AR1SgetA(x = H, c0 = c0, M = M, FFT = FFT, eps = eps))
        old <- distroptions("DistrResolution")
        distroptions("DistrResolution", 1e-12)
        H <- DiscreteDistribution(H)
        distroptions("DistrResolution", old)
    }
        
    return(list(A.H = A.H, M = M, b = A.H*c0, H = H))
}
