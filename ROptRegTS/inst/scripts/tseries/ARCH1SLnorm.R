###############################################################################
# ARCH(1) with lognormal innovations -- *=c, t=alpha=1
###############################################################################

require(distrEx)

###############################################################################
# Computation of optimal clipping bound
###############################################################################
.ARCH1Sgetc0 <- function(c0, H, M, rad){
    c.H <- c0/pmax(1e-10,abs(H-M))    

    return(rad^2*c0 - 2*mean(abs(H-M)*(dnorm(c.H) - c.H*pnorm(-c.H))))
}


###############################################################################
# Computation of optimal IC
###############################################################################
ARCH1SLnorm <- function(rad, alpha, scale = 1, shape = 1, H, M.start, 
                        eps = .Machine$double.eps^0.5, simn = 1e5, offs = 1e3, 
                        x.start = 0, upper = 1000){  
    if(length(rad) != 1)
        stop("'rad' has to be of length 1")
    if(rad <= 0)
        stop("'rad' has to be in (0,Inf)")
    if(length(alpha) != 1)
        stop("'alpha' has to be of length 1")

    fstat <- function(x, shape){log(abs(exp(shape*x)-exp(shape^2/2)))*dnorm(x)}
    int1 <- integrate(fstat, lower = -Inf, upper = shape/2, shape = shape)$value
    int2 <- distrExIntegrate(fstat, lower = shape/2+1e-10, upper = 10, shape = shape)
    MAX <- exp(-2*(int1 + int2))
    if((alpha <= 0)||(alpha >= MAX))
        stop("'alpha' has to be in (0,", MAX, ")")
    
    if(missing(M.start)) 
        M <- 0
    else
        M <- M.start

    if(missing(H)){
        x <- x.start
        X <- numeric(simn)
        for(i in 1:(simn+offs)){
            x <- scale*sqrt(1+alpha*x^2)*(rlnorm(1, sdlog = shape) - exp(0.5*shape^2))
            if(i > offs) X[i-offs] <- x
        }
        X1 <- c(X[simn],X[1:(simn-1)])
        rm(X)
        H <- X1^2/(2*(1 + alpha*X1^2))
        rm(X1)
    }

    c0 <- 0
    repeat{
        M.alt <- M
        c0.alt <- c0

        c0 <- uniroot(.ARCH1Sgetc0, lower = 1e-6, upper = upper, tol = eps, 
                      M = M, H = H, rad = rad)$root

        c.H <- c0/pmax(1e-10,abs(H-M))
        aa.H <- 2*pnorm(c.H)-1        
        M <- sum(H*aa.H)/sum(aa.H)
    
        prec <- max(abs(M.alt-M), abs(c0.alt-c0))
        cat("current precision:\t", prec, "\n")
        if(prec < eps) break
    }

    A.H <- 1/(mean(H^2*aa.H) - M*mean(H*aa.H))
    old <- distroptions("DistrResolution")
    distroptions("DistrResolution", 1e-12)
    H <- DiscreteDistribution(H)
    distroptions("DistrResolution", old)
    
    return(list(A.H = A.H, M = M, b = A.H*c0/shape, H = H))
}
