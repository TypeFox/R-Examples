# Gauﬂ-Legendre abscissas and weights
# cf. for example Numerical Recipies in C (1992), p. 152

#implementation in S:
.GLawOld <- function(n){                                                          
    if(n %% 2 == 1) stop("n has to be an even number")

    m <- (n + 1)/2
    A <- numeric(n)
    W <- numeric(n)
    for(i in 1:floor(m)){
        z <- cos(pi*(i - 0.25)/(n + 0.5))
        repeat{
            p1 <- 1
            p2 <- 0
            for(j in 1:n){
                p3 <- p2
                p2 <- p1
                p1 <- ((2*j - 1)*z*p2 - (j-1)*p3)/j
            }
            pp <- n*(z*p1 - p2)/(z^2 - 1)
            z1 <- z
            z <- z - p1/pp
            if(abs(z - z1) < .Machine$double.eps) break
        }
        A[i] <- -z
        A[n + 1 - i] <- z
        W[i] <- 2/((1 - z^2)*pp^2)
        W[n + 1 - i] <- W[i]
    }
    cbind(A, W)
}

## new: interface to C P.R. 01-04-06

.GLaw <- function(n){
    if(n %% 2 == 1) stop("n has to be an even number")

    A <- numeric(n)
    W <- numeric(n)

#    mm<-dyn.load("G:/rtest/GLaw.dll")
    erg<-.C("gauleg",n = as.integer(n),eps = as.double(.Machine$double.eps),
             A = as.double(A),W = as.double(W), PACKAGE = "distrEx") 
              ### PACKAGE ARGUMENT added P.R. 270507
#    dyn.unload("G:/rtest/GLaw.dll")
#
# P.R. 20140810: .Call interface instead of .C interface
#
#   erg0 <- .Call("Gauleg", n, eps, PACKAGE="distrEx")
#   erg <- matrix(erg0,n,2); colnames(erg) <- c("A","W")
#
    cbind(A=erg$A, W=erg$W)         
}


GLIntegrate <- function(f, lower, upper, order = 500, ...){
    if(order %in% c(100, 500, 1000))
        AW <- getFromNamespace(paste(".AW", as.character(order), 
                                     sep = "."), ns = "distrEx")
    else
        AW <- .GLaw(order)

    # transformation to [lower, upper]
    xl <- (upper - lower)/2
    W <- xl*AW[,2]
    A <- xl*AW[,1] + (lower + upper)/2

    res <- W*f(A, ...)
    sum(res)
}

distrExIntegrate <- function(f, lower, upper, subdivisions = 100, 
                             rel.tol = .Machine$double.eps^0.25, 
                             abs.tol = rel.tol, stop.on.error = TRUE, 
                             distr, order = .distrExOptions$GLIntegrateOrder, 
                             ...){
    res <- try(integrate(f, lower = lower, upper = upper, rel.tol = rel.tol, 
                  abs.tol = abs.tol, stop.on.error = stop.on.error, ...)$value, 
                  silent = TRUE)

    # if integrate fails => Gauﬂ-Legendre integration
    if(!is.numeric(res)){
        Zi <- 1
        if(lower >= upper){
            lo <- lower
            lower <- upper
            upper <- lo
            Zi <- -1
        }
        if(!is.finite(lower))
            if(missing(distr)) stop(res)
        else{
            value <- c(...)
            if(is(distr, "UnivariateCondDistribution"))
                lower <- q(distr)(.distrExOptions$GLIntegrateTruncQuantile, 
                                   cond = value$cond)
            else
                lower <- q(distr)(.distrExOptions$GLIntegrateTruncQuantile)
        }
        if(!is.finite(upper))
            if(missing(distr)) stop(res)
        else{
            value <- c(...)
            if(is(distr, "UnivariateCondDistribution")){
                if("lower.tail" %in% names(formals(distr@q)))
                  upper <- q(distr)(.distrExOptions$GLIntegrateTruncQuantile, 
                                      cond = value$cond, lower.tail = FALSE)  
                else    
                  upper <- q(distr)(1 - .distrExOptions$GLIntegrateTruncQuantile, 
                                      cond = value$cond)                                      
            }else{
                if("lower.tail" %in% names(formals(distr@q)))
                  upper <- q(distr)(.distrExOptions$GLIntegrateTruncQuantile, 
                                     lower.tail = FALSE)  
                else    
                  upper <- q(distr)(1 - .distrExOptions$GLIntegrateTruncQuantile)
            }
        }
        res <- Zi*GLIntegrate(f = f, lower = lower, upper = upper, 
                              order = order, ...)
    }
    
    return(res)
}
