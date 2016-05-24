EuclCondition <- function(dimension){
    new("EuclCondition", Range = EuclideanSpace(dimension = dimension))
}
                                                   
## access methods
setMethod("Range", "EuclCondition", function(object) object@Range)

LMParameter <- function(theta = 0, intercept = 0, scale = 1){
    if(any(!is.finite(theta)))
        stop("inifinite or missing values in 'theta'")
    if(length(intercept) != 1)
        stop("'intercept' has to be of length 1")
    if(!is.finite(intercept))
        stop("inifinite or missing value in 'intercept'")
    if(length(scale) != 1)
        stop("'scale' has to be of length 1")
    if(!is.finite(scale))
        stop("inifinite or missing value in 'scale'")

    LMP <- new("LMParameter")
    LMP@theta <- theta
    LMP@intercept <- intercept
    LMP@scale <- scale
    
    return(LMP)
}

## generating function
LMCondDistribution <- function(Error = Norm(), theta = 0, intercept = 0, 
                                scale = 1){
    if(!is(Error, "AbscontDistribution"))
        stop("distribution of 'Error' has to be of class 'AbscontDistribution'")
    param <- LMParameter(theta = theta, intercept = intercept, scale = scale)
    lth <- length(theta)
    cond <- EuclCondition(dimension = floor(lth))
    rfct <- r(Error)
    rfun <- function(n, cond, ...){}
    body(rfun) <- substitute({ if(length(cond) != lth) 
                                    stop("'cond' has wrong dimension")
                               r <- rfct
                               intercept + cond %*% theta + scale*r(n, ...) },
                        list(rfct = rfct, lth = lth, 
                             intercept = intercept, theta = theta, 
                             scale = scale))
    
    dfct <- d(Error)
    dfun <- function(x, cond, log = FALSE, ...){}
    body(dfun) <- substitute({ if(length(cond) != lth) 
                                    stop("'cond' has wrong dimension")
                               d <- dfct
                               if ("log" %in% names(formals(d)))
                                    d0 <- d((x - intercept - 
                                             as.vector(cond %*% theta))/scale, 
                                            log = log)
                               else {d0 <- d((x - intercept - 
                                             as.vector(cond %*% theta))/scale)
                                     if (log) d0 <- log(d0)}
                                            
                               if (log) 
                                   d0 <- d0 - log(scale)
                               else
                                   d0 <- d0 / scale    
                               return(d0)
                               },
                        list(dfct = dfct, lth = lth, 
                             intercept = intercept, theta = theta, 
                             scale = scale))

    pfct <- p(Error) 
    pfun <- function(q, cond, lower.tail = TRUE, log.p = FALSE, ...){}
    body(pfun) <- substitute({ if(length(cond) != lth) 
                                    stop("'cond' has wrong dimension")
                               p <- pfct
                               argList <- alist( (q - intercept - 
                                              as.vector(cond %*% theta))/scale ) 
                               if ("lower.tail" %in% names(formals(p)))
                                  argList <- c(argList, lower.tail = lower.tail)
                               if ("log.p" %in% names(formals(p)))
                                  argList <- c(argList, log.p = log.p)
                               dots <- alist(...)
                               if(length(dots)) argList <- c(argList, dots)
                               p0 <- do.call(p, argList)     
                               if (!("lower.tail" %in% names(formals(p))))
                                   if (!lower.tail) p0 <- 1 - p0  
                               if (!("log.p" %in% names(formals(p))))                                 
                                   if (log.p) p0 <- log( p0 )
                              return(p0)    
                               },
                        list(pfct = pfct, lth = lth, 
                             intercept = intercept, theta = theta, 
                             scale = scale))

    qfct <- q(Error)
    qfun <- function(p, cond, lower.tail = TRUE, log.p = FALSE, ...){}
    body(qfun) <- substitute({ if(length(cond) != lth) 
                                    stop("'cond' has wrong dimension")
                               q <- qfct
                               argList <- alist( p ) 
                               if ("lower.tail" %in% names(formals(q)))
                                   argList <- c(argList, lower.tail = lower.tail)
                               else if (log.p) p <- exp( p )   
                               if ("log.p" %in% names(formals(q)))
                                   argList <- c(argList, log.p = log.p)
                               else if (!lower.tail) p <- 1 - p
                               dots <- alist(...)
                               if(length(dots)) argList <- c(argList, dots)
                               scale*do.call(q, argList)  +  
                                       intercept + as.vector(cond %*% theta) }, 
                        list(qfct = qfct, lth = lth, 
                             intercept = intercept, theta = theta, 
                             scale = scale))

    CD1 <- new("AbscontCondDistribution")
    CD1@r <- rfun 
    CD1@d <- dfun
    CD1@p <- pfun
    CD1@q <- qfun
    CD1@param <- param
    CD1@cond <- cond
    CD1@img <- Reals()
    CD1@.withSim <- Error@.withSim
    CD1@.withArith <- TRUE
    CD1@.logExact <- FALSE
    CD1@.lowerExact <- Error@.lowerExact
    return(CD1)
}
