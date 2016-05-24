###############################################################################
## Class: PrognCondition                   
## condition in case of the simple prognose model 
## y = x + u
## i.e., conditioning by realizations of y
###############################################################################
setClass("PrognCondition", 
            representation(range = "EuclideanSpace"), 
            prototype(name = "condition in case of a linear regression model",
                      range = new("EuclideanSpace")),
            contains = "Condition")
PrognCondition <- function(range = EuclideanSpace()){
    new("PrognCondition", range = range)
}

setMethod("show", "PrognCondition",
    function(object){
        cat(gettextf("Name:\t%s\n", object@name))
        cat(gettextf("Range:\t%s with dimension %s\n", object@range@name, 
            object@range@dimension))
    })

## generating function
PrognCondDistribution <- function(Regr = Norm(), Error = Norm(),
             rel.tol= getdistrExOption("ErelativeTolerance"), 
             lowerTruncQuantile = getdistrExOption("ElowerTruncQuantile"), 
             upperTruncQuantile = getdistrExOption("EupperTruncQuantile"), 
             IQR.fac = getdistrExOption("IQR.fac")){
    if(!is(Error, "AbscontDistribution"))
        stop("Error has to be of type 'AbscontDistribution'")
    if(!is(Regr, "AbscontDistribution"))
        stop("Regr has to be of type 'AbscontDistribution'")
    param <- NULL
    cond <- PrognCondition(range = new("Reals"))
    rfun <- function(n, cond, ...){}
    body(rfun) <- substitute({ r <- rfun; cond - r(n,...)}, 
                             list(rfun = r(Error)))
                        
    dxfun <- d(Regr)
    dufun <- d(Error)
    qxfun <- q(Regr)

    Ib <- .getIntbounds(Error, low=-Inf, upp=Inf, lowerTruncQuantile, 
                       upperTruncQuantile, IQR.fac)
    low <- Ib["low"]
    upp <- Ib["upp"]


    eps <-  getdistrOption("TruncQuantile")
    dfun <- function(x, cond, log = FALSE){}
    body(dfun) <- substitute({ dx <- dxfun; du <- dufun; qx <- qxfun
                               dy <- function(cond){ 
                                  ix <- distrExIntegrate(f = function(x, cond){ 
                                                       dx <- dxfun 
                                                       du <- dufun
                                                       dx(x)*du(cond-x) }, 
                                        lower = low, 
                                        upper = upp, rel.tol=rel.tol,
                                        cond = cond)
                                  return(ix)
                                  }
                               if ("log" %in% names(formals(dx)) && log)
                                    d0 <- dx(x, log = TRUE) + 
                                          du(cond-x, log = TRUE) - log(dy)
                               else if (log) 
                                    d0 <- log(dx(x)) + log(du(cond-x)) - log(dy)
                               else d0 <- dx(x)*du(cond-x)/dy(cond)
                               return(d0)
                              },
                        list(dxfun = dxfun, dufun = dufun, qxfun = qxfun))
 
    pfun <- function(q, cond, lower.tail = TRUE, log.p = FALSE){} 

    body(pfun) <- substitute({ d <- dfun; qx <- qxfun
                               if (lower.tail)
                               p0 <- distrExIntegrate(f = d, lower = low, 
                                      upper = q, rel.tol = rel.tol, 
                                      cond = cond)
                               else 
                               p0 <- distrExIntegrate(f = d, lower = q, 
                                      upper = upp, rel.tol = rel.tol,
                                      cond = cond)
                               if (log.p) p0 <- log(p0)
                               },
                        list(dfun = dfun, qxfun = qxfun))

    qufun <- q(Error)
    qfun <- function(p, cond, lower.tail = TRUE, log.p = FALSE){}
    body(qfun) <- substitute({ qu <- qufun 
                               if (log.p) p <- exp(p)
                               p1 <- if (lower.tail) 1-p else p
                               q0 <- if( "lower.tail" %in% 
                                              names(formals(qu)))
                                   cond - qu(p, lower.tail = !lower.tail) else
                                   cond - qu(p1)
                                               },
                        list(qufun = qufun))
    
    return(new("AbscontCondDistribution", r = rfun, d = dfun, p = pfun, q = qfun, 
            param = param, cond = cond, .withSim = Error@.withSim | Regr@.withSim,
            .withArith = TRUE, .logExact = FALSE, 
            .lowerExact = Error@.lowerExact & Regr@.lowerExact))
}
