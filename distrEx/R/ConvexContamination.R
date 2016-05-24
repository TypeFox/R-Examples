###############################################################################
## Convex Contamination
## (1-size)*e1 + size*e2, size~Binom(1, size)
## e1: ideal distribution
## e2: contaminating distribution
## size: amout of contamination (gross errors)
###############################################################################
setMethod("ConvexContamination", signature(e1 = "AbscontDistribution",
                                           e2 = "AbscontDistribution",
                                           size = "numeric"),
    function(e1, e2, size){
        TruncQuantile <- getdistrOption("TruncQuantile")
        if(length(size) != 1)
            stop("length of 'size' has to be 1")
        if((size < 0)|(size > 1))
            stop("'size' has to be in [0,1]")
        rfun <- function(n){}
        body(rfun) <- substitute({ r1 <- r1fun; r2 <- r2fun
                                   ind <- rbinom(n, prob=size, size=1)
                                   (1-ind)*r1(n) + ind*r2(n)},
                            list(size = size, r1fun = r(e1), r2fun = r(e2)))

        dfun <- function(x, log = FALSE){} 
        body(dfun) <- substitute({ d1 <- d1fun; d2 <- d2fun
                                   d0 <- (1-size)*d1(x) + size*d2(x)
                                   if (log) d0 <- log(d0)
                                   return(d0)
                                   },
                            list(size = size, d1fun = d(e1), d2fun = d(e2)))

        pfun <- function(q, lower.tail = TRUE, log.p = FALSE){} 
        body(pfun) <- substitute({ 
                         p1 <- function(x){
                            if ("lower.tail" %in% names(formals(p1fun)))
                                 p1fun(x, lower.tail)
                            else {p0 <- p1fun(x)
                                  if(!lower.tail) p0 <- 1-p0
                                  return(p0)} }     
                         p2 <- function(x){
                            if ("lower.tail" %in% names(formals(p2fun)))
                                 p2fun(x, lower.tail)
                            else {p0 <- p2fun(x)
                                  if(!lower.tail) p0 <- 1-p0
                                  return(p0)} }
                         p0 <- (1-size)*p1(q) + size*p2(q)
                         if (log.p) p0 <- log(p0)
                         return(p0)
                                 },
                         list(size = size, p1fun = p(e1), p2fun = p(e2)))

        m1 <- min(q(e1)(TruncQuantile), q(e2)(TruncQuantile))
        m21 <- ifelse("lower.tail" %in% names(formals(e1@q)),
                      q(e1)(TruncQuantile, lower.tail = FALSE),
                      q(e1)(1-TruncQuantile))
        m22 <- ifelse("lower.tail" %in% names(formals(e2@q)),
                      q(e2)(TruncQuantile, lower.tail = FALSE),
                      q(e2)(1-TruncQuantile))
        m2 <- max(m21,m22); rm(m21,m22)

        qfun <- function(p, lower.tail = TRUE, log.p = FALSE){}
        body(qfun) <- substitute({ if (log.p) p <- exp(p)
                                   pfunx <- seq(from = m1, to = m2, length = 1e5)
                                   p0 <- pfun 
                                   pfuny <- p0(pfunx, lower.tail = lower.tail)
                                   qfun1 <- approxfun(x = pfuny, y = pfunx, 
                                                      rule = 2)
                                   y <- ifelse(p > 1, 
                                               NaN, 
                                               ifelse(p < 0, NaN, qfun1(p)))
                                   return(y)},
                            list(m1 = m1, m2 = m2, pfun = pfun)) 
    
        Symmetry <- NoSymmetry()
        if(is(e1@Symmetry,"SphericalSymmetry")&& 
           is(e2@Symmetry,"SphericalSymmetry"))
           if(.isEqual(SymmCenter(e1@Symmetry),SymmCenter(e2@Symmetry)))
              Symmetry <- SphericalSymmetry(SymmCenter(e1@Symmetry))   

        return(new("AbscontDistribution", r = rfun, d = dfun, p = pfun, 
                    q = qfun,  
                    .withSim = e1@.withSim|e2@.withSim, 
                    .withArith = e1@.withArith|e2@.withArith,
                    .logExact = FALSE, 
                    .lowerExact = e1@.lowerExact & e2@.lowerExact,
                    Symmetry = Symmetry))
    })

setMethod("ConvexContamination", signature(e1 = "DiscreteDistribution",
                                           e2 = "DiscreteDistribution",
                                           size = "numeric"),
    function(e1, e2, size){
        if(length(size) != 1)
            stop("length of 'size' has to be 1")
        if((size < 0)|(size > 1))
            stop("'size' has to be in [0,1]")
        supp <- union(support(e1), support(e2))
        len <- length(supp)
        if(length(usupp <- unique(supp)) < len){
            supp <- sort(usupp)
            len <- length(supp)
        }else{
            o <- order(supp)
            supp <- supp[o]
        }

        rfun <- function(n){}
        body(rfun) <- substitute({ r1 <- r1fun; r2 <- r2fun
                                   ind <- rbinom(n, prob=size, size=1)
                                   (1-ind)*r1(n) + ind*r2(n)},
                            list(size = size, r1fun = r(e1), r2fun = r(e2)))

        dfun <- function(x, log = FALSE){}
        body(dfun) <- substitute({ d1 <- d1fun; d2 <- d2fun
                                   d0 <- (1-size)*d1(x) + size*d2(x)
                                  if(log) d0 <- log(d0)
                                  return(d0) },
                            list(size = size, d1fun = d(e1), d2fun = d(e2)))

        pfun <- function(q, lower.tail = TRUE, log.p = FALSE){}
        body(pfun) <- substitute({ 
                         p1 <- function(x){
                            if ("lower.tail" %in% names(formals(p1fun)))
                                 p1fun(x, lower.tail)
                            else {p0 <- p1fun(x)
                                  if(!lower.tail) p0 <- 1-p0
                                  return(p0)} }     
                         p2 <- function(x){
                            if ("lower.tail" %in% names(formals(p2fun)))
                                 p2fun(x, lower.tail)
                            else {p0 <- p2fun(x)
                                  if(!lower.tail) p0 <- 1-p0
                                  return(p0)} }
                         p0 <- (1-size)*p1(q) + size*p2(q)
                         if (log.p) p0 <- log(p0)
                         return(p0)
                                 },
                         list(size = size, p1fun = p(e1), p2fun = p(e2)))

        cumprob.l <- pfun(supp)
        cumprob.u <- pfun(supp, lower.tail = FALSE)

        qfun <- function(p, lower.tail = TRUE, log.p = FALSE)
                         {if(log.p) p <- exp(p)
                          i01 <- (0<=p)&(p<=1)
                          p01 <- p[i01]
                          q0 <- p
                          q0[!i01] <- NaN
                          if (lower.tail)
                              q0[i01] <- sapply(p01, function(x){
                                       supp[ sum( cumprob.l < x ) + 1] })
                          else    
                              q0[i01] <- sapply(p01, function(x){ 
                                      (rev(supp))[sum( cumprob.u < x ) + 1] })
                          return(q0)
                          }        
    
        Symmetry <- NoSymmetry()
        if(is(e1@Symmetry,"SphericalSymmetry")&& 
           is(e2@Symmetry,"SphericalSymmetry"))
           if(.isEqual(SymmCenter(e1@Symmetry),SymmCenter(e2@Symmetry)))
              Symmetry <- SphericalSymmetry(SymmCenter(e1@Symmetry))   

        return(new("DiscreteDistribution", r = rfun, d = dfun, p = pfun, q = qfun, 
                    support = supp,  
                    .withSim = e1@.withSim|e2@.withSim, 
                    .withArith = e1@.withArith|e2@.withArith,
                    .logExact = FALSE, 
                    .lowerExact = e1@.lowerExact & e2@.lowerExact,
                    Symmetry = Symmetry))
    })

setMethod("ConvexContamination", signature(e1 = "UnivariateDistribution",
                                           e2 = "UnivariateDistribution",
                                           size = "numeric"),
    function(e1, e2, size){
        if(length(size) != 1)
            stop("length of 'size' has to be 1")
        if((size < 0)|(size > 1))
            stop("'size' has to be in [0,1]")
        rfun <- function(n){}
        body(rfun) <- substitute({ r1 <- r1fun; r2 <- r2fun
                                   ind <- rbinom(n, prob=size, size=1)
                                   (1-ind)*r1(n) + ind*r2(n)},
                            list(size = size, r1fun = r(e1), r2fun = r(e2)))

        pfun <- function(q, lower.tail = TRUE, log.p = FALSE){}
        body(pfun) <- substitute({ 
                         p1 <- function(x){
                            if ("lower.tail" %in% names(formals(p1fun)))
                                 p1fun(x, lower.tail)
                            else {p0 <- p1fun(x)
                                  if(!lower.tail) p0 <- 1-p0
                                  return(p0)} }     
                         p2 <- function(x){
                            if ("lower.tail" %in% names(formals(p2fun)))
                                 p2fun(x, lower.tail)
                            else {p0 <- p2fun(x)
                                  if(!lower.tail) p0 <- 1-p0
                                  return(p0)} }
                         p0 <- (1-size)*p1(q) + size*p2(q)
                         if (log.p) p0 <- log(p0)
                         return(p0)
                                 },
                         list(size = size, p1fun = p(e1), p2fun = p(e2)))

        TruncQuantile <- getdistrOption("TruncQuantile")
        m1 <- min(q(e1)(TruncQuantile), q(e2)(TruncQuantile))
        m21 <- ifelse("lower.tail" %in% names(formals(e1@q)),
                      q(e1)(TruncQuantile, lower.tail = FALSE),
                      q(e1)(1-TruncQuantile))
        m22 <- ifelse("lower.tail" %in% names(formals(e2@q)),
                      q(e2)(TruncQuantile, lower.tail = FALSE),
                      q(e2)(1-TruncQuantile))
        m2 <- max(m21,m22); rm(m21,m22)


        qfun <- function(p, lower.tail = TRUE, log.p = FALSE){}
        body(qfun) <- substitute({ if (log.p) p <- exp(p)
                                   pfunx <- seq(from = m1, to = m2, length = 1e5)
                                   p0 <- pfun 
                                   pfuny <- p0(pfunx, lower.tail = lower.tail)
                                   qfun1 <- approxfun(x = pfuny, y = pfunx, 
                                                      rule = 2)
                                   y <- ifelse(p > 1, 
                                               NaN, 
                                               ifelse(p < 0, NaN, qfun1(p)))
                                   return(y)},
                            list(m1 = m1, m2 = m2, pfun = pfun)) 
    
        Symmetry <- NoSymmetry()
        if(is(e1@Symmetry,"SphericalSymmetry")&& 
           is(e2@Symmetry,"SphericalSymmetry"))
           if(.isEqual(SymmCenter(e1@Symmetry),SymmCenter(e2@Symmetry)))
              Symmetry <- SphericalSymmetry(SymmCenter(e1@Symmetry))   

        return(new("UnivariateDistribution", img = img(e1), r = rfun, d = NULL, 
                    p = pfun, q = qfun, 
                    .withSim = e1@.withSim|e2@.withSim, 
                    .withArith = e1@.withArith|e2@.withArith,
                    .logExact = FALSE, 
                    .lowerExact = e1@.lowerExact & e2@.lowerExact,
                    Symmetry = Symmetry ))
    })


setMethod("ConvexContamination", signature(e1 = "AcDcLcDistribution",
                                           e2 = "AcDcLcDistribution",
                                           size = "numeric"),
    function(e1, e2, size){
        e1 <- as(e1, "UnivarLebDecDistribution")
        e2 <- as(e2, "UnivarLebDecDistribution")

        if(length(size) != 1)
            stop("length of 'size' has to be 1")
        if((size < 0)|(size > 1))
            stop("'size' has to be in [0,1]")
        
        return(flat.mix(UnivarMixingDistribution(e1, e2, 
                                    mixCoeff = c(1-size,size))))
    })
setMethod("ConvexContamination", signature(e1 = "LatticeDistribution", 
                                         e2 = "LatticeDistribution",
                                           size = "numeric"),
    getMethod("ConvexContamination", signature(e1 = "DiscreteDistribution", 
                                         e2 = "DiscreteDistribution",
                                           size = "numeric")))

setMethod("ConvexContamination", signature(e1 = "LatticeDistribution", 
                                         e2 = "DiscreteDistribution",
                                           size = "numeric"),
    getMethod("ConvexContamination", signature(e1 = "DiscreteDistribution", 
                                         e2 = "DiscreteDistribution",
                                           size = "numeric")))

setMethod("ConvexContamination", signature(e1 = "DiscreteDistribution", 
                                         e2 = "LatticeDistribution",
                                           size = "numeric"),
getMethod("ConvexContamination", signature(e1 = "DiscreteDistribution", 
                                         e2 = "DiscreteDistribution",
                                           size = "numeric")))
