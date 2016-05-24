###############################################################################
## Method: ContaminationSize
## size of contamination for two distributions
###############################################################################
setMethod("ContaminationSize", signature(e1 = "AbscontDistribution", 
                                         e2 = "AbscontDistribution"),
    function(e1, e2){
        ep <- getdistrOption("TruncQuantile")
        lower <- min(q(e1)(ep), q(e2)(ep))
        upper <- max(q(e1)(1-ep), q(e2)(1-ep))
        x <- seq(from = lower, to = upper, length = 1e5)
        
        d10  <- d(e1)(x); d1 <- d10[ d10>0 ]
        d20  <- d(e2)(x); d2 <- d20[ d10>0 ]
        res <- 1 - min(d2/d1)
        return(list(e1 = e1, e2 = e2, size.of.contamination = res))
    })

setMethod("ContaminationSize", signature(e1 = "DiscreteDistribution", 
                                         e2 = "DiscreteDistribution"),
    function(e1, e2){
        o.warn <- getOption("warn"); options(warn = -1)
        on.exit(options(warn=o.warn))
        x <- union(support(e1), support(e2))
        d10  <- d(e1)(x); d1 <- d10[ d10>0 ]
        d20  <- d(e2)(x); d2 <- d20[ d10>0 ]
        
        res <- min(1- min(d2/d1),1)



        return(list(e1 = e1, e2 = e2, size.of.contamination = res))
    })
#### new from version 2.0 on: Distance for Mixing Distributions
setMethod("ContaminationSize",  signature(e1 = "AcDcLcDistribution",
                                     e2 = "AcDcLcDistribution"),
           function(e1, e2){
           if( is(e1,"AbscontDistribution"))
               e1 <- as(as(e1,"AbscontDistribution"), "UnivarLebDecDistribution")
           if( is(e2,"AbscontDistribution"))
               e2 <- as(as(e2,"AbscontDistribution"), "UnivarLebDecDistribution")
           if(is(e1,"DiscreteDistribution"))
               e1 <- as(as(e1,"DiscreteDistribution"), "UnivarLebDecDistribution")
           if(is(e2,"DiscreteDistribution"))
               e2 <- as(as(e2,"DiscreteDistribution"), "UnivarLebDecDistribution")
              ac1 <- acPart(e1); ac2 <- acPart(e2)
              ac1d <- ac1@d; ac2d <- ac2@d
              ac1@d <- function(x) ac1d(x)*acWeight(e1)
              ac2@d <- function(x) ac2d(x)*acWeight(e2)
              dc1 <- discretePart(e1); dc2 <- discretePart(e2)
              dc1d <- dc1@d; dc2d <- dc2@d
              dc1@d <- function(x) dc1d(x)*discreteWeight(e1)
              dc2@d <- function(x) dc2d(x)*discreteWeight(e2)
              res <- max(ContaminationSize(ac1,ac2)$size.of.contamination,
                         ContaminationSize(dc1,dc2)$size.of.contamination)
              return(list(e1 = e1, e2 = e2, size.of.contamination = res))
              })

setMethod("ContaminationSize", signature(e1 = "LatticeDistribution", 
                                         e2 = "LatticeDistribution"),
    getMethod("ContaminationSize", signature(e1 = "DiscreteDistribution", 
                                         e2 = "DiscreteDistribution")))

setMethod("ContaminationSize", signature(e1 = "LatticeDistribution", 
                                         e2 = "DiscreteDistribution"),
    getMethod("ContaminationSize", signature(e1 = "DiscreteDistribution", 
                                         e2 = "DiscreteDistribution")))

setMethod("ContaminationSize", signature(e1 = "DiscreteDistribution", 
                                         e2 = "LatticeDistribution"),
    getMethod("ContaminationSize", signature(e1 = "DiscreteDistribution", 
                                         e2 = "DiscreteDistribution")))
