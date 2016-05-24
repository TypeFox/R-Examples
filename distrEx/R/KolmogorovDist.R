###############################################################################
## Method: KolmogorovDist
## Kolmogorov distance of two distributions
###############################################################################
setMethod("KolmogorovDist", signature(e1 = "AbscontDistribution",
                                      e2 = "AbscontDistribution"),
    function(e1, e2){
        TruncQuantile <- getdistrOption("TruncQuantile")  
        lower1 <- ifelse(!is.finite(q(e1)(0)), q(e1)(TruncQuantile), q(e1)(0))
        upper1 <- ifelse(!is.finite(q(e1)(1)), 
                         ifelse("lower.tail" %in% names(formals(e1@q)),
                                q(e1)(TruncQuantile, lower.tail = FALSE),
                                q(e1)(1-TruncQuantile)), 
                         q(e1)(1))
        lower2 <- ifelse(!is.finite(q(e2)(0)), q(e2)(TruncQuantile), q(e2)(0))
        upper2 <- ifelse(!is.finite(q(e2)(1)), 
                         ifelse("lower.tail" %in% names(formals(e2@q)),
                                q(e2)(TruncQuantile, lower.tail = FALSE),
                                q(e2)(1-TruncQuantile)), 
                         q(e2)(1))
        lower <- min(lower1, lower2)
        upper <- max(upper1, upper2)

        o.warn <- getOption("warn"); options(warn = -1)
        on.exit(options(warn=o.warn))
        x1 <- union(r(e1)(1e5), r(e2)(1e5))
        x2 <- seq(from=lower, to=upper, length=1e5)
        x <- union(x1, x2) 

        res <- max(abs(p(e1)(x)-p(e2)(x)))
        names(res) <- "Kolmogorov distance"


        return(res)
    })

setMethod("KolmogorovDist", signature(e1 = "DiscreteDistribution",
                                      e2 = "DiscreteDistribution"),
    function(e1, e2){
        o.warn <- getOption("warn"); options(warn = -1)
        on.exit(options(warn=o.warn))
        supp <- union(support(e1), support(e2))

        res <- max(abs(p(e1)(supp)-p(e2)(supp)))
        names(res) <- "Kolmogorov distance"

        return(res)
    })

setMethod("KolmogorovDist", signature(e1 = "DiscreteDistribution",
                                      e2 = "AbscontDistribution"),
    function(e1, e2){
        o.warn <- getOption("warn"); options(warn = -1)
        on.exit(options(warn=o.warn))
        x <- support(e1)
        res <- max(p(e1)(x)-p(e2)(x),p(e2)(x)-p.l(e1)(x))
        names(res) <- "Kolmogorov distance"

        return(res)
    })

setMethod("KolmogorovDist", signature(e1 = "AbscontDistribution",
                                      e2 = "DiscreteDistribution"),
    function(e1, e2){
        KolmogorovDist(e2, e1)
    })
## Kolmogorov distance
setMethod("KolmogorovDist", signature(e1 = "numeric",
                                      e2 = "UnivariateDistribution"),
    function(e1, e2){
        o.warn <- getOption("warn")
        options(warn = -1)
        emp <- DiscreteDistribution(e1)
        return(KolmogorovDist(emp,e2))
    })

setMethod("KolmogorovDist", signature(e1 = "UnivariateDistribution",
                                      e2 = "numeric"),
    function(e1, e2){
        return(KolmogorovDist(e2, e1))
    })

#### new from version 2.0 on: Distance for Mixing Distributions
setMethod("KolmogorovDist",  signature(e1 = "AcDcLcDistribution",
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
        if(is.null(e1@p)){
        e1.erg <- RtoDPQ(e1@r)
        e1 <- new("UnivariateDistribution", r=e1@r,
                   p = e1.erg$pfun, d = e1.erg$dfun, q = e1.erg$qfun,
                   .withSim = TRUE, .withArith = FALSE)}
        if(is.null(e2@p)){
        e2.erg <- RtoDPQ(e2@r)
        e2 <- new("UnivariateDistribution", r=e2@r,
                   p = e2.erg$pfun, d = e2.erg$dfun, q = e2.erg$qfun,
                   .withSim = TRUE, .withArith = FALSE)}
        TruncQuantile <- getdistrOption("TruncQuantile")
        lower1 <- ifelse(!is.finite(q(e1)(0)), q(e1)(TruncQuantile), q(e1)(0))
        upper1 <- ifelse(!is.finite(q(e1)(1)),
                         ifelse("lower.tail" %in% names(formals(e1@q)),
                                q(e1)(TruncQuantile, lower.tail = FALSE),
                                q(e1)(1-TruncQuantile)),
                         q(e1)(1))
        lower2 <- ifelse(!is.finite(q(e2)(0)), q(e2)(TruncQuantile), q(e2)(0))
        upper2 <- ifelse(!is.finite(q(e2)(1)),
                         ifelse("lower.tail" %in% names(formals(e2@q)),
                                q(e2)(TruncQuantile, lower.tail = FALSE),
                                q(e2)(1-TruncQuantile)),
                         q(e2)(1))
        lower <- min(lower1, lower2)
        upper <- max(upper1, upper2)

        o.warn <- getOption("warn"); options(warn = -1)
        on.exit(options(warn=o.warn))
        x1 <- union(r(e1)(1e5), r(e2)(1e5))
        x2 <- seq(from=lower, to=upper, length=1e5)
        x <- union(x1, x2)

        if( "support" %in% names(getSlots(class(e1))))
           x <- union(x,e1@support)
        if( "support" %in% names(getSlots(class(e2))))
           x <- union(x,e2@support)

        res <- max(abs(p(e1)(x)-p(e2)(x)))
        names(res) <- "Kolmogorov distance"


        return(res)
    })
setMethod("KolmogorovDist", signature(e1 = "LatticeDistribution", 
                                         e2 = "LatticeDistribution"),
    getMethod("KolmogorovDist", signature(e1 = "DiscreteDistribution", 
                                         e2 = "DiscreteDistribution")))

setMethod("KolmogorovDist", signature(e1 = "LatticeDistribution", 
                                         e2 = "DiscreteDistribution"),
    getMethod("KolmogorovDist", signature(e1 = "DiscreteDistribution", 
                                         e2 = "DiscreteDistribution")))

setMethod("KolmogorovDist", signature(e1 = "DiscreteDistribution", 
                                         e2 = "LatticeDistribution"),
    getMethod("KolmogorovDist", signature(e1 = "DiscreteDistribution", 
                                         e2 = "DiscreteDistribution")))
