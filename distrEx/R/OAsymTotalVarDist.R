###############################################################################
## Method: OAsymOAsymTotalVarDist
## asymmetric total variation distance of two distributions
##
##  newly introduced by P.R. 13 03 09 --- for use in collaboration with 
##  Matthias Brandl, Claudio Agostinelli
##
###############################################################################
setMethod("OAsymTotalVarDist", signature(e1 = "AbscontDistribution", 
                                        e2 = "AbscontDistribution"),
    function(e1, e2,  
             rel.tol = .Machine$double.eps^0.3, Ngrid = 10000,
             TruncQuantile = getdistrOption("TruncQuantile"),
             IQR.fac = 15){ 

        ## if we have to recall this method with a smaller TruncQuantile arg:
        mc <-  as.list(match.call(call = sys.call(sys.parent(1)))[-1])
        mc$TruncQuantile <- TruncQuantile * 1.8

        #block warnings:
        o.warn <- getOption("warn"); options(warn = -1)
        on.exit(options(warn=o.warn))
        
        ## find sensible lower and upper bounds for integration of q-cp
        # (a) quantile based
        low <- min(getLow(e1, eps = TruncQuantile), getLow(e2, eps = TruncQuantile))
        up  <- max(getUp(e1, eps = TruncQuantile), getUp(e2, eps = TruncQuantile))
        # (b) scale based
        s0 <- min(IQR(e1),IQR(e2))*IQR.fac
        low0 <- min(median(e1),median(e2))-s0
        up0 <- max(median(e1),median(e2))+s0
        # (a) & (b)
        low <- max(low,low0); up <- min(up,up0)
        #
        # store densities to avoid dispatch
        d1 <- d(e1); d2 <- d(e2)
        #
        ### integration as a function of c:
        Eip <- function(f, c00)
                distrExIntegrate(f, lower = low, upper = up, 
                                    rel.tol = rel.tol, c00 = c00)
       integ <- function(x,c00)  abs(d2(x)-c00*d1(x))

       fct <- function(c0) Eip(f=integ, c00=c0)
       
       ## find sensible search range for c-values
       ## goal: range of density quotient d2(x)/d1(x)
       ## x-range:
       x.range <- seq(low, up, length=Ngrid/3)
       x.range <- c(x.range, q(e1)(seq(TruncQuantile,1-TruncQuantile,length=Ngrid/3)))
       x.range <- c(x.range, q(e2)(seq(TruncQuantile,1-TruncQuantile,length=Ngrid/3)))
       ## to avoid division by 0:
       d1x.range <- d10x.range <- d1(x.range)
       d1x.range <- d1x.range+(d1x.range<1e-20)
       ## bound dx.range from 0 and from maximal value:
       d2x.range <- d2(x.range)
       dx.range <- (d2x.range/d1x.range*(d10x.range>=1e-10)+
                   d2x.range*1e10*(d10x.range<1e-10))*.9999999+1e-10
       ## gives range for c:
       low1 <- min(dx.range); up1 <- max(dx.range)
       
       c.opt <- try(optimize(fct, lower = low1, upper = up1, 
                            tol = rel.tol)$minimum, 
                    silent = TRUE)

       ## if does not give reasonable solution recall function with 
       #  smaller TruncQuantile
       if(!is.numeric(c.opt)) 
           return(do.call(getMethod("OAsymTotalVarDist", 
                              signature(e1 = "AbscontDistribution", 
                                        e2 = "AbscontDistribution")), 
                              args = mc))             
       ## else:
       res <- Eip(f=integ, c00=c.opt)
       names(res) <- "minimal asym. total variation distance"
       return(res)
    })

setMethod("OAsymTotalVarDist", signature(e1 = "DiscreteDistribution",
                                    e2 = "DiscreteDistribution"),
    function(e1, e2,  ...){
        o.warn <- getOption("warn"); options(warn = -1)
        on.exit(options(warn=o.warn))
        supp <- union(support(e1), support(e2))
        # store densities to avoid dispatch
        d1 <- d(e1); d2 <- d(e2)

        d2.range <- d2(supp)
        d1.range <- d1(supp)
        
        integ <- function(c00) abs(d2.range-c00*d1.range)

        fct <- function(c0) sum(integ(c0))
       
       d10.range <- d1.range
       d1e.range <- d1.range+(d1.range<1e-10)
       ## bound dx.range from 0 and from maximal value:
       d.range <- (d2.range/d1e.range*(d10.range>=1e-10)+
                   d2.range*1e10*(d10.range<1e-10))*.99999999+1e-10
       ## gives range for c:
       low1 <- min(d.range); up1 <- max(d.range)
       
       
       c.opt <- optimize(fct, lower=low1, upper=up1)$minimum
       res <- sum(integ(c.opt))
       names(res) <- "minimal asym. total variation distance"
       return(res)
    })
setMethod("OAsymTotalVarDist", signature(e1 = "DiscreteDistribution",
                                    e2 = "AbscontDistribution"),
    function(e1, e2,  ...){
        res <- 1
        names(res) <- "minimal asym. total variation distance"

        return(res)
    })
setMethod("OAsymTotalVarDist", signature(e1 = "AbscontDistribution",
                                    e2 = "DiscreteDistribution"),
    function(e1, e2,  ...){ 
        res <- 1
        names(res) <- "minimal asym. total variation distance"

        return(res)
    })
setMethod("OAsymTotalVarDist", signature(e1 = "numeric",
                                    e2 = "DiscreteDistribution"),
    function(e1, e2,  ...){
        t1 <- table(e1)
        d1 <- t1/length(e1)
        s1 <- as.numeric(names(t1))
        e11 <- DiscreteDistribution(supp=s1, prob=d1)
        return(OAsymTotalVarDist(e11,e2,  ...))
    })

setMethod("OAsymTotalVarDist", signature(e1 = "DiscreteDistribution",
                                    e2 = "numeric"),
    function(e1, e2, ...){
        return(OAsymTotalVarDist(e2, e1,  ...))
    })

## to avoid trivial distances (distance = 1)
## abs.cont. distributions may be discretized
## resp. empirical distributions may be smoothed 
## (by convolution with a normal distribution)
setMethod("OAsymTotalVarDist", signature(e1 = "numeric",
                                    e2 = "AbscontDistribution"),
     function(e1, e2, asis.smooth.discretize = "discretize", n.discr =
             getdistrExOption("nDiscretize"), low.discr = getLow(e2),
             up.discr = getUp(e2), h.smooth = getdistrExOption("hSmooth"),
             rel.tol = .Machine$double.eps^0.3,  Ngrid = 10000,
             TruncQuantile = getdistrOption("TruncQuantile"),
             IQR.fac = 15){
        .asis.smooth.discretize.distance(e1, e2, asis.smooth.discretize, n.discr,
                 low.discr, up.discr, h.smooth, OAsymTotalVarDist, 
                 rel.tol = rel.tol, Ngrid = Ngrid,
                 TruncQuantile = TruncQuantile, IQR.fac = IQR.fac)
     })
setMethod("OAsymTotalVarDist", signature(e1 = "AbscontDistribution",
                                     e2 = "numeric"),
    function(e1, e2, asis.smooth.discretize = "discretize", n.discr =
             getdistrExOption("nDiscretize"), low.discr = getLow(e1),
             up.discr = getUp(e1), h.smooth = getdistrExOption("hSmooth"),
             rel.tol = .Machine$double.eps^0.3,  Ngrid = 10000,
             TruncQuantile = getdistrOption("TruncQuantile"),
             IQR.fac = 15){
        return(OAsymTotalVarDist(e2, e1,
                  asis.smooth.discretize = asis.smooth.discretize, 
                  low.discr = low.discr, up.discr = up.discr, h.smooth = h.smooth,
                  rel.tol = rel.tol, Ngrid = Ngrid,
                  TruncQuantile = TruncQuantile, IQR.fac = IQR.fac))
    })

setMethod("OAsymTotalVarDist",  signature(e1 = "AcDcLcDistribution",
                                     e2 = "AcDcLcDistribution"),
           function(e1, e2, 
             rel.tol = .Machine$double.eps^0.3, Ngrid = 10000,
             TruncQuantile = getdistrOption("TruncQuantile"),
             IQR.fac = 15){
        ## if we have to recall this method with a smaller TruncQuantile arg:
        mc <-  as.list(match.call(call = sys.call(sys.parent(1)))[-1])
        mc$TruncQuantile <- TruncQuantile * 1.8

        #block warnings:
        o.warn <- getOption("warn"); options(warn = -1)
        on.exit(options(warn=o.warn))

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
           ac1.d <- function(x) ac1d(x)*acWeight(e1)
           ac2.d <- function(x) ac2d(x)*acWeight(e2)

           dc1 <- discretePart(e1); dc2 <- discretePart(e2)
           dc1d <- dc1@d; dc2d <- dc2@d
           dc1.d <- function(x) dc1d(x)*discreteWeight(e1)
           dc2.d <- function(x) dc2d(x)*discreteWeight(e2)

        ### continuous part
        
        ## find sensible lower and upper bounds for integration of q-cp
        # (a) quantile based
        low <- min(getLow(ac1, eps = TruncQuantile), getLow(ac2, eps = TruncQuantile))
        up  <- max(getUp(ac1, eps = TruncQuantile), getUp(ac2, eps = TruncQuantile))
        # (b) scale based
        s0 <- min(IQR(ac1),IQR(ac2))*IQR.fac
        low0 <- min(median(ac1),median(ac2))-s0
        up0 <- max(median(ac1),median(ac2))+s0
        # (a) & (b)
        low <- max(low,low0); up <- min(up,up0)
        #
        ### integration as a function of c:
        Eip <- function(f, c00)
                distrExIntegrate(f, lower = low, upper = up, 
                                    rel.tol = rel.tol, c00 = c00)
       integ.c <- function(x,c00)  abs(ac2.d(x)-c00*ac1.d(x))

       ### discrete part

       supp <- union(support(dc1), support(dc2))

        d2.range <- dc2.d(supp)
        d1.range <- dc1.d(supp)

        integ.d <- function(c00) abs(d2.range-c00*d1.range)


       fct <- function(c0){
           e.c <- Eip(f=integ.c, c00=c0)
           e.d <- sum(integ.d(c0))
           return(e.d+e.c)
           } 
       
       ## find sensible search range for c-values
       ## goal: range of density quotient d2(x)/d1(x)
       
       ### continuous part
       ## x-range:
       x.range <- seq(low, up, length=Ngrid/3)
       x.range <- c(x.range, q(ac1)(seq(TruncQuantile,1-TruncQuantile,length=Ngrid/3)))
       x.range <- c(x.range, q(ac2)(seq(TruncQuantile,1-TruncQuantile,length=Ngrid/3)))
       ## to avoid division by 0:
       d1x.range <- d10x.range <- ac1.d(x.range)
       d1x.range <- d1x.range+(d1x.range<1e-20)
       ## bound dx.range from 0 and from maximal value:
       d2x.range <- ac2.d(x.range)
       dx.range <- (d2x.range/d1x.range*(d10x.range>=1e-20)+
                   d2x.range*1e20*(d10x.range<1e-20))*.99+1e-5

       ### discrete part
       d10.range <- d1.range
       d1e.range <- d1.range+(d1.range<1e-7)
       ## bound dx.range from 0 and from maximal value:
       d.range <- (d2.range/d1e.range*(d10.range>=1e-10)+
                   d2.range*1e10*(d10.range<1e-10))*.999999+1e-10

       ## gives range for c:
       low1 <- min(d.range,dx.range); up1 <- max(d.range,dx.range)

       c.opt <- try(optimize(fct, lower = low1, upper = up1, tol = rel.tol)$minimum, 
                    silent = TRUE)

       ## if does not give reasonable solution recall function with 
       #  smaller TruncQuantile
       if(!is.numeric(c.opt)){ 
           return(do.call(getMethod("OAsymTotalVarDist", 
                           signature(e1 = "AcDcLcDistribution",
                                     e2 = "AcDcLcDistribution")), 
                              args = mc))             
       }
       res <- Eip(f=integ.c, c00=c.opt)+sum(integ.d(c.opt))
       names(res) <- "minimal asym. total variation distance"
       return(res)
              })

setMethod("OAsymTotalVarDist", signature(e1 = "LatticeDistribution", 
                                         e2 = "LatticeDistribution"),
    getMethod("OAsymTotalVarDist", signature(e1 = "DiscreteDistribution", 
                                         e2 = "DiscreteDistribution")))

setMethod("OAsymTotalVarDist", signature(e1 = "LatticeDistribution", 
                                         e2 = "DiscreteDistribution"),
    getMethod("OAsymTotalVarDist", signature(e1 = "DiscreteDistribution", 
                                         e2 = "DiscreteDistribution")))

setMethod("OAsymTotalVarDist", signature(e1 = "DiscreteDistribution", 
                                         e2 = "LatticeDistribution"),
    getMethod("OAsymTotalVarDist", signature(e1 = "DiscreteDistribution", 
                                         e2 = "DiscreteDistribution")))


