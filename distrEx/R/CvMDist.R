###############################################################################
## Method: CvMDist
## Cramer - von Mises distance of two distributions
###############################################################################
setMethod("CvMDist", signature(e1 = "UnivariateDistribution",
                                    e2 = "UnivariateDistribution"),
    function(e1, e2, mu = e1, useApply = FALSE, ... ){
        o.warn <- getOption("warn"); options(warn = -1)
        on.exit(options(warn=o.warn))
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
        res <- E(mu, fun = function(t) {(p(e1)(t)-p(e2)(t))^2}, useApply = useApply, ...)^.5
        names(res) <- "CvM distance"
        return(res)
    })

## CvM distance
setMethod("CvMDist", signature(e1 = "numeric",
                                    e2 = "UnivariateDistribution"),
    function(e1, e2, mu = e1, ...)
        { o.warn <- getOption("warn"); options(warn = -1)
          on.exit(options(warn=o.warn))
          if(identical(mu,e2))
             return(.newCvMDist(e1,e2))
          e10 <- DiscreteDistribution(e1)       
          if(identical(mu,e1)) mu <- e10
          CvMDist(e1 = e10, e2 = e2, mu = mu, ...)
         }
    )

### new Method if mu=e2: explicit integration...
.newCvMDist <- function(e1,e2){
 ### use that 
 ##  int (P_n(t)-P(t))^2 P(dt) = 
 ###      1/n^2 sum_ij (1-P(max(xi, xj)) - 
 ###      1/n sum_i (1-P(xi)^2) + 
 ###      (P(infty)^3-P(-infty)^3)/3 = 
 ###      1/n^2 sum_i (2i-1) (1-P(x(i))) -  # x(i) ordnungsstatistik
 ###      1/n sum_i (1-P(xi)^2) + 1/3
 ### on my Laptop ~ 30 times faster!!
 x1 <- sort(e1)
 p <- p(e2)
 p0 <- p(x1)
 p1 <- if("lower" %in% names(formals(p))) 
          p(e2)(x1,lower=FALSE) else 1-p0
 p2 <- 1-p0^2
 n <- length(x1)
 i1 <- 2*(1:n)-1
 d <- (mean(i1*p1)/n-mean(p2)+1/3)^.5
 names(d) <- "CvM distance" 
 return(d)
}
