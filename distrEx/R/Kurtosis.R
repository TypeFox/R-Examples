###################################################################################
#kurtosis  --- code due to G. Jay Kerns, gkerns@ysu.edu
###################################################################################


###################################################################################
#kurtosis
###################################################################################
setMethod("kurtosis", signature(x = "UnivariateDistribution"),
    function(x, fun = function(t) {t}, cond, withCond = FALSE, useApply = TRUE, ...){
        dots <- match.call(call = sys.call(sys.parent(1)), 
                        expand.dots = FALSE)$"..."
        low <- -Inf; upp <- Inf
        if(hasArg(low)) low <- dots$low
        if(hasArg(upp)) upp <- dots$upp
        LowIsUpp <- if(low == -Inf) 
                    low == -upp else .isEqual(low,upp)

        if(LowIsUpp && missing(cond)&&missing(fun)){
           if(is(Symmetry(x),"SphericalSymmetry")){
                m2 <- 2*E(x, fun = function(t)t^2, low=0, useApply = useApply, ...)
                m4 <- 2*E(x, fun = function(t)t^4, low=0, useApply = useApply, ...)
                return(m4/m2^2 -3)
           }
        }
        f2 <- function(t) {fun(t)^2}
        f3 <- function(t) {fun(t)^3}
        f4 <- function(t) {fun(t)^4}
        if(missing(cond))
            {
            m <- E(x, fun = fun, useApply = useApply, ...)
            m2 <- E(x, fun = f2, useApply = useApply, ...)
            m3 <- E(x, fun = f3, useApply = useApply, ...)
            m4 <- E(x, fun = f4, useApply = useApply, ...)

            return( (m4-4*m3*m+6*m2*m^2-3*m^4)/(var(x, fun = fun, 
                                                    useApply = TRUE, ...))^2 -3)
            }
        else{
            m <- E(x, cond = cond, fun = fun, withCond  = withCond, useApply = useApply, ...)
            m2 <- E(x, cond = cond, fun = f2, withCond  = withCond, useApply = useApply, ...)
            m3 <- E(x, cond = cond, fun = f3, withCond  = withCond, useApply = useApply, ...)
            m4 <- E(x, cond = cond, fun = f4, withCond  = withCond, useApply = useApply, ...)

            return( (m4-4*m3*m+6*m2*m^2-3*m^4)/(var(x, fun = fun, cond = cond, 
                                  withCond = FALSE, useApply = TRUE,...))^2  -3)

            }

    })

setMethod("kurtosis", signature(x = "AffLinDistribution"),
    function(x, fun = function(t) {t}, cond, withCond = FALSE, useApply = TRUE, ...){
        if (missing(fun) && missing(cond)){

            return( kurtosis(x@X0, withCond = withCond, useApply = useApply, 
                             ...))

            }
        else return(kurtosis( x = as(x, sub("AffLin","",class(x))), 
                    fun = fun, cond = cond, withCond = withCond, 
                    useApply = useApply, ... ))
    })

setMethod("kurtosis", signature(x = "AffLinAbscontDistribution"),
           getMethod("kurtosis", signature(x = "AffLinDistribution")))    
setMethod("kurtosis", signature(x = "AffLinDiscreteDistribution"),
           getMethod("kurtosis", signature(x = "AffLinDistribution")))    
setMethod("kurtosis", signature(x = "AffLinLatticeDistribution"),
           getMethod("kurtosis", signature(x = "AffLinDistribution")))    
  
    
###
# some exact kurtoses:
###
setMethod("kurtosis", signature(x = "Norm"),
    function(x,...){ 
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp))
       return(kurtosis(as(x,"AbscontDistribution"),...))
    else
        return(0)
    })
#
setMethod("kurtosis", signature(x = "Binom"),
    function(x,  ...){
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp))
       return(kurtosis(as(x,"DiscreteDistribution"),...))
    else
        p <- prob(x)
        return((1-6*p*(1-p))/(size(x)*p*(1-p)))
    })
### source: http://mathworld.wolfram.com/BinomialDistribution.html

#
setMethod("kurtosis", signature(x = "Cauchy"),
    function(x,...){    
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp))
      return(kurtosis(as(x,"AbscontDistribution"),...))
    else
        return(NA)
    })
### source http://mathworld.wolfram.com/CauchyDistribution.html

#
setMethod("kurtosis", signature(x = "Chisq"),
    function(x,...){    
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp))
       return(kurtosis(as(x,"AbscontDistribution"),...))
    else
        return(12*(df(x)+4*ncp(x))/(df(x)+2*ncp(x))^2)
    })
### source http://mathworld.wolfram.com/Chi-SquaredDistribution.html

#
setMethod("kurtosis", signature(x = "Dirac"),
    function(x, ...){return(NA)})

#
setMethod("kurtosis", signature(x = "DExp"),
    function(x, ...){    
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp)) 
         return(kurtosis(as(x,"AbscontDistribution"),...))
    else
        return(3)
    })
### source http://mathworld.wolfram.com/LaplaceDistribution.html

#
setMethod("kurtosis", signature(x = "Exp"),
    function(x, ...){    
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp)) 
         return(kurtosis(as(x,"AbscontDistribution"),...))
    else
        return(6)
    })
 ### source http://mathworld.wolfram.com/ExponentialDistribution.html

#
setMethod("kurtosis", signature(x = "Fd"),
    function(x, ...){
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp)) {
         return(kurtosis(as(x,"AbscontDistribution"),...))
    }else {
        if (df2(x)>8){
          m <- df1(x)
          n <- df2(x)
          d <- ncp(x)
          m2 <- var(x)
          m1 <- E(x)
          m3 <- (n/m)^3/(n-2)/(n-4)/(n-6)*
                  (m^3+6*m^2+8*m+3*d*(m^2+6*m+8)+3*d^2*(m+4)+d^3)
          mm1 <- m-1
          mm2 <- mm1 * (m+1)
          mm3 <- mm2 * (m+3)
          mm4 <- mm3 * (m+5)
          mmd1 <- d+1
          mmd2 <- 3 + 6*d + d^2
          mmd3 <- 15 + 45*d + 15*d^2 + d^3
          mmd4 <- 105 + 420*d + 210*d^2 + 28*d^3 + d^4
          mm <- mm4 + 4*mm3*mmd1 + 6*mm2*mmd2 + 4*mm1*mmd3+ mmd4          
          m4 <- (n/m)^4/(n-2)/(n-4)/(n-6)/(n-8)*mm
          return((m4-4*m3*m1+6*m2*m1^2+3*m1^4)/m2^2-3)
#          L <- d/m
#          m2 <- 2*n^2*(m+n-2)/m/(n-2)^2/(n-4)*(1+2*L+m*L^2/(m+n-2))
#          a <-  12*n^4*(m+n-2)/m^3/(n-2)^4/(n-4)/(n-6)/(n-8)
#          b <-  (1+4*L)*(2*(3*m+n-2)*(2*m+n-2)+(m+n-2)*(n-2)*(m+2))
#          c <-  2*m*(3*m+2*n-4)*(n+10)*L^2
#          d <-  4*m^2*(n+10)*L^3
#          e <-  m^3*(n+10)*L^4/(m+n-2)
#          m4 <- a*(b+c+d+e)
#          return(m4/m2^2-3)
        } else {
          return(NA)
        }
    }
    })
### source (without ncp) http://mathworld.wolfram.com/F-Distribution.html
#
setMethod("kurtosis", signature(x = "Gammad"),
    function(x, ...){    
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp)) 
         return(kurtosis(as(x,"AbscontDistribution"),...))
    else
        return(6/shape(x))
    })

### source http://mathworld.wolfram.com/GammaDistribution.html
#
setMethod("kurtosis", signature(x = "Geom"),
    function(x, ...){    
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp)) 
         return(kurtosis(as(x,"DiscreteDistribution"),...))
    else
        return(6+ prob(x)^2/(1-prob(x)))
    })
### source http://mathworld.wolfram.com/GeometricDistribution.html
#
setMethod("kurtosis", signature(x = "Hyper"),
    function(x, ...){    
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp)) 
         return(kurtosis(as(x,"DiscreteDistribution"),...))
    else
       {k <- k(x);
        m <- m(x); 
        n <- n(x);
        a <- (m+n)^2*(m+n-1)/(k*m*n*(m+n-k)*(m+n-2)*(m+n-3));
        return(
                a*((m+n)*(m+n+1-6*k)+3*m*n*(k-2)+6*k^2+3*m*n*k*(6-k)/(m+n)
                -18*m*n*k^2/(m+n)^2)-3
              )
        }
    })
### source http://mathworld.wolfram.com/HypergeometricDistribution.html
#
setMethod("kurtosis", signature(x = "Logis"),
    function(x, ...){
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp)) 
        return(kurtosis(as(x,"AbscontDistribution"),...))
    else
        return(6/5)
    })
### source http://mathworld.wolfram.com/LogisticDistribution.html
#
setMethod("kurtosis", signature(x = "Lnorm"),
    function(x, ...){
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp)) {
        return(kurtosis(as(x,"AbscontDistribution"),...))
    } else {
        w <- exp(sdlog(x)^2)
        return( w^4+2*w^3+3*w^2-6)
    }
    })
### source http://mathworld.wolfram.com/LogNormalDistribution.html
#
setMethod("kurtosis", signature(x = "Nbinom"),
    function(x, ...){    
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp)) 
         return(kurtosis(as(x,"DiscreteDistribution"),...))
    else
        return(6/size(x)+prob(x)^2/(size(x)*(1-prob(x))))
    })
### source http://mathworld.wolfram.com/NegativeBinomialDistribution.html
#
setMethod("kurtosis", signature(x = "Pois"),
    function(x, ...){
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp)) 
        return(kurtosis(as(x,"DiscreteDistribution"),...))
    else
        return(1/lambda(x))
    })
### source http://mathworld.wolfram.com/PoissonDistribution.html
#
setMethod("kurtosis", signature(x = "Td"),
    function(x, ...){
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp)){ 
        return(kurtosis(as(x,"AbscontDistribution"),...))
    } else {
        if (df(x)>4){
          n <- df(x)
          d <- ncp(x)
          m2 <- var(x)
          m1 <- E(x)
          m3 <- (n/2)^1.5*(3*d+d^3)*exp(lgamma((n-3)/2)-lgamma(n/2))
          m4 <- n^2*(3+6*d^2+d^4)/(n-2)/(n-4)
          return((m4-4*m3*m1+6*m2*m1^2+3*m1^4)/m2^2-3)
        } else {
          return(NA)
        }
    }
    })
### source http://mathworld.wolfram.com/NoncentralStudentst-Distribution.html

#
setMethod("kurtosis", signature(x = "Unif"),
    function(x, ...){
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp)) 
        return(kurtosis(as(x,"AbscontDistribution"),...))
    else
        return(-6/5)
    })
### source http://mathworld.wolfram.com/UniformDistribution.html
#
setMethod("kurtosis", signature(x = "Weibull"),
    function(x, ...){
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp)) 
        return(kurtosis(as(x,"AbscontDistribution"),...))
    else
        g1 <- gamma(1+1/shape(x))
        g2 <- gamma(1+2/shape(x))
        g3 <- gamma(1+3/shape(x))
        g4 <- gamma(1+4/shape(x))
        v <- (g2-g1^2)^2
        return( (g4-4*g3*g1+6*g2*g1^2-3*g1^4)/v - 3 )
    })
### source http://mathworld.wolfram.com/WeibullDistribution.html
#    
setMethod("kurtosis", signature(x = "Beta"),
    function(x, ...){
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if((hasArg(fun))||(hasArg(cond))||(!isTRUE(all.equal(ncp(x),0)))) 
        return(kurtosis(as(x,"AbscontDistribution"),...))
    else
        {a<-shape1(x); b<- shape2(x)
        return(6*(a^3-a^2*(2*b-1)+b^2*(b+1)-2*a*b*(b+2))/(a*b*(a+b+2)*(a+b+3)) )}
    })
## source: http://mathworld.wolfram.com/BetaDistribution.html

###################################################################################
#kurtosis --- code P.R.:
###################################################################################

setMethod("kurtosis", signature(x = "Arcsine"),
    function(x, ...){
    dots <- match.call(call = sys.call(sys.parent(1)), 
                       expand.dots = FALSE)$"..."
    fun <- NULL; cond <- NULL; low <- NULL; upp <- NULL
    if(hasArg(low)) low <- dots$low
    if(hasArg(upp)) upp <- dots$upp
    if(hasArg(fun)||hasArg(cond)||!is.null(low)||!is.null(upp)) 
        return(kurtosis(as(x,"AbscontDistribution"),...))
    else    return(-3/2)
    })

