## Generating function

EllipticalDistribution <- function(radDistr = sqrt(Chisq(df = length(loc))),
                            loc = c(0,0), scale = diag(length(loc)), p = NULL, q = NULL){

   ldscale <- as.numeric(determinant(as.matrix(scale),
                         logarithm = TRUE)$modulus)
   Iscale <- solve(scale)

   dim0 <- length(loc)

   param <- new("EllipticalParameter", loc=loc, scale=scale)

   if(!is(radDistr,"UnivariateDistribution"))
       stop("must be a univariate Distribution")
   if(p.l(radDistr)(0)>0)
      stop("distr must have pos. support")

   dr <- d(radDistr)
   dlog <- if(distr:::.inArgs("log", dr))
           quote(dr(r, log = TRUE)) else quote(log(dr(r)))

   if(is(radDistr,"AbscontDistribution")){
     dfun <- function(x, log = FALSE){}
     body(dfun) <- substitute({
          x0 <- x-loc0
          x1 <- Iscale0 %*% x0
          r <- colSums(x0*x1)^.5
          lg <- dlog0
          lg <- lg + (1-dim1)*log(r) + lgamma(dim1/2) -
                dim1/2*log(pi)-log(2)- ldscale0/2
          return(if(log) lg else exp(lg))},
          list(loc0 = loc, Iscale0 = Iscale, ldscale0 = ldscale,
               dlog0 = dlog, dim1=dim0))
    }else dfun <- NULL

    rfun <- function(n){}
    body(rfun) <- substitute({
         r0 <- r(radDistr)(n)
         u0 <- matrix(rnorm(n*dim1),ncol=dim1)
         u0n <- rowSums(u0^2)^.5
         un <- t(u0/u0n*r0)
         scale0 %*% un + loc0
      }, list(scale0=scale, loc0 = loc, dim1 = dim0))

    img0 <- new("EuclideanSpace", dimension = round(dim0,0))

    new("EllipticalDistribution", 
        r=rfun, d=dfun, p=p, q=q,
        radDistr = radDistr,
        img = img0, param = param,
        .withSim = radDistr@.withSim,
        .withArith = radDistr@.withArith,
        .logExact = radDistr@.logExact,
        .lowerExact = radDistr@.lowerExact,
        Symmetry = EllipticalSymmetry(loc))
   }

## Parameter for Elliptically symmetric Distribution
# accessors

setMethod("scale", "EllipticalParameter",
           function(x,  center, scale) x@scale)
setMethod("location", "EllipticalParameter",
           function(object) object@loc)

# replacements
setReplaceMethod("scale", "EllipticalParameter",
      function(object, value){ new("EllipticalParameter",
                                                loc = object@loc,
                                                scale = as.matrix(value))})
setReplaceMethod("location", "EllipticalParameter",
      function(object, value) new("EllipticalParameter", loc = value,
                                   scale = object@scale))

## Elliptically symmetric Distribution
# accessors
setMethod("scale", "EllipticalDistribution",
           function(x,  center, scale) (x@param)@scale)
setMethod("location", "EllipticalDistribution",
           function(object) (object@param)@loc)
# replacements
setReplaceMethod("scale", "EllipticalDistribution",
      function(object, value){   param <- new("EllipticalParameter",
                                              loc = object@param@loc,
                                              scale= as.matrix(value))
                          object@param <- param; object})
setReplaceMethod("location", "EllipticalDistribution",
      function(object, value){   param <- new("EllipticalParameter",
                                              loc = value,
                                              scale = object@param@scale)
                          object@param <- param; object})


setAs("UnivariateDistribution", "EllipticalDistribution",
      function(from){
        if(!is(Symmetry(from),"SphericalSymmetry"))
            return(from)
        else{ sc <- SymmCenter(Symmetry(from))
              radDistr <- abs(from-sc)
              ell <- EllipticalDistribution (radDistr = radDistr,
                            loc = sc, scale = 1, p = from@p, q = from@q)
              ell@r <- from@r
              ell@d <- from@d              
              ell@.withArith <- from@.withArith  
              ell@.lowerExact <- from@.lowerExact  
              ell@.logExact <- from@.logExact 
              return(ell)}
})



setAs("EllipticalDistribution", "UnivariateDistribution", 
      function(from){
        if(dimension(from)>1) return(from)
        radD <- radDistr(from)
        sca <- scale(from)
        loc <- location(from)
        if(!is(radD,"AcDcLcDistribution")){
            rfun <- function(n) sca * r(radD)(n) * 
                            sample(c(-1,1),n,replace=TRUE) + loc
            D <- new("UnivariateDistribution", r = rfun)
        }
        else{
             D <- radD * DiscreteDistribution(sca*c(-1,1)) + loc 
        }    
        D@Symmetry <- SphericalSymmetry(loc) 
        return(D)
      })

## functionals:
setMethod("E", signature(object = "EllipticalDistribution",
                        fun = "missing", cond = "missing"),
           function(object,...) location(object))
setMethod("var", signature(x = "EllipticalDistribution"),
           function(x,...) scale(x)%*%t(scale(x)) *
                    E(radDistr(x),fun=function(y)y^2,...)/dimension(x)
           )

setMethod("+", c("EllipticalDistribution","numeric"),
           function(e1,e2){ if(dimension(e1)!=length(e2))
                               stop("Dimension mismatch of operands in '+'")
                            location(e1) <- location(e1)+e2
                            return(e1)})   
setMethod("*", c("EllipticalDistribution","numeric"),
           function(e1,e2){ if((length(e2)!=1)&&length(e2)!=dimension(e1))
                               warning("Dimension mismatch of operands in '*'; using trimming/recycling rules.")
                             e2 <- rep(e2, length.out= dimension(e1))   
                             e2 <- if(length(e2)==1) matrix(e2) else diag(e2) 
                            scale(e1) <- e2 %*% scale(e1)
                            return(e1)})   

setMethod("%*%", signature(x="matrix",y="EllipticalDistribution"),
           function(x,y){ if(ncol(x)!=dimension(y))
                               stop("Dimension mismatch of operands in '%*%'.")
                            scale(y) <- x %*% scale(y)
                            return(y)})

