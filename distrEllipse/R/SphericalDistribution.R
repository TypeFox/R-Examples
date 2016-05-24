## Generating function

SphericalDistribution <- function(radDistr = sqrt(Chisq(df=dim)), dim = 2,
                                  p = NULL, q = NULL){
   ell <- EllipticalDistribution(radDistr = radDistr,
               loc = numeric(dim), scale = diag(dim), p = p, q = q)
   sph <- as(ell, "SphericalDistribution")
   sph@Symmetry <- SphericalSymmetry(SymmCenter(Symmetry(ell)))
   sph
   }

## accessors

setMethod("dimension", "SphericalDistribution",
           function(object) dimension(object@img))
setMethod("dim", "SphericalDistribution",
           function(x) dimension(x@img))

setMethod("radDistr", "SphericalDistribution",
           function(object) object@radDistr)

setMethod("scale", "SphericalDistribution",
           function(x,  center, scale) as.matrix(diag(dimension(x))))

setMethod("location", "SphericalDistribution",
           function(object) numeric(dimension(object)))

#not necessary:
##setMethod("Symmetry", "SphericalDistribution", function(object) object@Symmetry)

## replacements

setReplaceMethod("radDistr", "SphericalDistribution",
           function(object,value){
              if(!is(value,"UnivariateDistribution"))
                  stop("RHS must be a univariate Distribution")
              if(p.l(value)(0)>0)
                 stop("RHS must have pos. support")
              object@radDistr <- value
              object }
           )

## wrappers:

setMethod("plot.rd", "SphericalDistribution",
           function(x, ... ) plot(x@radDistr,...))
setMethod("r.rd", "SphericalDistribution", function(object) r(object@radDistr))
setMethod("d.rd", "SphericalDistribution", function(object) d(object@radDistr))
setMethod("p.rd", "SphericalDistribution", function(object) p(object@radDistr))
setMethod("q.rd", "SphericalDistribution", function(object) q(object@radDistr))

## functionals:

setMethod("E", signature(object = "SphericalDistribution",
                        fun = "missing", cond = "missing"),
           function(object,...) numeric(dimension(object)))

setMethod("var", signature(x = "SphericalDistribution"),
           function(x,...) diag(dimension(x)) *
                    E(radDistr(x), fun = function(y)y^2,...)/dimension(x)
           )


setAs("SphericalDistribution", "EllipticalDistribution",
            function(from){
               sc <- SymmCenter(Symmetry(from))
               slotNames <- slotNames(from)
               lst <- sapply(slotNames, function(x) slot(from,x))
               names(lst) <- slotNames
               dim <- dimension(from)
               lst$"param" <- new("EllipticalParameter",
                                   loc = numeric(dim),
                                   scale = diag(dim))
               lst$"Symmetry" <- EllipticalSymmetry(numeric(dim))
               ell <- new("EllipticalDistribution")
               for (i in 1: length(lst))
                   slot(ell, name = names(lst)[i]) <- lst[[i]]
               ell})
               

setMethod("plot", signature(x = "SphericalDistribution", y = "missing"),
      function(x, Nsim = getdistrEllipseOption("Nsim"), ...,
               withED = getdistrEllipseOption("withED"),
               lwd.Ed = getdistrEllipseOption("lwd.Ed"),
               col.Ed = getdistrEllipseOption("col.Ed"),
               withMean = getdistrEllipseOption("withMean"),
               cex.mean = getdistrEllipseOption("cex.mean"),
               pch.mean = getdistrEllipseOption("pch.mean"),
               col.mean = getdistrEllipseOption("col.mean")){
      dots <- match.call(call = sys.call(sys.parent(1)),
                         expand.dots = FALSE)$"..."
      cex <- 0.5
      if(hasArg(cex)) cex <- dots$cex
      col <- "black"
      if(hasArg(col)) col <- dots$col
      
      qchs <- qchisq(.95, df = 2)^.5
      col.Ed <- rep(col.Ed, length.out = 2)

      X <- r(x)(2000)

      if(hasArg(panel))
         pairs(t(X), ...)
      else
         pairs(t(X), ...,
               panel = function(x,y, cex = cex, col = col, ...){
                  dots$col <- NULL
                  dots$cex <- NULL
                  do.call(points, c(list(x = x, y = y), dots))
                  if(withED){
                     co <- var(cbind(x,y))
                     eig <- eigen(co)
                     ev <- eig$values^.5
                     x1 <- eig$vectors[,1]*ev[1]*qchs
                     x2 <- eig$vectors[,2]*ev[2]*qchs
                     
                     lines(x = c(0,x1[1]) + mean(x),
                           y = c(0,x1[2]) + mean(y), lwd = lwd.Ed,
                        col = col.Ed[1])
                     lines( x = c(0,x2[1]) + mean(x),
                         y = c(0,x2[2]) + mean(y), lwd = lwd.Ed,
                         col = col.Ed[2])
                     }
                     if(withMean)
                     points(mean(x),mean(y), col = col.mean,
                         cex = cex.mean , pch = pch.mean)
                  } )
      return(invisible(NULL))
      } )


setMethod("+", c("SphericalDistribution","numeric"),
           function(e1,e2) as(e1, "EllipticalDistribution")+e2)
setMethod("*", c("SphericalDistribution","numeric"),
           function(e1,e2)as(e1, "EllipticalDistribution")*e2)
setMethod("%*%", signature(x="matrix",y="SphericalDistribution"),
           function(x,y) x %*% as(y, "EllipticalDistribution"))

setMethod("+", c("numeric", "SphericalDistribution"),
          function(e1, e2){
            e2 + e1
          })


setMethod("*", c("numeric", "SphericalDistribution"),
          function(e1, e2){
            e2 * e1
          })

setMethod("-", c("SphericalDistribution", "missing"),
          function(e1){
            e1*(-1)
          })


setMethod("-", c("SphericalDistribution", "numeric"),
          function(e1, e2){
            return(e1 + (-e2))
          })

setMethod("-", c("numeric", "SphericalDistribution"),
          function(e1, e2){
            -1*e2 + e1
          })
