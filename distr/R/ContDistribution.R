###############################################################################
# Methods for Absolutely Continuous Distributions
###############################################################################

## (c) P.R. 300408

AbscontDistribution <- function(r = NULL, d = NULL, p = NULL, q = NULL,
                   gaps = NULL, param = NULL, img = new("Reals"),
                   .withSim = FALSE, .withArith = FALSE,
                    .lowerExact = FALSE, .logExact = FALSE,
                   withgaps = getdistrOption("withgaps"),
                   low1 = NULL, up1 = NULL, low = -Inf, up =Inf,
                   withStand = FALSE,
                   ngrid = getdistrOption("DefaultNrGridPoints"),
                   ep = getdistrOption("TruncQuantile"),
                   e = getdistrOption("RtoDPQ.e"),
                   Symmetry = NoSymmetry() 
                  )
{ if(missing(r) && missing(d) && missing(p) && missing(q))
    stop("At least one of arg's r,d,p,q must be given")

  d1 <-  d
  wS <- .withSim
  wA <- .withArith
  if(is.null(r)){
      if(is.null(q)){
          if(is.null(p)){
              if(is.null(low1)){
                  i <- 0; x0 <- -1
                  while(d(x0)> ep && i < 20) x0 <- x0 * 2
                  low1 <- x0
              }
              if(is.null(up1)){
                  i <- 0; x0 <- 1
                  while(d(x0)> ep && i < 20) x0 <- x0 * 2
                  up1 <- x0
              }
              ### new: allow for non standardized d functions i.e. 
              ### possibly, int d(x) dx != 1
              if(withStand){
                  if(is(try(
                     stand <- integrate(d,-Inf,Inf)$value,
                     silent=TRUE),
                     "try-error")){
                     if(is(try(
                         stand <- integrate(d, low1, up1)$value,
                         silent = TRUE),
                         "try-error")){
                          n1 <- 2*ngrid+1
                          n1.odd <- seq(1, n1, by = 2)
                          n1.even <- (1:n1)[-n1.odd]
                          x <- seq(low1,up1, length = n1)
                          h <- diff(x)[1]
                          dx <- d(x)
                          stand <- (4*sum(dx[n1.even])+
                                    2*sum(dx[n1.odd])-
                                    dx[1]-rev(dx)[1])*h/6 
                     }
                  }
                  d0 <- d
                  if(.inArgs("log",d0))
                     d1 <- function(x, log = FALSE){
                                 d00 <- d0(x, log = log)
                                 d00 <- if(log) d00 - log(stand) else d00 / stand
                                 return(d00)
                           }      
                  else
                     d1 <- function(x, log = FALSE){
                                 d00 <- d0(x)/stand
                                 if(log) d00 <- log(d00) 
                                 return(d00)
                           }      
              }
              p <- .D2P(d = d1, ql = low1, qu = up1,  ngrid = ngrid)
              q <- .P2Q(p = p, ql = low1, qu = up1,  ngrid = ngrid,
                        qL = low, qU = up)
              r <- function(n) q(runif(n)) 
          }else{ 
              if(is.null(low1)){
                  i <- 0; x0 <- -1
                  while(p(x0)> ep && i < 20) x0 <- x0 * 2
                  low1 <- x0
              }
              if(is.null(up1)){
                  i <- 0; x0 <- 1
                  while(p(x0)< 1-ep && i < 20) x0 <- x0 * 2
                  up1 <- x0
              }

              q <- .P2Q(p = p, ql = low1, qu = up1,  ngrid = ngrid,
                       qL = low, qU = up)
              r <- function(n) q(runif(n))
              if( is.null(d))
                 d <- .P2D(p = p, ql = low1, qu = up1,  ngrid = ngrid)
          }
      }else{
          if(is.null(p))
             p <- .Q2P(q, ngrid = ngrid)
          r <- function(n) q(runif(n))
          if( is.null(d)){
              if(is.null(low1))
                 low1 <- q(ep)
              if(is.null(up1))
                 up1 <- q(1-ep)
              d <- .P2D(p = p, ql = low1, qu = up1,  ngrid = ngrid)
              }
      }
  }else{
      if(is.null(d)){
          if(is.null(p)){
              if(is.null(q)){
                  erg <- RtoDPQ(r = r, e = e, n = ngrid)
                  wS <- TRUE
                  d <- erg$d; p <- erg$p; q<- erg$q
              }else{
                  p <- .Q2P(q, ngrid = ngrid)
                  if( is.null(d)){
                      if(is.null(low1))
                         low1 <- q(ep)
                      if(is.null(up1))
                         up1 <- q(1-ep)
                      d <- .P2D(p = p, ql = low1, qu = up1,  ngrid = ngrid)
                  }
              }
          }else{
              if(is.null(q)){
                  if(is.null(low1)){
                      i <- 0; x0 <- -1
                      while(p(x0)> ep && i < 20) x0 <- x0 * 2
                      low1 <- x0
                  }
                  if(is.null(up1)){
                      i <- 0; x0 <- 1
                      while(p(x0)< 1-ep && i < 20) x0 <- x0 * 2
                      up1 <- x0
                  }
                  q <- .P2Q(p = p, ql = low1, qu = up1,  ngrid = ngrid,
                           qL = low, qU = up)
                  d <- .P2D(p = p, ql = low1, qu = up1,  ngrid = ngrid)
              }
          }
      }else{
          if(is.null(p)){
              if(is.null(q)){
                  if(is.null(low1)){
                      i <- 0; x0 <- -1
                      while(d(x0)> ep && i < 20) x0 <- x0 * 2
                      low1 <- x0
                      }
                  if(is.null(up1)){
                      i <- 0; x0 <- 1
                      while(d(x0)> ep && i < 20) x0 <- x0 * 2
                      up1 <- x0
                      }

              ### new: allow for non standardized d functions i.e. 
              ### possibly, int d(x) dx != 1
              
                  if(withStand){
                      if(is(try(
                          stand <- integrate(d, -Inf, Inf)$value,
                                silent=TRUE),
                            "try-error")){
                          if(is(try(
                              stand <- integrate(d, low1, up1)$value,
                                    silent = TRUE),
                                "try-error")){
                              n1 <- 2*ngrid+1
                              n1.odd <- seq(1,n1, by=2)
                              n1.even <- (1:n1)[-n1.odd]
                              x <- seq(low1,up1, length = n1)
                              h <- diff(x)[1]
                              dx <- d(x)
                              stand <- (4*sum(dx[n1.even])+
                                        2*sum(dx[n1.odd])-
                                        dx[1]-rev(dx)[1])*h/6 
                          }
                      }
                      d0 <- d
                      if(.inArgs("log",d0))
                         d1 <- function(x, log = FALSE){
                                     d00 <- d0(x, log = log)
                                     d00 <- if(log) d00 - log(stand) else d00 / stand
                                     return(d00)
                               }      
                      else
                         d1 <- function(x, log = FALSE){
                                     d00 <- d0(x)/stand
                                     if(log) d00 <- log(d00) 
                                     return(d00)
                               }      
                  }
        
                  p <- .D2P(d = d1, ql = low1, qu=up1,  ngrid = ngrid)
                  q <- .P2Q(p = p, ql = low1, qu=up1,  ngrid = ngrid,
                            qL = low, qU = up)
              }else
                  p <- .Q2P(q, ngrid = ngrid)
          }else{
              if(is.null(q)){
                  if(is.null(low1)){
                      i <- 0; x0 <- -1
                      while(p(x0)> ep && i < 20) x0 <- x0 * 2
                      low1 <- x0
                  }
                  if(is.null(up1)){
                      i <- 0; x0 <- 1
                      while(p(x0)< 1-ep && i < 20) x0 <- x0 * 2
                      up1 <- x0
                  }
                  q <- .P2Q(p = p, ql = low1, qu=up1,  ngrid = ngrid,
                            qL = low, qU = up)
              }          
          }
      }
  }
  obj <- new("AbscontDistribution", r = r, d = d1, p = p, q = q, 
      gaps = gaps, param = param, img = img, .withSim = wS,
      .withArith = wA, .lowerExact = .lowerExact, .logExact = .logExact,
      Symmetry = Symmetry)

  if(is.null(gaps) && withgaps) setgaps(obj)
  if(!is.null(obj@gaps)) 
     obj@q <- .modifyqgaps(pfun = obj@p, qfun = obj@q, gaps = obj@gaps)
  return(obj)
}

## Access Methods

setMethod("gaps", signature(object = "AbscontDistribution"),  
           function(object) object@gaps)

## ReplaceMethods

setReplaceMethod("gaps", signature(object = "AbscontDistribution"),  
          function(object, value)
                  {dsvalue <- deparse(substitute(value))
                   if(!is.null(value)){
                      if(!is.matrix(value))
                         stop("value must either be a matrix or NULL")
                      if(!ncol(value)==2)
                         stop("if matrix, value must have 2 columns")
                      l <- length(value)
                      if(!identical(1:l, order(c(t(value)))))
                         stop(gettextf("c(t(%s)) must be increasing",
                              dsvalue))
                      colnames(value) <- c("from", "to")
                      }
                   object@gaps <- value; object})
                                                             

setMethod("setgaps", signature(object = "AbscontDistribution"), 
function(object, exactq = 6, ngrid = 50000, ...){
       object1 <- object
       lower <- getLow(object, eps = getdistrOption("TruncQuantile")*2)
       upper <- getUp(object, eps = getdistrOption("TruncQuantile")*2)
       #lower <- 0 ; upper <- 8
       dist <- upper - lower
       low1 <- max(q(object)(0),lower-0.1*dist)
       upp1 <- min(q(object)(1),upper+0.1*dist)
       grid <- seq(from = low1, to = upp1, length = ngrid) 
       dxg <- d(object)(grid)
       
       ix <-  1:ngrid
       ix0 <- (dxg < 1/10^exactq)&(grid>=lower)&(grid<=upper)     
       if(any(ix0)){
          ixc <- ix[ix0]
          dixc <- c(2,diff(ixc))
          dixc0 <- dixc>1
          l2 <- length(ixc)
          ixc2 <- seq(l2)
          ixcl <- ixc2[dixc0]
          ixcr <- c(ixcl[-1]-1,l2)
          gridl <- grid[ixc[ixcl]]
          gridr <- grid[ixc[ixcr]]
          
          mattab.d <- cbind(gridl, gridr)
          
          ox <- order(mattab.d[,1])
          mattab.d <- matrix(mattab.d[ox,], ncol = 2)
          mattab.d <- .consolidategaps(mattab.d)
          if(nrow(mattab.d)==0) mattab.d <- NULL
          if(length(mattab.d)==0) mattab.d <- NULL
          } else mattab.d <- NULL
          eval(substitute( "slot<-"(object,'gaps', value = mattab.d)))
       return(invisible())
})
 
## Arithmetics

setMethod("+", c("AbscontDistribution","AbscontDistribution"),
function(e1,e2){
            ### Step 1 : Truncation
            
            lower <- min(getLow(e1), getLow(e2))
            upper <- max(getUp(e1) , getUp(e2))

            ### Step 2 : Discretizing

            n <- getdistrOption("DefaultNrFFTGridPointsExponent")
            h <- (upper-lower)/2^n

            dpe1 <- .discretizeP(e1, lower, upper, h)
            dpe2 <- .discretizeP(e2, lower, upper, h)

            x <- seq(from = 2*lower, to = 2*upper, by = h)

            ### Step 3 : Zero-Padding

            dpe1 <- c(dpe1, numeric(2^n))
            dpe2 <- c(dpe2, numeric(2^n))

            ## Step 4: computation of DFT

            ftpe1 <- fft(dpe1); ftpe2 <- fft(dpe2)
            ## convolution theorem for DFTs
            d2 <- c(0,Re(fft(ftpe1*ftpe2, inverse = TRUE)) / length(ftpe1))

            ## density & cdf (steps 5--7)
            dfun <- .makeDNew(x, d2, h)
            pfun <- .makePNew(x, d2, h, .notwithLArg(e1)||.notwithLArg(e2) )


            ## quantile function
            yL <-  if ((q(e1)(0) == -Inf)||(q(e2)(0) == -Inf))
                 -Inf else getLow(e1)+getLow(e2)
            yR <-  if ((q(e1)(1) ==  Inf)||(q(e2)(1) ==  Inf))
                  Inf else getUp(e1)+getUp(e2)

            px.l <- pfun(x + 0.5*h)
            px.u <- pfun(x + 0.5*h, lower.tail = FALSE)
            
            qfun <- .makeQNew(x + 0.5*h, px.l, px.u,
                              .notwithLArg(e1)||.notwithLArg(e2), yL, yR)
            
            rfun <- function(n){}
            body(rfun) <- substitute({ f(n) + g(n) },
                                     list(f = e1@r, g = e2@r))

            object <- AbscontDistribution(r = rfun, d = dfun, p = pfun,
                          q = qfun, .withSim = FALSE, .withArith = TRUE)

            rm(d2, dpe1,dpe2, ftpe1,ftpe2)
            rm(h, px.l, px.u, rfun, dfun, qfun, pfun, upper, lower)
            if(is(e1@Symmetry,"SphericalSymmetry")&& 
               is(e2@Symmetry,"SphericalSymmetry"))
               object@Symmetry <- SphericalSymmetry(SymmCenter(e1@Symmetry)+
                                                     SymmCenter(e2@Symmetry))   
            object
          })


###setMethod("m1df", "AbscontDistribution",
###   function(object){
###     lower <- q(object)(TruncQuantile)
###     upper <- q(object)(1 - TruncQuantile)
###     
###     gitter.x <- seq(from = lower, to = upper, length = DefaultNrGridPoints)
###     
###    integrand <- function(x) x * d(object)(x)
###     
###     tmp <- function(t) integrate(integrand, lower = lower, upper = t)$value
###     
###     gitter.y <- sapply(gitter.x, tmp)
###     
###     approxfun(gitter.x, gitter.y, rule = 2)
###   })


###setMethod("m2df", "AbscontDistribution", 
###   function(object){
###     lower <- q(object)(TruncQuantile)
###     upper <- q(object)(1 - TruncQuantile)
###     
###     gitter.x <- seq(from = lower, to = upper, length = DefaultNrGridPoints)
###     
###     integrand <- function(x) x^2 * d(object)(x)
###     
###     tmp <- function(t) integrate(integrand, lower = lower, upper = t)$value
###     
###     gitter.y <- sapply(gitter.x, tmp)
###     
###     approxfun(gitter.x, gitter.y, rule = 2)
###   })

## binary operators for absolut continuous distributions


setMethod("*", c("AbscontDistribution","numeric"),
          function(e1, e2) {Distr <-  .multm(e1,e2, "AbscontDistribution")                               
                            if(is(Distr, "AffLinDistribution"))
                                 Distr@X0 <- e1
                            if(is(e1@Symmetry,"SphericalSymmetry"))
                               Distr@Symmetry <- 
                                 SphericalSymmetry(SymmCenter(e1@Symmetry)*e2)

                            Distr})
setMethod("+", c("AbscontDistribution","numeric"),
           function(e1, e2) {Distr <-  .plusm(e1,e2, "AbscontDistribution")                               
                            if(is(Distr, "AffLinDistribution"))
                                 Distr@X0 <- e1
                            if(is(e1@Symmetry,"SphericalSymmetry"))
                               Distr@Symmetry <- 
                                 SphericalSymmetry(SymmCenter(e1@Symmetry)+e2)
                            Distr})                            
setMethod("*", c("AffLinAbscontDistribution","numeric"),
          function(e1, e2){Distr <-  .multm(e1,e2, "AffLinAbscontDistribution")
                           if(is(e1@Symmetry,"SphericalSymmetry"))
                               Distr@Symmetry <- 
                                 SphericalSymmetry(SymmCenter(e1@Symmetry)*e2)
                           Distr                         })
setMethod("+", c("AffLinAbscontDistribution","numeric"),
           function(e1, e2){Distr <-  .plusm(e1,e2, "AffLinAbscontDistribution")
                            if(is(e1@Symmetry,"SphericalSymmetry"))
                               Distr@Symmetry <- 
                                 SphericalSymmetry(SymmCenter(e1@Symmetry)+e2)
                            Distr                         })

## Group Math for absolutly continuous distributions
setMethod("Math", "AbscontDistribution",
          function(x){            

            rnew <- function(n, ...){}           
            body(rnew) <- substitute({ f(g(n, ...)) },
                              list(f = as.name(.Generic), g = x@r))
            
            n <- 10^getdistrOption("RtoDPQ.e")+1
            u <- seq(0,1,length=n+1); u <- (u[1:n]+u[2:(n+1)])/2
            y <- callGeneric(q(x)(u))
            DPQnew <- RtoDPQ(r=rnew, y=y)
                       
            object <- AbscontDistribution(d = DPQnew$d, p = DPQnew$p, 
                           r = rnew, q = DPQnew$q,
                           .withSim = TRUE, .withArith = TRUE)
            object
          })

## exact: abs for absolutly continuous distributions
setMethod("abs", "AbscontDistribution",
    function(x){
       if (.isEqual(p(x)(0),0)) return(x)
       xx <- x
       rnew <- function(n, ...){}
       body(rnew) <- substitute({ abs(g(n, ...)) }, list(g = xx@r))
       
       isSym0 <- FALSE
       if(is(Symmetry(xx),"SphericalSymmetry"))
          if(.isEqual(SymmCenter(Symmetry(xx)),0))
             isSym0 <- TRUE  
       
       if(isSym0){
          if (is.null(gaps(xx)))
              gapsnew <- NULL
          else {gapsnew <- gaps[gaps[,2]>=0,]
                VZW <- gapsnew[,1] <= 0 
                gapsnew[VZW,1] <- 0
                gapsnew <- .consolidategaps(gapsnew)}
          dOx <- d(xx)

          dxlog <- if("log" %in% names(formals(dOx))) 
                        quote({dOx(x, log = TRUE)})
                   else quote({log(dOx(x))})
          pxlog <- if("log.p" %in% names(formals(p(x))) && 
                       "lower.tail" %in% names(formals(p(x)))) 
                        quote({p(xx)(q, lower.tail = FALSE, log.p = TRUE)})
                   else
                        quote({log(1-p(xx)(q))})

          qxlog <- if("lower.tail" %in% names(formals(q(xx)))) 
                          quote({qx <- if(lower.tail)
                                          q(xx)((1+p1)/2)
                                       else
                                          q(xx)(p1/2,lower.tail=FALSE)}) 
                      else
                          quote({qx <- q(xx)(if(lower.tail) (1+p1)/2 else 1-p1/2)})
          if("lower.tail" %in% names(formals(q(xx)))&& 
             "log.p" %in% names(formals(q(xx))))           
              qxlog <- quote({qx <- if(lower.tail) q(xx)((1+p1)/2)
                                       else
                                          q(xx)(if(log.p)p-log(2)
                                               else p1/2,lower.tail=FALSE,log.p=log.p)}) 
          dnew <- function(x, log = FALSE){}
          body(dnew) <- substitute({
                    dx <- (dxlog0 + log(2))*(x>=0)
                    if (!log) dx <- exp(dx)
                    dx[x<0] <- if(log) -Inf else 0
                    return(dx)
                    }, list(dxlog0 = dxlog))
            
          pnew <- function(q, lower.tail = TRUE, log.p = FALSE){}
          body(pnew) <- substitute({
                    if (!lower.tail){
                        px <- (log(2) + pxlog0)*(q>=0)
                        if(!log.p) px <- exp(px)
                    }else{
                        px <- pmax(2 * p(x)(q) - 1,0)
                        if(log.p) px <- log(px)
                    }
                    return(px)            
            }, list(pxlog0 = pxlog))

          qnew <- function(p, lower.tail = TRUE, log.p = FALSE){}
          body(qnew) <- substitute({
                   p1 <- if(log.p) exp(p) else p 
                   qxlog0
                   qx[p1<0] <- NaN
                   if (any((p1 < -.Machine$double.eps)|(p1 > 1+.Machine$double.eps)))
                   warning(gettextf("q method of %s produced NaN's ", objN))
                   return(qx)
            }, list(qxlog0 = qxlog, objN= quote(.getObjName(1))))
                   
       }else{
            if (is.null(gaps(xx)))
                gapsnew <- NULL
            else {VZW <- gaps(xx)[,1] <= 0 & gaps(xx)[,2] >= 0
                  gapsnew <- t(apply(abs(gaps(xx)), 1, sort))
                  gapsnew[VZW,2] <- pmin(-gaps(xx)[VZW,1], gaps(x)[VZW,2])
                  gapsnew[VZW,1] <- 0
                  gapsnew <- .consolidategaps(gapsnew)}
            
            lower <- max(0, getLow(xx))
            upper <- max(-getLow(xx) , abs(getUp(xx)))

            n <- getdistrOption("DefaultNrFFTGridPointsExponent")
            h <- (upper-lower)/2^n

            x.g <- seq(from = lower, to = upper, by = h)

            dnew <- function(x, log = FALSE){
                    o.warn <- getOption("warn"); 
                    on.exit(options(warn=o.warn))
                    options(warn = -1)
                    dx <- (x>=0) * (d(xx)(x) + d(xx)(-x))                     
                    if (log) dx <- log(dx)
                    return(dx)
            }
            
            pxlow <- if("lower.tail" %in% names(formals(p(xx))))
                        substitute({p(xx)(q, lower=FALSE)})
                   else
                        substitute({1-p(xx)(q)})

            pnew <- function(q, lower.tail = TRUE, log.p = FALSE){}
            body(pnew) <- substitute({
                    px <- if (lower.tail)
                            (q>=0) * (p(xx)(q) - p(xx)(-q))                    
                          else pxlow0 + p(xx)(-q)
                    if (log.p) px <- log(px)
                    return(px)
            }, list(pxlow0 = pxlow))

            px.l <- pnew(x.g + 0.5*h)
            px.u <- pnew(x.g + 0.5*h, lower.tail = FALSE)
            
            yR <- max(q(xx)(1), abs(q(xx)(0)))

            qnew <- .makeQNew(x.g + 0.5*h, px.l, px.u,
                              notwithLLarg = FALSE,  lower, yR)
    
            lowerExact <- FALSE
 
    }
    object <- AbscontDistribution( r = rnew, p = pnew, q = qnew, d = dnew, 
                     gaps = gapsnew,  .withSim = xx@.withSim, .withArith = TRUE,
                     .lowerExact = .lowerExact(x), .logExact = FALSE)
    object
    })

## exact: exp for absolutly continuous distributions
setMethod("exp", "AbscontDistribution",
           function(x) .expm.c(x))


### preliminary to export special functions
if (getRversion()>='2.6.0'){ 

setMethod("log", "AbscontDistribution",
           function(x, base = exp(1)) {
           xs <- as.character(deparse(match.call(
                 call = sys.call(sys.parent(1)))$x))
           ep <- getdistrOption("TruncQuantile")
           basl <- log(base)
           if(p(x)(0)>ep) 
                stop(gettextf("log(%s) is not well-defined with positive probability ", xs))
           else return(.logm.c(x)/basl)})
                       
                       
setMethod("log10", "AbscontDistribution",
          function(x) log(x = x, base = 10))

setMethod("sign", "AbscontDistribution",
          function(x){ 
       DiscreteDistribution(supp=c(-1,0,1), 
              prob=c(p(x)(0), 0, p(x)(0, lower=FALSE)))                     
          })



setMethod("digamma", "AbscontDistribution",
          function(x){
            rnew <-  function(n, ...){}
            body(rnew) <- substitute({ digamma(g(n, ...)) }, list(g = x@r))
            px0 <- p(x)(0)
            if(px0>0) stop("argument of 'digamma' must be concentrated on positive values")
            xx <- x
                    
            pnew <- function(q, lower.tail = TRUE, log.p = FALSE){
                    iq <- igamma(q) 
                    px <- p(xx)(iq, lower.tail = lower.tail, log.p = log.p)
                    return(px)
            }
            dnew <- function(x, log = FALSE){
                    ix <- igamma(x)
                    dx <- d(xx)(ix, log = log)
                    nx <- trigamma(ix)
                    if(log) dx <- dx - log(nx)
                    else dx <- dx/nx
                    return(dx)
            }
            
            .x <- sort(c(qexp(unique(pmin(seq(0,1,length=5e4)+1e-10,1-1e-10))),
                       -abs(rnorm(1e4)),
                       qcauchy(seq(0.999,1-1e-10,length=5e3),lower.tail=FALSE)))
            i <- 0; x0 <- 1
            while(pnew(x0,lower.tail = FALSE)>  getdistrOption("TruncQuantile") && i < 20) 
                 x0 <- x0 * 2
             up1 <- x0
            i <- 0; x0 <- -1
            while(pnew(x0)> getdistrOption("TruncQuantile") && i < 20) 
                 x0 <- x0 * 2
             low1 <- x0
          

            
            qnew <- .P2Q(p = pnew, xx =.x,
                         ql = low1, qu=up1,  
                         ngrid = getdistrOption("DefaultNrGridPoints"),
                            qL = -Inf, qU = Inf)
 
            
            object <- AbscontDistribution( r = rnew, d = dnew, p = pnew, q=qnew,
                           .withSim = TRUE, .withArith = TRUE, .logExact = FALSE)
            object
          })

setMethod("lgamma", "AbscontDistribution",
          function(x){
            rnew <- function(n, ...){}
            body(rnew) <- substitute({ lgamma(g(n, ...)) }, list(g = x@r))

            n <- 10^getdistrOption("RtoDPQ.e")+1
            u <- seq(0,1,length=n+1); u <- (u[1:n]+u[2:(n+1)])/2
            y <- lgamma(q(x)(u))
            DPQnew <- RtoDPQ(r=rnew, y=y)
            
            object <- AbscontDistribution( r = rnew, d = DPQnew$d, p = DPQnew$p,
                           q=DPQnew$q, .withSim = TRUE, .withArith = TRUE)
            object
          })

setMethod("gamma", "AbscontDistribution",
          function(x){
            rnew <- function(n, ...){}
            body(rnew) <- substitute({ gamma(g(n, ...)) }, list(g = x@r))
            n <- 10^getdistrOption("RtoDPQ.e")+1
            u <- seq(0,1,length=n+1); u <- (u[1:n]+u[2:(n+1)])/2
            y <- gamma(q(x)(u))
            DPQnew <- RtoDPQ(r=rnew, y=y)

            object <- AbscontDistribution( r = rnew, d = DPQnew$d, p = DPQnew$p,
                           q=DPQnew$q, .withSim = TRUE, .withArith = TRUE)
            object
          })
          
setMethod("sqrt", "AbscontDistribution",
            function(x) x^0.5)

}

#------------------------------------------------------------------------
# new p.l, q.r methods
#------------------------------------------------------------------------

setMethod("p.l", signature(object = "AbscontDistribution"),  
           function(object) p(object))

setMethod("q.r", signature(object = "AbscontDistribution"),  
           function(object){
                if(!is.null(gaps(object))) 
                   .modifyqgaps(pfun = p(object), qfun = q(object), 
                                gaps = gaps(object), leftright = "right")
                else
                    q(object)
            })
