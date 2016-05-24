###############################################################################
# Methods for Discrete Distributions
###############################################################################

## (c) Matthias Kohl: revised P.R. 030707

DiscreteDistribution <- function(supp, prob, .withArith = FALSE,
     .withSim = FALSE, .lowerExact = TRUE, .logExact = FALSE,
     .DistrCollapse =
                  getdistrOption("DistrCollapse"),
     .DistrCollapse.Unique.Warn =
                  getdistrOption("DistrCollapse.Unique.Warn"),
     .DistrResolution = getdistrOption("DistrResolution"),
     Symmetry = NoSymmetry()){
    if(!is.numeric(supp))
        stop("'supp' is no numeric vector")
    if(any(!is.finite(supp)))   # admit +/- Inf?
        stop("infinite or missing values in supp")
    len <- length(supp)
    if(missing(prob)){
        prob <- rep(1/len, len)
    }else{
        if(len != length(prob))
            stop("'supp' and 'prob' must have equal lengths")
        if(any(!is.finite(prob)))
            stop("infinite or missing values in prob")
        if(!identical(all.equal(sum(prob), 1,
                          tolerance = 2*getdistrOption("TruncQuantile")), TRUE))
            stop("sum of 'prob' has to be (approximately) 1")
        if(!all(prob >= 0))
            stop("'prob' contains values < 0")
    }

    o <- order(supp)
    supp <- supp[o]
    prob <- prob[o]
    rm(o)

    if(.DistrCollapse){
       if (len>1 && min(diff(supp))< .DistrResolution){
           erg <- .DistrCollapse(supp, prob, .DistrResolution)
           if (len>length(erg$prob) && .DistrCollapse.Unique.Warn)
               warning("collapsing to unique support values")
           prob <- erg$prob
           supp <- erg$supp
       }
    }else{
       usupp <- unique(supp)
       if(length(usupp) < len){
          if(.DistrCollapse.Unique.Warn)
             warning("collapsing to unique support values")
          prob <- as.vector(tapply(prob, supp, sum))
          supp <- sort(usupp)
          len <- length(supp)
          rm(usupp)
       }
       if(len > 1){
          if(min(diff(supp))< .DistrResolution)
             stop("grid too narrow --> change DistrResolution")
       }
    }
    rm(len)

    if(length(supp) == 1){
      rfun <- function(n){
        rep(supp, n)
      }
    }else{
      rfun <- function(n){
          sample(x = supp, size = n, replace = TRUE, prob = prob)
      }
    }

    dfun <- .makeDNew(supp, prob, Cont = FALSE)
    pfun <- .makePNew(supp, prob, .withSim, Cont = FALSE)
    qfun <- .makeQNew(supp, cumsum(prob), rev(cumsum(rev(prob))),
                      .withSim, min(supp), max(supp), Cont = FALSE)

    object <- new("DiscreteDistribution", r = rfun, d = dfun, q = qfun, p=pfun,
         support = supp, .withArith = .withArith, .withSim = .withSim,
         .lowerExact = .lowerExact, .logExact = .logExact, Symmetry = Symmetry)
}


setMethod("support", "DiscreteDistribution", function(object) object@support)

### left continuous cdf

setMethod("p.l", "DiscreteDistribution", function(object){
       if (.inArgs("lower.tail", p(object))){
           function(q, lower.tail = TRUE, log.p = FALSE){
                px <- p(object)(q, lower.tail = lower.tail)
                o.warn <- getOption("warn");
                on.exit(options(warn=o.warn))
                options(warn = -2)
                dx <- d(object)(.setEqual(q, support(object)))
                options(warn = o.warn)
                px0 <- pmax(px + if(lower.tail) -dx else  dx,0)
                if (log.p) px0 <- log(px0)
                return(px0)
                }
       }else{
           function(q, lower.tail = TRUE, log.p = FALSE){
                px <- p(object)(q)
                o.warn <- getOption("warn")
                on.exit(options(warn=o.warn))
                options(warn = -2)
                dx <- d(object)(.setEqual(q, support(object)))
                options(warn = o.warn)
                px0 <- pmax(if(lower.tail) px - dx else 1 - px + dx, 0)
                if (log.p) px0 <- log(px0)
                return(px0)
                }
       }
})

### right continuous quantile function

setMethod("q.r", "DiscreteDistribution", function(object){
    if (.inArgs("log.p", q(object))){
       if (.inArgs("lower.tail", q(object))){
           function(p, lower.tail = TRUE, log.p = FALSE){
                s <- support(object)
                psx <- p(object)(s, lower.tail = lower.tail,
                                 log.p = log.p)
                ps0 <- .setEqual(p, psx)

                o.warn <- getOption("warn"); options(warn = -2)
                on.exit(options(warn=o.warn))
                qx0 <- q(object)(ps0, lower.tail = lower.tail,
                                 log.p = log.p)
                options(warn = o.warn)

                m <- match(ps0, psx)
                n.ina.m <- !is.na(m)
                if(any(n.ina.m))
                   { M.n.ina.m  <- m[n.ina.m]
                     qx0[n.ina.m] <- (support(object))[pmin(M.n.ina.m+1,
                                                       length(s))]
                   }
                if(any(is.nan(qx0)))
                   warning("NaN's produced")
                return(qx0)
                }
       }else{
           function(p, lower.tail = TRUE, log.p = FALSE){
                s <- support(object)
                psx <- p(object)(s, log.p = log.p)
                if (lower.tail) p <- 1 - p
                ps0 <- .setEqual(p, psx)

                o.warn <- getOption("warn"); options(warn = -2)
                on.exit(options(warn=o.warn))
                qx0 <- q(object)(ps0, lower.tail = lower.tail,
                                 log.p = log.p)
                options(warn = o.warn)

                m <- match(ps0, psx)
                n.ina.m <- !is.na(m)
                if(any(n.ina.m))
                   { M.n.ina.m  <- m[n.ina.m]
                     qx0[n.ina.m] <- (support(object))[pmin(M.n.ina.m+1,
                                                       length(s))]
                   }
                if(any(is.nan(qx0)))
                   warning("NaN's produced")
                return(qx0)
                }
       }
    }else{
       if (.inArgs("lower.tail", q(object))){
           function(p, lower.tail = TRUE, log.p = FALSE){
                if (log.p) p <- exp(p)
                s <- support(object)
                psx <- p(object)(s, lower.tail = lower.tail)
                ps0 <- .setEqual(p, psx)

                o.warn <- getOption("warn"); options(warn = -2)
                on.exit(options(warn=o.warn))
                qx0 <- q(object)(ps0, lower.tail = lower.tail,
                                 log.p = log.p)
                options(warn = o.warn)

                m <- match(ps0, psx)
                n.ina.m <- !is.na(m)
                if(any(n.ina.m))
                   { M.n.ina.m  <- m[n.ina.m]
                     qx0[n.ina.m] <- (support(object))[pmin(M.n.ina.m+1,
                                                       length(s))]
                   }
                if(any(is.nan(qx0)))
                   warning("NaN's produced")
                return(qx0)
                }
       }else{
           function(p, lower.tail = TRUE, log.p = FALSE){
                if (log.p) p <- exp(p)
                s <- support(object)
                psx <- p(object)(s)
                if (lower.tail) p <- 1 - p
                ps0 <- .setEqual(p, psx)

                o.warn <- getOption("warn"); options(warn = -2)
                on.exit(options(warn=o.warn))
                qx0 <- q(object)(ps0, lower.tail = lower.tail,
                                 log.p = log.p)
                options(warn = o.warn)

                m <- match(ps0, psx)
                n.ina.m <- !is.na(m)
                if(any(n.ina.m))
                   { M.n.ina.m  <- m[n.ina.m]
                     qx0[n.ina.m] <- (support(object))[pmin(M.n.ina.m+1,
                                                       length(s))]
                   }
                if(any(is.nan(qx0)))
                   warning("NaN's produced")
                }
       }
    }
})



## Convolution Discrete Distributions

setMethod("+", c("DiscreteDistribution","DiscreteDistribution"),
function(e1,e2){

            if(length(support(e1))==1) return(e2+support(e1))
            if(length(support(e2))==1) return(e1+support(e2))
            e1.L <- as(e1, "LatticeDistribution")
            e2.L <- as(e2, "LatticeDistribution")
            if(is(e1.L, "LatticeDistribution") & is(e2.L, "LatticeDistribution"))
                  {w1 <- width(lattice(e1.L))
                   w2 <- width(lattice(e2.L))
                   W <- sort(abs(c(w1,w2)))
                   if (abs(abs(w1)-abs(w2))<getdistrOption("DistrResolution") ||
                       W[2] %% W[1] < getdistrOption("DistrResolution") )
                       return(e1.L + e2.L)
                  }
            .convDiscrDiscr(e1,e2)})

setMethod("+", c("Dirac","DiscreteDistribution"),
      function(e1,e2){e2+location(e1)})


## binary operators for discrete distributions

setMethod("*", c("DiscreteDistribution","numeric"),
           function(e1, e2) { Distr <- .multm(e1,e2, "DiscreteDistribution")
                              if(is(Distr, "AffLinDistribution"))
                                 Distr@X0 <- e1

                              if(is(e1@Symmetry,"SphericalSymmetry"))
                                 Distr@Symmetry <-
                                   SphericalSymmetry(SymmCenter(e1@Symmetry)*e2)

                              Distr
                             })
setMethod("+", c("DiscreteDistribution","numeric"),
           function(e1, e2) { Distr <- .plusm(e1,e2, "DiscreteDistribution")
                              if(is(Distr, "AffLinDistribution"))
                                 Distr@X0 <- e1

                              if(is(e1@Symmetry,"SphericalSymmetry"))
                                 Distr@Symmetry <-
                                   SphericalSymmetry(SymmCenter(e1@Symmetry)+e2)

                              Distr
                             })

setMethod("*", c("AffLinDiscreteDistribution","numeric"),
           function(e1, e2) {
                Distr <- .multm(e1,e2, "AffLinDiscreteDistribution")
                if(is(e1@Symmetry,"SphericalSymmetry"))
                      Distr@Symmetry <-
                        SphericalSymmetry(SymmCenter(e1@Symmetry)*e2)
                Distr
                })
setMethod("+", c("AffLinDiscreteDistribution","numeric"),
           function(e1, e2) {
                Distr <- .plusm(e1,e2, "AffLinDiscreteDistribution")
                if(is(e1@Symmetry,"SphericalSymmetry"))
                      Distr@Symmetry <-
                        SphericalSymmetry(SymmCenter(e1@Symmetry)*e2)
                Distr
                })

## Group Math for discrete distributions
setMethod("Math", "DiscreteDistribution",
          function(x){
            rnew <- function(n, ...){}
            body(rnew) <- substitute({ f(g(n, ...)) },
                                         list(f = as.name(.Generic), g = x@r))
            object <- new("DiscreteDistribution", r = rnew,
                           .withSim = TRUE, .withArith = TRUE)
            object
          })
setMethod("Math", "Dirac",
          function(x){ loc <- location(x)
                       lc <- callGeneric(loc)
                       Dirac(lc)})

## exact: abs for discrete distributions
setMethod("abs", "DiscreteDistribution",function(x){

       rnew <- function(n, ...){}
       body(rnew) <- substitute({ abs(g(n, ...)) },
                                      list(g = x@r))

       xx <- x
       supportnew <- support(x)

       isSym0 <- FALSE
       if(is(Symmetry(x),"SphericalSymmetry"))
          if(.isEqual(SymmCenter(Symmetry(x)),0))
             isSym0 <- TRUE

       if(isSym0){
          supportnew <- supportnew[supportnew>=0]

          .lowerExact = .lowerExact(x)

          dxlog <- if("log" %in% names(formals(d(x))))
                        quote({dx <- d(xx)(x, log = TRUE)})
                   else quote({dx <-  log(d(xx)(x))})
          pxlog <- if("log.p" %in% names(formals(p(x))) &&
                       "lower.tail" %in% names(formals(p(x))))
                        quote({p(x)(q, lower.tail = FALSE, log.p = TRUE)})
                   else
                        quote({log(1-p(x)(q))})


          qxlog <- if("lower.tail" %in% names(formals(q(x))))
                          quote({qx <- if(lower.tail)
                                          q(x)((1+p1)/2)
                                       else
                                          q(x)(p1/2,lower.tail=FALSE)})
                      else
                          quote({qx <- q(x)(if(lower.tail) (1+p1)/2 else 1-p1/2)})
          if("lower.tail" %in% names(formals(q(x)))&&
             "log.p" %in% names(formals(q(x))))
              qxlog <- quote({qx <- if(lower.tail) q(x)((1+p1)/2)
                                       else
                                          q(x)(if(log.p)p-log(2)
                                               else p1/2,lower.tail=FALSE,log.p=log.p)})


          dnew <- function(x, log = FALSE){}
          body(dnew) <- substitute({
                    dxlog0
                    dx[x>0] <- dx+log(2)
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
            if (.isEqual(p.l(x)(0),0)) return(x)

            supportnew <- sort(unique(abs(supportnew)))

            dnew <- function(x, log = FALSE){
                    o.warn <- getOption("warn"); options(warn = -1)
                    on.exit(options(warn=o.warn))
                    dx <- (x>=0) * d(xx)(x) + (x>0) * d(xx)(-x)
                    options(warn = o.warn)
                    if (log) dx <- log(dx)
                    return(dx)
                 }

            pxlow <- if("lower.tail" %in% names(formals(p(x))))
                        substitute({p(x)(q, lower=FALSE)})
                   else
                        substitute({1-p(x)(q)})

            pnew <- function(q, lower.tail = TRUE, log.p = FALSE){}
            body(pnew) <- substitute({
                    px <- if (lower.tail)
                            (q>=0) * (p(x)(q) - p.l(x)(-q))
                          else pxlow0 + p.l(x)(-q)
                    if (log.p) px <- log(px)
                    return(px)
            }, list(pxlow0=pxlow))

            prob <- dnew(supportnew)

            qnew <- .makeQNew(supportnew, cumsum(prob),
                            rev(cumsum(rev(prob))), notwithLLarg = x@.withSim,
                            min(supportnew), max(supportnew), Cont = FALSE)

         }
         object <- new("DiscreteDistribution", r = rnew, p = pnew,
                        q = qnew, d = dnew, support = supportnew,
                        .withSim = x@.withSim, .withArith = TRUE,
                        .lowerExact = .lowerExact(x))
         object
})

## exact: abs for discrete distributions
setMethod("exp", "DiscreteDistribution",
           function(x) .expm.d(x))


### preliminary to export special functions
if (getRversion()>='2.6.0'){

setMethod("log", "DiscreteDistribution",
           function(x, base = exp(1)) {
           xs <- as.character(deparse(match.call(
                 call = sys.call(sys.parent(1)))$x))
           ep <- getdistrOption("TruncQuantile")
           basl <- log(base)
           if(p(x)(0)>ep)
                stop(gettextf("log(%s) is not well-defined with positive probability ", xs))
           else return(.logm.d(x)/basl)})

setMethod("log", "Dirac",
          function(x, base = exp(1)){
                       xs <- as.character(deparse(match.call(
                             call = sys.call(sys.parent(1)))$x))
                       loc <- location(x)
                       ep <- getdistrOption("TruncQuantile")
                       basl <- log(base)
                       if(loc < ep)
                          stop(gettextf("log(%s) is not well-defined with positive probability ", xs))
                       Dirac(log(loc)/basl)})

setMethod("log10", "DiscreteDistribution",
          function(x) log(x = x, base = 10))

setMethod("sign", "DiscreteDistribution",
          function(x){
          d0 <- d(x)(0)
          DiscreteDistribution(supp=c(-1,0,1),
              prob=c(p(x)(-getdistrOption("TruncQuantile")),
                     d0,
                     p(x)(getdistrOption("TruncQuantile"), lower=FALSE)))
          })


setMethod("digamma", "DiscreteDistribution",
          function(x){
            px0 <- p(x)(0)
            if(px0>0) stop("argument of 'digamma' must be concentrated on positive values")
            rnew <-  function(n, ...){}
            body(rnew) <- substitute({ digamma(g(n, ...)) }, list(g = x@r))

            object <- DiscreteDistribution(
                     supp=digamma(support(x)),
                     prob=prob(x), .withArith = TRUE)
            object
          })

setMethod("lgamma", "DiscreteDistribution",
          function(x){
            rnew = function(n, ...){}
            body(rnew) <- substitute({ lgamma(g(n, ...)) }, list(g = x@r))
            object <- new("DiscreteDistribution", r = rnew,
                           .withSim = TRUE, .withArith = TRUE)
            object
          })

setMethod("gamma", "DiscreteDistribution",
          function(x){
            rnew = function(n, ...){}
            body(rnew) <- substitute({ gamma(g(n, ...)) }, list(g = x@r))
            object <- new("DiscreteDistribution", r = rnew,
                           .withSim = TRUE, .withArith = TRUE)
            object
          })
setMethod("sqrt", "DiscreteDistribution",
            function(x) x^0.5)

}
setMethod("prob", "DiscreteDistribution",
function(object) {sp <- object@support
                  pr <- object@d(sp)
                  names(pr) <- paste(sp)
                  return(pr)
                  })
## Replace Methods
setReplaceMethod("prob", "DiscreteDistribution",
                  function(object, value){
                  return(DiscreteDistribution(supp = object@support,
                             prob = value,
                            .withArith = object@.withArith,
                            .withSim = object@.withSim,
                            .lowerExact = .lowerExact(object),
                            .logExact = .logExact(object)))}
                  )


