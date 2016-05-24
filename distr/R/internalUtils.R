#------------------------------------------------------------------------------
### internal help function to check whether numeric vectors are integers
#------------------------------------------------------------------------------

.isInteger  <- function(x, tol = .Machine$double.eps) abs(as.integer(x)-x)< tol
.isNatural  <- function(x, tol = .Machine$double.eps) .isInteger(x, tol) & (x>0)
.isNatural0 <- function(x, tol = .Machine$double.eps) .isInteger(x, tol) & (x>=0)
setAs("numeric","Integer",function(from) new("Integer",as.integer(from)))

#------------------------------------------------------------------------------
### internal help function to check if a vector can be made a lattice
#------------------------------------------------------------------------------

.is.vector.lattice <- function(x)
  {    ### is x equally spaced?
    all( sapply(diff(x), function(y)
         isTRUE(all.equal(y, diff(x)[1],
                tolerance = getdistrOption("DistrResolution"),
                check.attributes = FALSE))))
  }

### internal help function to check consistency of lattice with support:
.is.consistent <- function(lattice, support, eq.space = TRUE)
  {
   p <- pivot(lattice); l <- Length(lattice); w <- width(lattice)
   ds <- diff(support); ms <- min(support); Ms <- max(support)
   ### is support equally spaced?
   if (! .is.vector.lattice(support)  && eq.space)
      return(FALSE)
   ### are width of lattice and support consistent
   if (! isTRUE(all.equal(min(ds), abs(w), check.attributes = FALSE)))
      return(FALSE)
   ### pivot is left or right endpoint of support
   if ( isTRUE(all.equal(ms, p, check.attributes = FALSE)) || isTRUE(all.equal(Ms, p, check.attributes = FALSE)) )
      return(TRUE)

   if (isTRUE(all.equal(min((support[1]-p)%%w,w-(support[1]-p)%%w),0,
                        tolerance = getdistrOption("TruncQuantile"),
                        check.attributes = FALSE)))
      return(TRUE)
  return(FALSE)
  }

 ## generate a lattice from equally spaced vector
.make.lattice.es.vector <- function(x){
  new("Lattice", pivot = x[1], width = x[2]-x[1],
       Length = length(x))
}

#------------------------------------------------------------------------------
### internal help function to determin common lattice width
#------------------------------------------------------------------------------

.EuclidAlgo <- function(n1,n2){
   r<- 2
   m <- round(max(n1,n2))
   n <- round(min(n1,n2))
   if (n==1) return(1)
   while (r>1){
      r <- m %% n
      m <- n
      if(r==1) return(1)
      if(r==0) break
      n <- r
      }
   return(n)
}
.getCommonWidth <- function(x1,x2, tol=.Machine$double.eps){
  rI <- function(x) .isInteger(x, tol = x * tol)
  if(rI(x1)&&rI(x2)) return(.EuclidAlgo(x1,x2))
  if(rI(10000*x1) && rI(10000*x2)) return(.EuclidAlgo(10000*x1,10000*x2)/10000)
  n <- max(x1,x2); m <- min(x1,x2)
  if(rI(n/m)) return(m)
  vc <- 1:10000
  vecT <- sapply(vc,function(x) rI(x*(n/m)))
  vec <- m/vc
#  print(min(vc[vecT]))
  if(any(vecT)) return((vec/min(vc[vecT]))[1])
  return(NULL)  
}

#------------------------------------------------------------------------------
### %in% for numerics with tolerance
#------------------------------------------------------------------------------
.inWithTol <- function(x,y,tol=.Machine$double.eps){
   sapply(x, function(u) any(abs(u-y)<tol))
}                                                                         

#------------------------------------------------------------------------------
### brute force convolution
#------------------------------------------------------------------------------
.convDiscrDiscr <- function(e1,e2){
            convolutedsupport <- rep(support(e1), each = length(support(e2))) +
                                 support(e2)

            gridvalues1 <- d(e1)(support(e1)); gridvalues2 <- d(e2)(support(e2))
            convolutedvalues <- rep(gridvalues1, each = length(support(e2))) *
                                gridvalues2
            rm(gridvalues1,gridvalues2)

            tmptable <- data.frame(x = convolutedsupport, dx = convolutedvalues)
            rm(convolutedsupport,convolutedvalues)
            tmp <- tapply(tmptable$dx, tmptable$x, sum)
            rm(tmptable)

            supp.u <- as.numeric(names(tmp))
            prob.u <- as.numeric(tmp)

            o <- order(supp.u)
            supp <- supp.u[o]
            prob <- prob.u[o]

            #supp.u <- unique(supp)

            len <- length(supp)

            if(len > 1){
               if (min(diff(supp))< getdistrOption("DistrResolution")){
                   if (getdistrOption("DistrCollapse")){
                       erg <- .DistrCollapse(supp, prob, 
                                   getdistrOption("DistrResolution"))
                       if ( len > length(erg$prob) && 
                                getdistrOption("DistrCollapse.Unique.Warn") )
                            warning("collapsing to unique support values")         
                       prob <- erg$prob
                       supp <- erg$supp
                   }else
                    stop("grid too narrow --> change DistrResolution")
                }
            }

            rm(tmp, len)

            .withSim <- e1@.withSim || e2@.withSim

            rfun <- function(n) {}
            body(rfun) <- substitute({ f(n) + g(n) },
                                         list(f = e1@r, g = e2@r))

            dfun <- .makeDNew(supp, prob, Cont = FALSE)
            pfun <- .makePNew(supp, prob, .withSim, Cont = FALSE)
            qfun <- .makeQNew(supp, cumsum(prob), rev(cumsum(rev(prob))),
                      .withSim, min(supp), max(supp), Cont = FALSE)

            object <- new("DiscreteDistribution", r = rfun, d = dfun, p = pfun,
                           q = qfun, support = supp,
                           .withSim = .withSim, .withArith = TRUE)
            rm(rfun, dfun, qfun, pfun)

            if(is(e1@Symmetry,"SphericalSymmetry")&& 
               is(e2@Symmetry,"SphericalSymmetry"))
               object@Symmetry <- SphericalSymmetry(SymmCenter(e1@Symmetry)+
                                                     SymmCenter(e2@Symmetry))   

            object
          }

#------------------------------------------------------------------------------
### .fm, .fM, .fM2 functions
#------------------------------------------------------------------------------

.fM <- function(x,f){
   owarn <- getOption("warn")
   options("warn"=-1)
   on.exit(options("warn"=owarn))
   if(!.inArgs("log.p", f))
      f0 <- f
   else f0 <- function(x) f(log(x), log.p = TRUE)   
   xo <-  x
   x1 <- (1+xo)/2
   i <- 1
   while( i < 30)
   { i <- i+1
     while(!is.na(f1 <- f0(x1)) && f1 < Inf)
        {xo <- x1
         x1 <- (1+x1)/2
         }
     x1 <- (x1+xo)/2}
   log(xo)}
   
.fM2 <- function(x,f){
   owarn <- getOption("warn")
   options("warn"=-1)
   on.exit(options("warn"=owarn))
   if(!.inArgs("log.p", f))
      f0 <- function(x) f(x, lower.tail = FALSE)
   else f0 <- function(x) f(log(x), lower.tail = FALSE, log.p = TRUE)   
   xo <- x
   x1 <- xo/2
   i <- 1
   while( i < 30)
   { i <- i+1
     while(!is.na(f1 <- f0(x1)) && f1 < Inf)
        {xo <- x1
         x1 <- x1/2
         }
     x1 <- (x1+xo)/2}   
   log(xo)}

.fm <- function(x,f){
   owarn <- getOption("warn")
   options("warn"=-1)
   on.exit(options("warn"=owarn))
   if(!.inArgs("log.p", f))
      f0 <- f
   else f0 <- function(x) f(log(x), log.p = TRUE)   
   
   xo <- x
   x1 <- xo/2
   i <- 1
   while( i < 30)
   { i <- i+1
     while(!is.na(f1 <- f0(x1)) && f1 > -Inf)
        {xo <- x1
         x1 <- x1/2
         }
     x1 <- (x1+xo)/2}
   log(xo)}
   
#------------------------------------------------------------------------------
# .presubs : for titles etc
#------------------------------------------------------------------------------

.presubs <- function(inp, frompat, topat){
### replaces in an expression or a string all frompat patterns to topat patterns

logic <- FALSE
inCx <- sapply(inp,
   function(inpx){
      inC <- deparse(inpx)
      l <- length(frompat)
      for(i in 1:l)
         { if (is.language(topat[[i]])){
               totxt <- deparse(topat[[i]])
               totxt <- gsub("expression\\(", "\", ", gsub("\\)$",", \"",totxt))
               if (length(grep(frompat[i],inC))) logic <<- TRUE
               inC <- gsub(frompat[i],totxt,inC)
           }else inC <- gsub(frompat[i], topat[[i]], inC)
         }
      return(inC)
    })
if(length(grep("expression",inCx))>0)
   inCx <- gsub("expression\\(", "", gsub("\\)$","",inCx))
if (length(inCx) > 1) {
   inCx <- paste(inCx, c(rep(",", length(inCx)-1), ""),
                 sep = "", collapse = "\"\\n\",")
   if ( any(as.logical(c(lapply(inp,is.language)))) | logic )
      inCx <- paste("expression(paste(", gsub("\\\\n"," ", inCx), "))", sep ="")
   else
      inCx <- paste("paste(",inCx,")", sep ="")
}else inCx <- paste("expression(paste(",inCx,"))",sep="")
outC <- eval(parse(text = eval(inCx)))
return(outC)
}

#------------------------------------------------------------------------------
# help check functions
#------------------------------------------------------------------------------

.inArgs <- function(arg, fct)
          {as.character(arg) %in% names(formals(fct))}

.isEqual <- function(p0, p1, tol = min( getdistrOption("TruncQuantile")/2,
                                          .Machine$double.eps^.7
                                          ))
                abs(p0-p1)< tol

.isIn <- function(p0, pmat, tol = min( getdistrOption("TruncQuantile")/2,
                                          .Machine$double.eps^.7
                                          ))
                  {list1 <- lapply(1:nrow(pmat), function(x){ 
                            (p0+tol > pmat[x,1]) & (p0-tol < pmat[x,2]) })
                   apply(matrix(unlist(list1), ncol = nrow(pmat)), 1, any)}           


.isEqual01<- function(x) .isEqual(x,0)|.isEqual(x,1)

.setEqual <- function(x, y, tol = 1e-7){ 
### all elements of x equal to some element of y up tol are set to exactly
###     the respective element of y
   x1 <- round(2*x/tol,0)
   y1 <- round(2*y/tol,0)
   z  <- x
   m  <- match(x1,y1)
   n.ina.m <- !is.na(m)
   z[n.ina.m] <- y[m[n.ina.m]]
   z
}


.notwithLArg <- function(D)  D@.withSim||!.inArgs("lower.tail",p(D))

#------------------------------------------------------------------------------
# other helpers
#------------------------------------------------------------------------------

.getObjName <- function(i = 1)
     {Ca <- sys.call(-4);  as.character(as.list(Ca[[1]]))[i+1]}

.discretizeP <- function(D, lower, upper, h){
   h0 <- 40*(getUp(D)-getLow(D)) /
         2^getdistrOption("DefaultNrFFTGridPointsExponent")
   if(h > h0 )
      warning(paste("Grid for approxfun too wide, ",
                     "increase DefaultNrFFTGridPointsExponent", sep =""))

   x <- seq(from = lower, to = upper, by = h)
   if(TRUE){#.notwithLArg(D)){
      return(diff(p(D)(x)))
#      return((diff(p(D)(x))+diff(rev(p(D)(x,lower=FALSE))))/2)
   }else{
      M <- q(D)(0.5);   L <- length(x)
      x.l <- x [ x <= M ];  x.u <- x [ x >= M ]
      L.l <- length(x.l);   L.u <- length(x.u)
      if(L.u+L.l-L>0)
         return(c(diff(p(D)(x.l, lower.tail = TRUE)),
               -diff(p(D)(x.u, lower.tail = FALSE))))
      else
         return(c(diff(p(D)(x.l, lower.tail = TRUE)),
                p(D)(x.u[1])-p(D)(x.l[L.l]),
               -diff(p(D)(x.u, lower.tail = FALSE))))
            }
   }



#------------------------------------------------------------------------------
# .makeD, .makeP, .makeQ
#------------------------------------------------------------------------------

.makeD <- function(object, argList, stand = NULL, fac = NULL)
         { d <- function(x, log = FALSE, ...){}
           if(inA <- .inArgs("log", object@d))
                     argList <- substitute(c(AL, log = log, ...),
                                           list(AL = argList))
           else      argList <- substitute(c(AL, ...),
                                           list(AL = argList))
           myCall <- substitute(d0 <- do.call("@"(objC,d),AL) ,
                                list(objC = quote(object),
                                AL = argList))
           myBod0 <- myCall
           if (!inA) myBod0 <- substitute({myCall
                                           if (log) d0 <- log(d0)})
           myBod1 <- myBod0
           if (!is.null(stand))
                myBod1 <- substitute(
                     {myBod0
                      d0 <- if (log) d0 - log(stand) else d0 / stand },
                      list(stand = stand, myBod0 = myBod0))
           myBod <- myBod1
           if (!is.null(fac))
                myBod <- substitute(
                     {myBod1
                      d0 <- if (log) d0 + log(fac) else d0 * fac },
                      list(fac = fac, myBod1 = myBod1))
           body(d) <- substitute({myBod
                                  return(d0)})
           return(d)
         }

.makeP <- function(object, argList, sign = TRUE, correct = NULL, 
          fac = NULL, fac2 = NULL)
         {
           p <- function(q, lower.tail = TRUE, log.p = FALSE, ...){}
           siY <- substitute(lower.tail)
           siN <- substitute(!lower.tail)
           lowS <-  if( sign) siY else siN
           lowS1 <- if(!sign) siY else siN
           if(inA1 <- .inArgs("lower.tail", object@p))
                      argList <- substitute(c(aL, lower.tail = lowT),
                                list( lowT = lowS,
                                      aL = argList )
                                      )
           else  argList <- substitute(aL, list(aL = argList ))
           if(inA2 <- (.inArgs("log.p", object@p) && inA1))
                     argList <- substitute(c(argList, log.p = log.p, ...))
           else      argList <- substitute(c(argList, ...))

           myCall <- substitute(p0 <- facC * do.call("@"(objC,p),AL) + facD,
                                list(objC = quote(object),
                                     AL = argList,
                                     facC = if(is.null(fac)) 1 else fac,
                                     facD = if(is.null(fac2)) 0 else fac2)
                                     )
           if (!inA1)
               myBod0 <- substitute({myC
                                     if (lowT) p0 <- 1 - p0},
                         list(myC=myCall, lowT = lowS1))
           else
               myBod0 <- substitute(myCall)

           if (!sign && !is.null(correct))
               myBod1 <- substitute({myC
                                     corC},
                         list(myC = myBod0, corC = correct))
           else myBod1 <- substitute(myBod0)

           if (!inA2) myBod <- substitute({myBod1
                                           if (log.p) p0 <- log(p0)
                                           return(p0)}  )

           else myBod <- substitute({myBod1
                                     return(p0)})
           body(p) <- substitute(myBod)
           return(p)
         }


.makeQ <- function(object, lastCall, sign = TRUE, Cont = TRUE)
     ## lastCall of form e1 %op% e2
         {
           siY <- substitute(lower.tail)
           siN <- substitute(!lower.tail)
           lowS <-  if( sign) siY else siN
           lowS1 <- if(!sign) siY else siN

           q <- function(p, lower.tail = TRUE, log.p = FALSE, ...){}
           argList <- alist(p)
           if(inA1 <- .inArgs("lower.tail", object@q))
                     argList1 <- substitute(alist(lower.tail = lowT),
                                list( lowT = lowS))
           else argList1 <- NULL

           if(inA2 <- (.inArgs("log.p", object@q) && inA1))
                      argList1 <- substitute(c(argList1,
                                               alist(log.p = log.p), ...))
           else       argList1 <- substitute(c(argList1, alist(...)))


           argList <- substitute(c(AL,AL1),list(AL=argList,AL1=argList1))
           myCall <- substitute({q0 <- do.call("@"(objC,q),AL)
                                 q0 <- lasC
                                 return(q0)},
                                list(objC = quote(object), AL=argList,
                                     lasC = lastCall))

           if (!Cont && !sign){
               if (inA1){
                     dqfS <- substitute({
                             dq <- getdistrOption("DistrResolution")
                             if (lower.tail) dq <- -dq
                                      })
               }else{
                     dqfS <- substitute({
                             dq <- getdistrOption("DistrResolution")
                                        })
               }
               if (inA2){ indS <- substitute(if(log.p){ ind <- (p>-Inf)&(p<0)
                                               p0[!ind] <- +1 }
                                             else { ind <- (p>0)&(p<1)
                                               p0[!ind] <- -1 } )
               }else    { indS <- substitute({ ind <- (p>0)&(p<1)
                                               p0[!ind] <- -1 })
               }
               myBod1 <- substitute({
                   q01 <- do.call("@"(object,q),AL)
                   p0 <- do.call("@"(object,p), c(alist(q=q01),AL1))
                   indC
                   ind2 <- .isEqual(p,p0)
                   dqfC # dq * if (with lower.tail-arg & lower.tail == FALSE)
                        # -1 else 1
                   p[ind2] <- p[ind2] + dq
                   my0C
               }, list(my0C = substitute(myCall),
                       AL = substitute(argList), AL1 = substitute(argList1),
                       indC = substitute(indS), dqfC = substitute(dqfS) )
                       )
                 ## discrete correction
           }else{
               myBod1 <- substitute(myCall)
           }

           if (!inA1) myBod2 <- substitute({if (lowT) p <- 1 - p
                                            myC},
                                    list(lowT = lowS1, myC = myBod1 ))
           else myBod2 <- substitute(myBod1)

           if (!inA2){
               myBod <- substitute({if (log.p) p <- exp(p)
                                    myC   },
                         list(myC = myBod2))
               }else{ myBod <- substitute(myBod2) }

           body(q) <- substitute(myBod)
           return(q)
         }

#------------------------------------------------------------------------------
# .plusm, .multm
#------------------------------------------------------------------------------

.plusm <- function(e1, e2, Dclass = "DiscreteDistribution"){
            if (length(e2)>1) stop("length of operator must be 1")
            if (.isEqual(e2, 0)) return(e1)

            if(Dclass %in% c("AffLinDiscreteDistribution",
                              "AffLinAbscontDistribution"))
               if(.isEqual(e1@a,1)&&.isEqual(e1@b+e2,0)) return(e1@X0)

            if ((Dclass == "DiscreteDistribution")||
                (Dclass == "AffLinDiscreteDistribution"))
                supportnew <- e1@support + e2

            if ((Dclass == "AbscontDistribution")||
                (Dclass == "AffLinAbscontDistribution"))
                 gapsnew <- if(is.null(e1@gaps)) NULL else e1@gaps + e2

            rnew <- function(n, ...){}
            body(rnew) <- substitute({ f(n, ...) + g },
                                         list(f = e1@r, g = e2))
            dnew <- .makeD(e1, substitute(alist(x = x - e2), list(e2 = e2)))
            pnew <- .makeP(e1, substitute(alist(q = q - e2), list(e2 = e2)))
            qnew <- .makeQ(e1, substitute(q0 + e2, list(e2 = e2)))

            if (Dclass == "AffLinDiscreteDistribution"){
                object <- new("AffLinDiscreteDistribution", 
                              r = rnew, d = dnew, p = pnew,
                              q = qnew, support = supportnew,
                              a = e1@a, b = e1@b + e2, X0 = e1@X0,
                             .withSim = FALSE, .withArith = TRUE,
                    .logExact = .logExact(e1), .lowerExact = .lowerExact(e1))
                rm(supportnew)

            }else if (Dclass == "DiscreteDistribution"){
                object <- new("AffLinDiscreteDistribution", 
                              r = rnew, d = dnew, p = pnew,
                              q = qnew, support = supportnew,
                              a = 1, b = e2, X0 = e1,
                             .withSim = FALSE, .withArith = TRUE,
                    .logExact = .logExact(e1), .lowerExact = .lowerExact(e1))
                rm(supportnew)

            }else if (Dclass == "AffLinAbscontDistribution"){
                object <- new("AffLinAbscontDistribution", 
                              r = rnew, d = dnew, p = pnew, q = qnew, 
                              gaps = gapsnew, a = e1@a, b = e1@b + e2, 
                              X0 = e1@X0, .withSim = FALSE, .withArith = TRUE,
                    .logExact = .logExact(e1), .lowerExact = .lowerExact(e1))

            }else if (Dclass == "AbscontDistribution"){  
                object <- new("AffLinAbscontDistribution", r = rnew, 
                              d = dnew, p = pnew, q = qnew, gaps = gapsnew, 
                              a = 1, b = e2, X0 = e1, .withSim = FALSE, 
                              .withArith = TRUE,
                    .logExact = .logExact(e1), .lowerExact = .lowerExact(e1))
            }
            rm(pnew, qnew, dnew, rnew)
            object
          }

.multm <- function(e1, e2, Dclass = "DiscreteDistribution"){

            withg <- getdistrOption("withgaps")
            on.exit(distroptions(withgaps=withg))
            distroptions(withgaps=FALSE)

            if (length(e2)>1) stop("length of operator must be 1")

            if (.isEqual(e2, 1)) return(e1)
            if (.isEqual(e2, 0))  return(new("Dirac", location = 0))

            if(Dclass %in% c("AffLinDiscreteDistribution",
                              "AffLinAbscontDistribution"))
               if(.isEqual(e1@a*e2,1)&&.isEqual(e1@b,0)) return(e1@X0)

            rnew <- function(n, ...){}
            body(rnew) <- substitute({ f(n, ...) * g },
                                         list(f = e1@r, g = e2))

            if (Dclass == "AffLinDiscreteDistribution"){
                 supportnew <- e1@support * e2
                 if (e2 < 0) supportnew <- rev(supportnew)

                 coR <- substitute({
                             o.warn <- getOption("warn"); options(warn = -1)
                             on.exit(options(warn=o.warn))
                             #
                             x0 <- .setEqual(q / e2C, support(object))
                             d0 <- object@d(x = x0)*(x0 %in% support(object))
                             #
                             options(warn = o.warn)
                             if (!lower.tail) d0 <- -d0
                             p0 <- p0 + d0},
                             list(e2C = e2)
                             )

                 dnew <- .makeD(substitute(e1, list(e1 = e1)),
                                substitute(alist(x = x / e2), list(e2 = e2)))
                 pnew <- .makeP(substitute(e1, list(e1 = e1)),
                                substitute(alist(q = q / e2), list(e2 = e2)),
                                sign = e2>0, correct = coR)
                 qnew <- .makeQ(substitute(e1, list(e1 = e1)),
                                substitute(q0 * e2, list(e2 = e2)),
                                sign = e2>0, Cont = FALSE)
                 object <- new(Dclass, r = rnew, d = dnew, p = pnew,
                               q = qnew, support = supportnew,
                               a = e1@a * e2, b = e2 * e1@b, X0 = e1@X0,                              
                              .withSim = FALSE, .withArith = TRUE,
                    .logExact = .logExact(e1), .lowerExact = .lowerExact(e1))
                 rm(supportnew)

            }else if (Dclass == "DiscreteDistribution"){
                 supportnew <- e1@support * e2
                 if (e2 < 0) supportnew <- rev(supportnew)

                 coR <- substitute({
                             o.warn <- getOption("warn"); options(warn = -1)
                             on.exit(options(warn=o.warn))
                             d0 <- object@d(x = q / e2C)
                             options(warn = o.warn)
                             if (!lower.tail) d0 <- -d0
                             p0 <- p0 + d0},
                             list(e2C = e2)
                             )

                 dnew <- .makeD(substitute(e1, list(e1 = e1)),
                                substitute(alist(x = x / e2), list(e2 = e2)))
                 pnew <- .makeP(substitute(e1, list(e1 = e1)),
                                substitute(alist(q = q / e2), list(e2 = e2)),
                                sign = e2>0, correct = coR)
                 qnew <- .makeQ(substitute(e1, list(e1 = e1)),
                                substitute(q0 * e2, list(e2 = e2)),
                                sign = e2>0, Cont = FALSE)
                 object <- new("AffLinDiscreteDistribution", r = rnew, d = dnew, 
                               p = pnew, q = qnew, support = supportnew,
                               a = e2, b = 0, X0 = e1,                              
                              .withSim = FALSE, .withArith = TRUE,
                    .logExact = .logExact(e1), .lowerExact = .lowerExact(e1))
                 rm(supportnew)

            }else if (Dclass == "AffLinAbscontDistribution"){
                 trY <- try(
                 if(is.null(e1@gaps))
                    gapsnew <- NULL
                 else {gapsnew <- e1@gaps
                       if(is.numeric(gapsnew)&&length(gapsnew)>0){
                          gapsnew <- matrix(gapsnew * e2, ncol=2)
                          if (e2 < 0) gapsnew <-
                             gapsnew[rev(seq(nrow(gapsnew))),c(2,1),drop = FALSE]
                       }
                 }, silent=TRUE)
                 if(is(trY,"try-error")) gapsnew <- NULL
                 dnew <- .makeD(substitute(e1, list(e1 = e1)),
                                substitute(alist(x = x / e2), list(e2 = e2)),
                                stand = abs(e2))
                 pnew <- .makeP(substitute(e1, list(e1 = e1)),
                                substitute(alist(q = q / e2), list(e2 = e2)),
                                sign = e2>0)
                 qnew <- .makeQ(substitute(e1, list(e1 = e1)),
                                substitute(q0 * e2, list(e2 = e2)),
                                sign = e2>0)
                 object <- new(Dclass, r = rnew, d = dnew, p = pnew, q = qnew, 
                               gaps = gapsnew, a = e1@a * e2, b = e2 * e1@b, 
                               X0 = e1@X0, .withSim = FALSE, .withArith = TRUE,
                    .logExact = .logExact(e1), .lowerExact = .lowerExact(e1))
 
            }else if (Dclass == "AbscontDistribution"){
                 trY <- try(
                 if(is.null(e1@gaps))
                    gapsnew <- NULL
                 else {gapsnew <- e1@gaps
                       if(is.numeric(gapsnew)&&length(gapsnew)>0){
                             gapsnew <- matrix(gapsnew * e2, ncol=2)
                             if (e2 < 0) gapsnew <-
                                    gapsnew[rev(seq(nrow(gapsnew))),
                                            c(2,1),drop = FALSE]
                             }
                 }, silent=TRUE)
                 if(is(trY,"try-error")) gapsnew <- NULL

                 dnew <- .makeD(substitute(e1, list(e1 = e1)),
                                substitute(alist(x = x / e2), list(e2 = e2)),
                                stand = abs(e2))
                 pnew <- .makeP(substitute(e1, list(e1 = e1)),
                                substitute(alist(q = q / e2), list(e2 = e2)),
                                sign = e2>0)
                 qnew <- .makeQ(substitute(e1, list(e1 = e1)),
                                substitute(q0 * e2, list(e2 = e2)),
                                sign = e2>0)
                 object <- new("AffLinAbscontDistribution", r = rnew, d = dnew, 
                               p = pnew, q = qnew, gaps = gapsnew, a = e2, 
                               b = 0, X0 = e1, .withSim = FALSE, 
                               .withArith = TRUE,
                    .logExact = .logExact(e1), .lowerExact = .lowerExact(e1))
            }
            rm(pnew, qnew, dnew, rnew)
            object
          }
          

#------------------------------------------------------------------------------
# .makeDd, .makePd, .makeQd
#------------------------------------------------------------------------------



.makeDd <- function(x,y, yleft, yright){
   intervall <- getdistrOption("DistrResolution") / 2
   supp.grid <- c(matrix(rbind(x - intervall, x + intervall), nrow = 1))
   prob.grid <- c(matrix(rbind(0, y), nrow = 1), 0)
   df0 <- stepfun(x = supp.grid, y = prob.grid)
   rm(intervall, supp.grid, prob.grid)
   return(df0)
}

.makePd <- function(x,y, yleft, yright){
   stepfun(x = x, y = c(yleft, y))
}

.makeQd <- function(x,y, yleft, yright){
force(y)
force(x)
f <- function(u) {
               q0 <- sapply(u, 
                       function(z) y[min(sum(x < z-.Machine$double.eps) + 1,
                                         length(y)) ] )
               q0[.isEqual(u,0)] <- yleft
               q0[.isEqual(u,1)] <- yright
               return(q0)}
### individualize local variables
#environment(f) <- new.env()
#assign("y", y, environment(f))
#assign("x", x, environment(f))
#assign("yright", yright, environment(f))
#assign("yleft", yleft, environment(f))
return(f)
}

.makeQc <- function(x,y, yleft, yright){
#f0 <- function(u) {
#               q0 <- sapply(u, 
#                       function(z) y[min(sum(x < z-.Machine$double.eps) + 1,
#                                         length(y)) ] )
#               return(q0)}
#eps <- getdistrOption("TruncQuantile")
#x00 <- seq(eps, 1-eps, length = getdistrOption("DefaultNrGridPoints"))
#y00 <- f0(x00)
#idx <- cumsum(rle(y00)[[1]]) ### use only unique y's and corresponding 
#x0 <- x00[idx]               ### maximal x's
#y0 <- y00[idx]
#f1 <- approxfun(x = x0, y = y0, yleft = y0[1], yright = y0[length(y0)])
yleft <- yleft[1]
yright <- yright[1]

isna <- is.na(x)|is.na(y)
x <- x[!isna]
y <- y[!isna]

l0 <- length(unique(x[!.isEqual01(x)]))
if(l0 > 1){
   yl <- if(!is.na(yleft) && is.finite(yleft))  yleft  else y[1]
   yr <- if(!is.na(yright)&& is.finite(yright)) yright else y[length(y)]

   f1 <- approxfun(x = x, y = y, yleft = yl, yright = yr)
}else{ 
   i0 <- (1:length(x))[x==unique(x[!.isEqual01(x)])]
   y0 <- if(l0 ==1) y[min(i0)] else yleft
   f1 <- function(x) return(y0)
}
f <- function(x) 
   {y1 <- f1(x)
    y1[.isEqual(x,0)] <- yleft
    y1[.isEqual(x,1)] <- yright
    return(y1)
    }
return(f)
}


#------------------------------------------------------------------------------
# .makeDNew, .makePNew, .makeQNew
#------------------------------------------------------------------------------


.makeDNew <- function(x, dx, h = NULL, Cont = TRUE, standM = "sum"){
            dx <- (dx >= .Machine$double.eps)*dx
            if( length(dx) < length(x) ) dx <- c(0,dx)

            if (is.null(h)) h <- 1

            dx1 <- dx / h
            mfun <- if (Cont) approxfun else .makeDd

            ## density
            df1 <- mfun(x = x, y = dx1, yleft = 0, yright = 0)

            if (standM == "sum")
                   stand <- sum(dx)
            else   {
            stand <- try(integrate(df1, -Inf, Inf)$value, TRUE)
            if (is(stand,"try-error")){
               if(getdistrOption("warn.makeDNew"))
                  warning("'integrate()' threw an error ---result may be inaccurate.")
               stand <- sum(df1(x))*h*(x[2]-x[1])
               }
            }
            dfun <- function(x, log = FALSE)
                    {if (log)
                          d0 <-    log(df1(x))-log(stand)
                     else d0 <- df1(x) / stand
                     return (d0)}
            rm(x,dx1,h)
            return(dfun)
}

.primefun <- function(f,x, nm = NULL){
 
 h <- diff(x)
 l <- length(x)

 xm <- (x[-l]+x[-1])/2

 fxm <- f(xm)
 fx <- f(x)
 
 
 fxs  <- 2 * cumsum(fx) - fx - fx[1]
 fxsm <- 4 * cumsum(fxm)

 fxx <- c(0, (fxs[-1]+fxsm)* h / 6 )

 if (is.null(nm)) nm <- fxx[l]

 fx1 <- approxfun(x, fxx, yright = nm, yleft = 0)

 ffx <- function(u){
      ffy <- fx1(u) 
      ffy[u > max(x)] <- nm 
      ffy[u < min(x)] <- 0
      return(ffy)
     }

 return(ffx)
}

.csimpsum <- function(fx){
 l <- length(fx)
 l2 <- l%/%2
 if (l%%2 == 0) {
     fx <- c(fx[1:l2],(fx[l2]+fx[l2+1])/2,fx[(l2+1):l])
     l <- l+1}
 f.even <- fx[seq(l) %% 2 == 0]
 f.odd  <- fx[seq(l) %% 2 == 1]
 fs    <- 2 * cumsum(f.odd) - f.odd - f.odd[1]
 fsm   <- 4 * cumsum(f.even)
 ff <- c(0,(fs[2:(l2+1)]+fsm)/3 )
 ff
}

.makePNew <- function(x, dx, h = NULL, notwithLLarg = FALSE,
                      Cont = TRUE, myPf = NULL, pxl = NULL, pxu = NULL){

  if (is.null (h)) h <- 0

  x.u <- x.l <- x
  if (Cont){
         mfun <- if (is.null (myPf)) approxfun else myPf
         l <- length(x)
         if ((l%%2==0)&& is.null(myPf)){
               l2 <- l/2
               if (is.null(pxl))
                   x.l <- c(x[1:l2],(x[l2]+x[l2+1])/2,x[(l2+1):l])
               if (is.null(pxu))
                   x.u <- c(x[1:l2],(x[l2]+x[l2+1])/2,x[(l2+1):l])
               l <- l+1
               }
         cfun <- .csimpsum
         if (is.null(pxl)&& is.null(myPf))
             x.l <- x.l[seq(l)%%2==1]
         if (is.null(pxu)&& is.null(myPf))
             x.u <- x.u[seq(l)%%2==1]
  }else    {
         mfun <- .makePd
         cfun <- cumsum
  }       

  p.l <- if(!is.null(pxl)) pxl else cfun(dx)
  
  nm <- max(p.l)
  p1.l <- mfun(x = x.l, y = p.l, yleft = 0, yright = nm)
  nm <- p1.l(max(x))

  if(notwithLLarg){
      ifElsePS <- substitute(if (lower.tail) p1.l(q) else 1 - p1.l(q))
  }else{
      p.u <- if(!is.null(pxu)) pxu else rev(cfun(rev(dx)))
    ## continuity correction by h/2
      if (!Cont) p.u <- c(p.u[-1],0)
      p1.u <- mfun(x = x.u, y = p.u, yright = 0, yleft = nm)
      rm(p.u)
      ifElsePS <- substitute(if (lower.tail) p1.l(q) else p1.u(q))
  }
  pfun <- function(q, lower.tail = TRUE, log.p = FALSE){}
  body(pfun) <- substitute(
              { p0 <- ifElsePC
                p0 <- if (log.p) log(p0)-log(nm) else p0/nm
                return(p0)
              }, list(ifElsePC = ifElsePS))
  rm(dx, p.l, notwithLLarg)
  return(pfun)
}


.makeQNew <- function(x, px.l, px.u, notwithLLarg = FALSE, yL , yR,
                      Cont = TRUE){
  o.warn <- getOption("warn"); options(warn = -1)
  on.exit(options(warn=o.warn))
  ix <- .isEqual01(px.l)
  mfun <- if (Cont) .makeQc else  .makeQd
  if(!is.finite(yR)||Cont)
     {xx <- px.l[!ix]; yy <- x[!ix]}
  else  
     {xx <- px.l; yy <- x}
  q.l <- mfun(x = xx, y = yy, yleft = yL, yright = yR)
  rm(xx,yy)
  if(notwithLLarg){
     ifElseQS <- quote(if (lower.tail) q.l(p01) else q.l(1-p01))
  }else{
#         px.u <- rev(px.u);
         x <- rev(x)
         if (Cont) px.u <- rev(px.u)
         ix <- .isEqual01(px.u)
         xx <- px.u[!ix]
         yy <- if (Cont) x[!ix] else x[rev(!ix)]
         q.u <- mfun(x = xx, y = yy, yleft = yR, yright = yL)
         rm(xx,yy)
     ifElseQS <- quote(if (lower.tail) q.l(p01) else q.u(p01))
  }
  options(warn = o.warn)
  qfun <- function(p, lower.tail = TRUE, log.p = FALSE){}
  body(qfun) <- substitute({
          if (log.p) p <- exp(p)
          if (any((p < -.Machine$double.eps)|(p > 1+.Machine$double.eps)))
              warning(gettextf("q method of %s produced NaN's ", objN))
              i01 <- (-.Machine$double.eps<=p)&(p<=1+.Machine$double.eps)
              p01 <- p[i01] ## only values in [0,1] are used
              q0  <- p*0
              q0[!i01] <- NaN
              q0[ i01] <- ifElseQC
              return(as.numeric(q0))
              }, list(ifElseQC = ifElseQS, objN = quote(.getObjName(1))))
  return(qfun)
}




#setMethod("+", c("LatticeDistribution","AbscontDistribution"),plusDAC)

#------------------------------------------------------------------------------
# .expm.[dc], .logm.[dc] (exact calculation of exp /log) 
#    [ similarly possible : other monotone, smooth trafos, e.g. tan / atan ]
#------------------------------------------------------------------------------
.expm.d <- function(e1){
            supportnew <- exp(e1@support)

            rnew <- function(n, ...){}
            body(rnew) <- substitute({ exp(f(n, ...)) },
                                         list(f = e1@r))
            dnew <- .makeD(substitute(e1, list(e1 = e1)),
                           substitute(alist(x = log(x*(x>0)+(x<=0)))),
                           fac = substitute((x>0)))
            pnew <- .makeP(substitute(e1, list(e1 = e1)),
                           substitute(alist(q = log(q*(q>0)+(q<=0)))),
                           fac = substitute((q>0)),
                           fac2 = substitute((q<=0)& !lower.tail))
            qnew <- .makeQ(substitute(e1, list(e1 = e1)),
                           substitute(exp(q0)),
                           Cont = FALSE)

            object <- new("DiscreteDistribution",
                           r = rnew, d = dnew, p = pnew,
                           q = qnew, support = supportnew,
                          .withSim = FALSE, .withArith = TRUE,
                    .logExact = .logExact(e1), .lowerExact = .lowerExact(e1))
            rm(supportnew)
            rm(pnew, qnew, dnew, rnew)
            object
          }
.expm.c <- function(e1){
            gapsnew <- if(is.null(e1@gaps)) NULL else exp(e1@gaps)

            rnew <- function(n, ...){}
            body(rnew) <- substitute({ exp(f(n, ...)) },
                                         list(f = e1@r))
                                    
            
            dnew0 <- .makeD(substitute(e1, list(e1 = e1)),
                            substitute(alist(x = log(x*(x>0)+(x<=0)))),
                            stand = substitute(x+(x==0)),
                            fac = substitute((x>0)))

            # extrapolation for x=0
            x0a <- 10^(-(1:10)/4) 
            f0a <- dnew0(x = x0a, log = FALSE)
            f0  <- max(spline(x0a,f0a, xout = 0)$y,0)
            lf0 <- log(f0)

            dnew <- function(x, log = FALSE, ...){
               d0 <- dnew0(x, log = log, ...)
               i0 <- .isEqual(x,0)
               if(!any(i0)) return(d0)
               if(log) d0[i0] <- lf0 else d0[i0] <- f0
               return(d0)
            }
            pnew <- .makeP(substitute(e1, list(e1 = e1)),
                           substitute(alist(q = log(q*(q>0)+(q<=0)))),
                           fac = substitute((q>0)),
                           fac2 = substitute((q<=0)& !lower.tail))
            qnew <- .makeQ(substitute(e1, list(e1 = e1)),
                           substitute(exp(q0)))

            object <- AbscontDistribution(
                           r = rnew, d = dnew, p = pnew,
                           q = qnew, gaps = gapsnew,
                          .withSim = FALSE, .withArith = TRUE,
                    .logExact = .logExact(e1), .lowerExact = .lowerExact(e1))
            if(exists("gapsnew")) rm(gapsnew)
            rm(pnew, qnew, dnew, rnew)
            object
          }
.logm.d <- function(e1){
            supportnew <- log(e1@support)

            rnew <- function(n, ...){}
            body(rnew) <- substitute({ log(f(n, ...)) },
                                         list(f = e1@r))
            dnew <- .makeD(substitute(e1, list(e1 = e1)),
                           substitute(alist(x = exp(x))))
            pnew <- .makeP(substitute(e1, list(e1 = e1)),
                           substitute(alist(q = exp(q))))
            qnew <- .makeQ(substitute(e1, list(e1 = e1)),
                           substitute(log(q0)),
                           Cont = FALSE)

            object <- new("DiscreteDistribution",
                           r = rnew, d = dnew, p = pnew,
                           q = qnew, support = supportnew,
                          .withSim = FALSE, .withArith = TRUE,
                    .logExact = .logExact(e1), .lowerExact = .lowerExact(e1))
            rm(supportnew)
            rm(pnew, qnew, dnew, rnew)
            object
          }
.logm.c <- function(e1){
            gapsnew <- if(is.null(e1@gaps)) NULL else log(e1@gaps)

            rnew <- function(n, ...){}
            body(rnew) <- substitute({ log(f(n, ...)) },
                                         list(f = e1@r))
            dnew <- .makeD(substitute(e1, list(e1 = e1)),
                           substitute(alist(x = exp(x))),
                           fac = substitute(exp(x)))
            pnew <- .makeP(substitute(e1, list(e1 = e1)),
                           substitute(alist(q = exp(q))))
            qnew <- .makeQ(substitute(e1, list(e1 = e1)),
                           substitute(log(q0)))

            object <- AbscontDistribution(
                           r = rnew, d = dnew, p = pnew,
                           q = qnew, gaps = gapsnew,
                          .withSim = FALSE, .withArith = TRUE,
                    .logExact = .logExact(e1), .lowerExact = .lowerExact(e1))
            if(exists("gapsnew")) rm(gapsnew)
            rm(pnew, qnew, dnew, rnew)
            object
          }


#------------------------------------------------------------------------------
# .P2D, .D2P, .P2Q, .Q2P casting functions
#------------------------------------------------------------------------------

###Functions for AbscontDistribution 

#determines slot d from p
.P2D <- function(p, xx, ql, qu, ngrid = getdistrOption("DefaultNrGridPoints"))
{if(missing(xx))
    xx <- seq(ql, qu, length = ngrid)
 px <- p(xx)
 dx <- D1ss(xx,px)
 return(.makeDNew(xx, dx, h = NULL, Cont = TRUE, standM = "integrate"))
}


#determines slot q from p

.P2Q <- function(p, xx, ql,qu, ngrid = getdistrOption("DefaultNrGridPoints"), 
                qL = -Inf, qU = Inf){
if(missing(xx))
   {xx <- seq(ql, qu, length = ngrid)
     h <- (qu-ql)/ngrid}
else h <- c(diff(xx),(rev(xx)[1]-rev(xx)[2]))
px.l <- p(xx + 0.5*h)
px.u <- p(xx + 0.5*h, lower.tail = FALSE)
return(.makeQNew(xx + 0.5*h, px.l, px.u, FALSE, qL, qU))
}

#determines slot p from d
.D2P <- function(d, xx, ql, qu,  ngrid = getdistrOption("DefaultNrGridPoints"))
{if(missing(xx))
   { if(ngrid%%2==0) ngrid  <- ngrid+1
     xx <- seq(ql, qu, length = ngrid)
     h <- (qu-ql)/ngrid}
 else h <- c(diff(xx),(rev(xx)[1]-rev(xx)[2]))

 dx <- d(xx)
 return(.makePNew(xx, dx, h = h))
}

#determines slot p from q

.Q2P <- function(q, ngrid = getdistrOption("DefaultNrGridPoints")){
ep <- getdistrOption("TruncQuantile")^2
xx0 <- seq(ep, 1-ep, length = ngrid)
qx0 <- q(xx0)
qx1 <- unique(qx0)
x1 <- tapply(xx0, qx0, max)
p0 <- approxfun(x = qx1, y = x1, yleft = 0, yright = 1, rule = 2)
p01 <- function(x) {
       p11 <- p0(x) 
       p11[x>max(qx1)] <- 1
       p11[x<min(qx1)] <- 0
       return(p11)}
return(function(q, lower.tail = TRUE, log.p = FALSE){
  p1 <- p01(q)
  p1[q==Inf] <- 1
  p1[q==-Inf] <- 0
  if(!lower.tail) p1 <- 1-p1
  if(log.p) p1 <- log(p1)
  return(p1)
})
}

#------------------------------------------------------------------------------
# modify slot q for AbscontDistribution if there are gaps
#------------------------------------------------------------------------------
.modifyqgaps <- function(pfun, qfun, gaps, leftright = "left"){
  if(length(gaps)==0) return(qfun)
  p.gaps <- pfun(gaps[,1]) 
  p.gaps.l <- pfun(gaps[,1], lower.tail = FALSE)
  dP <- deparse(body(qfun))
  syC <- paste(sys.calls())
  if(length(grep("q\\.r\\(",syC[length(syC)-2])) == 0)
     if(length(grep("q0 <- qfun", dP[2]))>0) 
        return(qfun)
  if(pmatch(leftright, table = c("left","right"), nomatch = 1) == 1){
     qnew <- function(p, lower.tail = TRUE, log.p = FALSE) {
             q0 <- qfun(p, lower.tail = lower.tail, log.p = log.p)
             i0 <- seq(length=length(p))
             if(lower.tail){
                if(log.p) p.gaps <- log(p.gaps)
                for(i in 1:nrow(gaps)){
                    i1 <- .isEqual(p,p.gaps[i])
                    i2 <- i0[p<p.gaps[i]]
                    i3 <- i0[p>p.gaps[i]]
                    q0[i1] <- gaps[i,1]
                    if(length(i2)>0)
                       q0[i2] <- pmin(q0[i2],gaps[i,1])
                    if(length(i3)>0)
                       q0[i3] <- pmax(q0[i3],gaps[i,2])
                }        
             }else{
                if(log.p) p.gaps.l <- log(p.gaps.l)
                for(i in 1:nrow(gaps)){
                    i1 <- .isEqual(p,p.gaps.l[i])
                    i2 <- i0[p<p.gaps.l[i]]
                    i3 <- i0[p>p.gaps.l[i]]
                    q0[i1] <- gaps[i,2]
                    if(length(i2)>0)
                       q0[i2] <- pmax(q0[i2],gaps[i,2])    
                    if(length(i3)>0)
                       q0[i3] <- pmin(q0[i3],gaps[i,1])    
                }        
             }
             return(q0)
     }
  }else{
     qnew <- function(p, lower.tail = TRUE, log.p = FALSE) {
             q0 <- qfun(p, lower.tail = lower.tail, log.p = log.p)
             i0 <- seq(length=length(p))
             if(lower.tail){
                if(log.p) p.gaps <- log(p.gaps)
                for(i in 1:nrow(gaps)){
                    i1 <- .isEqual(p,p.gaps[i])
                    i2 <- i0[p<p.gaps[i]]
                    i3 <- i0[p>p.gaps[i]]
                    q0[i1] <- gaps[i,2]
                    if(length(i2)>0)
                       q0[i2] <- pmin(q0[i2],gaps[i,1])
                    if(length(i3)>0)
                       q0[i3] <- pmax(q0[i3],gaps[i,2])
                }        
             }else{
                if(log.p) p.gaps.l <- log(p.gaps.l)
                for(i in 1:nrow(gaps)){
                    i1 <- .isEqual(p,p.gaps.l[i])
                    i2 <- i0[p<p.gaps.l[i]]
                    i3 <- i0[p>p.gaps.l[i]]
                    q0[i1] <- gaps[i,1]
                    if(length(i2)>0)
                       q0[i2] <- pmax(q0[i2],gaps[i,2])    
                    if(length(i3)>0)
                       q0[i3] <- pmin(q0[i3],gaps[i,1])    
                }        
             }
             return(q0)
     }
  }
  return(qnew)           
}

#------------------------------------------------------------------------------
# issue warnings in show / print as to Arith or print
#------------------------------------------------------------------------------
.IssueWarn <- function(Arith,Sim){
    msgA1 <- msgA2 <- msgS1 <- msgS2 <- NULL
    if(Arith && getdistrOption("WarningArith")){ 
      msgA1 <- gettext(
       "arithmetics on distributions are understood as operations on r.v.'s\n")
      msgA2 <- gettext(
       "see 'distrARITH()'; for switching off this warning see '?distroptions'")
       }
    if(Sim && getdistrOption("WarningSim")){ 
      msgS1 <- gettext(
       "slots d,p,q have been filled using simulations; ")
      msgS2 <- gettext(
       "for switching off this warning see '?distroptions'")
       }
    return(list(msgA=c(msgA1,msgA2), msgS = c(msgS1,msgS2)))  
    }

#------------------------------------------------------------------------------
# fill a list acc. recycling rules
#------------------------------------------------------------------------------
.List <- function(list0) if(is.list(list0)) list0 else list(list0)

.fillList <- function(list0, len = length(list0)){
            if(is.null(list0)) return(vector("list",len))
            list0 <- .List(list0)
            if(len == length(list0)) 
               return(list0)
            i <- 0
            ll0 <- length(list0)
            li0 <- vector("list",len)
            if(ll0){
              while(i < len){
                 j <- 1 + ( i %% ll0)
                 i <- i + 1
                 li0[[i]] <- list0[[j]]
              }
            }
            return(li0)
}

#------------------------------------------------------------------------------
# .DistributionAggregate by Jacob van Etten, jacobvanetten@yahoo.com
# on a mail on Feb 27th 2009
#------------------------------------------------------------------------------
.DistrCollapse <- function(support, prob, 
                              eps = getdistrOption("DistrResolution")){
    supp <- support
    prob <- as.vector(prob)
    suppIncr <- diff(c(supp[1]-2*eps,supp)) < eps
    groups <- cumsum(!suppIncr)
    prob <- as.vector(tapply(prob, groups, sum))
    supp <- as.vector(tapply(supp, groups, quantile, probs = 0.5, type = 1)) 
           ### in order to get a "support member" take the leftmost median
    return(list(supp = supp, prob = prob))
#    newDistribution <- DiscreteDistribution(supp=supp,prob=prob)
#    return(newDistribution)
}


.panel.mingle <- function(dots, element){
  pF <- dots[[element]]
  if(is.list(pF)) return(pF)
  pFr <- if(typeof(pF)=="symbol") eval(pF) else{
     pFc <- as.call(pF)
     if(as.list(pFc)[[1]] == "list"){
        lis <- vector("list",length(as.list(pFc))-1)
        for(i in 1:length(lis)){
            lis[[i]] <- pFc[[i+1]]
        }
        lis
     }else pF
  }
  return(pFr)
}
