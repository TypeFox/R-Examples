## extra methods
## binary operators


if(!isGeneric("simplifyr")) 
    setGeneric("simplifyr", 
                function(object, 
                         size = 10^getdistrOption("RtoDPQ.e")) 
                         standardGeneric("simplifyr")
              )

setMethod("simplifyr", "UnivariateDistribution", 
          function(object, size = 10^getdistrOption("RtoDPQ.e")){
            Sample <- r(object)(size)       
            rneu <- function(n) sample(x = Sample, size = n, replace = TRUE)
            eval.parent(substitute(object@r<-rneu))           
          })


## function to automatically generate, starting from simulations, density, 
## quantile function and cdf
## first version for absolutely continuous, second for discrete distributions

## we use 10^RtoDPQExponent random numbers to generate new distr
## density should use DefaultNrGridPoints equally spaced points for evaluation

RtoDPQ <- function(r, e = getdistrOption("RtoDPQ.e"),
                      n = getdistrOption("DefaultNrGridPoints"), y = NULL){
  zz <- if(!is.null(y)) y else r(10^e)
  zz <- zz[!is.na(zz)]
  
  dxy <-  xy.coords(density(zz, n = n))
  dfun <- .makeDNew(dxy$x, dxy$y, standM = "int")

  pf0 <- function(x, y, yleft, yright) ecdf(x)
  pfun <- .makePNew(x=zz, dx=0, notwithLLarg=TRUE, myPf = pf0)
            ## quantile function

  yL <-  min(zz);   yR <-  max(zz); rm(zz)
  px.l <- pfun(dxy$x);   px.u <- pfun(dxy$x, lower.tail = FALSE)
  qfun <- .makeQNew(dxy$x, px.l, px.u, TRUE, yL, yR)

  rm(px.l, px.u, dxy, pf0)
  list(dfun = dfun, pfun = pfun, qfun = qfun)}


RtoDPQ.d <- function(r, e = getdistrOption("RtoDPQ.e")){
  zz <- r(10^e)
  X <- table(zz)
  rm(zz)

  supp <- as.numeric(names(X))
  prob <- X/(10^e)
  rm(X)

  len = length(supp)

  if(len > 1){
    if(min(diff(supp)) <
           getdistrOption("DistrResolution"))
       stop("grid too narrow --> change DistrResolution")
  }

  dfun <- .makeDNew(supp, prob, Cont = FALSE)
  pfun <- .makePNew(supp, prob, TRUE, Cont = FALSE)
  qfun <- .makeQNew(supp, cumsum(prob), rev(cumsum(rev(prob))),
                      TRUE, min(supp), max(supp), Cont = FALSE)

  list(dfun = dfun, pfun = pfun, qfun = qfun)
}

### new from 2.0:

RtoDPQ.LC <- function(r, e = getdistrOption("RtoDPQ.e"),
                      n = getdistrOption("DefaultNrGridPoints"), y = NULL){

  zz <- if(!is.null(y)) y else r(10^e)
  hasDis <- FALSE
  zz.nr <- zz

  zz.T <- table(zz)
  zz.T1 <- zz.T[zz.T>1]
  zz.replic <- as.numeric(names(zz.T1))
  w.d <- sum(zz %in% zz.replic)/10^e
  rm(zz.T)

  f.d <- Dirac(0)
  if(w.d)
  {hasDis <- TRUE
   zz.nr <- zz[! zz %in% zz.replic]
   d.r <- zz.T1/sum(zz.T1)
   f.d <- DiscreteDistribution(supp = zz.replic, prob = d.r,
                     .withSim = TRUE, .withArith = TRUE,
                     .lowerExact = FALSE, .logExact = FALSE)
   rm(d.r,zz.replic,zz.T1)
  }
  rm(zz)
  
  if(1-w.d){
  dxy <-  xy.coords(density(zz.nr, n = n))
  dcfun <- .makeDNew(dxy$x, dxy$y, standM = "int")

  pf0 <- function(x, y, yleft, yright) ecdf(x)
  pcfun <- .makePNew(x=zz.nr, dx=0, notwithLLarg=TRUE, myPf = pf0)
            ## quantile function

  yL <-  min(zz.nr);   yR <-  max(zz.nr); rm(zz.nr)
  px.l <- pcfun(dxy$x);   px.u <- pcfun(dxy$x, lower.tail = FALSE)
  qcfun <- .makeQNew(dxy$x, px.l, px.u, TRUE, yL, yR)

  rm(px.l, px.u, dxy, pf0)
  f.c <- AbscontDistribution( r= function(n) qcfun(runif(n)),
             d=dcfun, p = pcfun, q = qcfun, .withSim = TRUE,
             .withArith = TRUE, .lowerExact = FALSE, .logExact = FALSE)
  }
  else f.c <-Norm()
  UnivarLebDecDistribution(discretePart = f.d, acPart = f.c,
                           discreteWeight = w.d)
  }

####################################################################################


###########################################################


