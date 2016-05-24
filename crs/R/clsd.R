## These functions are for (currently univariate) logspline density
## estimation written by racinej@mcmaster.ca (Jeffrey S. Racine). They
## make use of spline routines in the crs package (available on
## CRAN). The approach involves joint selection of the degree and
## knots in contrast to the typical approach (e.g. Kooperberg and
## Stone) that sets the degree to 3 and optimizes knots only. Though
## more computationally demanding, the estimators are more efficient
## on average.

par.init <- function(degree,segments,monotone,monotone.lb) {

  ## This function initializes parameters for search along with upper
  ## and lower bounds if appropriate.

  dim.p <- degree+segments

  ## The weights for the linear tails must be non-positive. The lower
  ## bound places a maximum bound on how quickly the tails are allowed
  ## to die off. Trial and error suggests the values below seem to be
  ## appropriate for a wide range of (univariate)
  ## distributions. Kooperberg suggests that in order to get the
  ## constraint theta < 0 use theta <= -epsilon for some small epsilon
  ## > 0. We therefore use sqrt machine epsilon.

  par.ub <- - sqrt(.Machine$double.eps)
  par.lb <- monotone.lb
  par.init <- c(runif(1,-10,par.ub),rnorm(dim.p-2),runif(1,-10,par.ub))
  par.upper <- c(par.ub,rep(Inf,dim.p-2),par.ub)
  par.lower <- if(monotone){rep(-Inf,dim.p)}else{c(par.lb,rep(-Inf,dim.p-2),par.lb)}

  return(list(par.init=par.init,
              par.upper=par.upper,
              par.lower=par.lower))

}

gen.xnorm <- function(x=NULL,
                      xeval=NULL,
                      lbound=NULL,
                      ubound=NULL,
                      er=NULL,
                      n.integrate=NULL) {

  er <- extendrange(x,f=er)
  if(!is.null(lbound)) er[1] <- lbound
  if(!is.null(ubound)) er[2] <- ubound
  if(min(x) < er[1] | max(x) > er[2]) warning(" data extends beyond the range of `er'")
  xint <- sort(as.numeric(c(seq(er[1],er[2],length=round(n.integrate/2)),
                 quantile(x,seq(sqrt(.Machine$double.eps),1-sqrt(.Machine$double.eps),length=round(n.integrate/2))))))
  if(is.null(xeval)) {
    xnorm <- c(x,xint)
  } else {
    xnorm <- c(xeval,x,xint)
    if(min(xeval) < er[1] | max(xeval) > er[2]) warning(" evaluation data extends beyond the range of `er'")
  }
  ## Either x will be the first 1:length(x) elements in
  ## object[rank.xnorm] or xeval will be the first 1:length(xeval)
  ## elements in object[rank.xnorm]
  return(list(xnorm=xnorm[order(xnorm)],rank.xnorm=rank(xnorm)))

}

density.basis <- function(x=NULL,
                          xeval=NULL,
                          xnorm=xnorm,
                          degree=NULL,
                          segments=NULL,
                          basis="tensor",
                          knots="quantiles",
                          monotone=TRUE) {

  ## To obtain the constant of integration for B-spline bases, we need
  ## to compute log(integral exp(P%*%beta)) so we take an equally
  ## spaced extended range grid of length n plus the sample
  ## realizations (min and max of sample therefore present for what
  ## follows), and evaluation points xeval if they exist.

  ## Charles Kooperberg has a manuscript "Statistical Modeling with
  ## Spline Functions', Jan 5 2006 on his web page. Chapter 6, page
  ## 286, figure 6.7 reveals a hybrid spline basis that in essence
  ## appears to drop two columns from my B-spline with 2 segments
  ## added artificially. This has the effect of removing the two bases
  ## that were delivering weight in tails leading to `kinks'. Hat-tip
  ## to Charles for his clear descriptions. Note we require the same
  ## for the derivatives below. Note that this logspline basis does
  ## not have the B-spline property that the pointwise sum of the
  ## bases is 1 everywhere.

  suppressWarnings(Pnorm <- prod.spline(x=x,
                                        xeval=xnorm,
                                        K=cbind(degree,segments+if(monotone){2}else{0}),
                                        knots=knots,
                                        basis=basis))

  if(monotone) Pnorm <- Pnorm[,-c(2,degree+segments+1)]

  ## Compute the normalizing constant so that the estimate integrates
  ## to one. We append linear splines to the B-spline basis to
  ## generate exponentially declining tails (K=cbind(1,1) creates the
  ## linear basis).

  suppressWarnings(P.lin <- prod.spline(x=x,
                                        xeval=xnorm,
                                        K=cbind(1,1),
                                        knots=knots,
                                        basis=basis))

  ## We append the linear basis to the left and rightmost polynomial
  ## bases. We match the slope of the linear basis to that of the
  ## polynomial basis at xmin/xmax (note that
  ## Pnorm[xnorm==max(x),-ncol(Pnorm)] <- 0 is there because the
  ## gsl.bspline values at the right endpoint are very small but not
  ## exactly zero but want to rule out any potential issues hence set
  ## them correctly to zero)

  Pnorm[xnorm<min(x),] <- 0
  Pnorm[xnorm>max(x),] <- 0
  Pnorm[xnorm==max(x),-ncol(Pnorm)] <- 0
  P.left <- as.matrix(P.lin[,1])
  P.right <- as.matrix(P.lin[,2])

  ## We want the linear segment to have the same slope as the
  ## polynomial segment it connects with and to match at the joint
  ## hence conduct some carpentry at the left boundary.

  index <- which(xnorm==min(x))
  index.l <- index+1
  index.u <- index+5
  x.l <- xnorm[index.l]
  x.u <- xnorm[index.u]
  slope.poly.left <- as.numeric((Pnorm[index.u,1]-Pnorm[index.l,1])/(x.u-x.l))
  index.l <- index+1
  index.u <- index+5
  x.l <- xnorm[index.l]
  x.u <- xnorm[index.u]
  slope.linear.left <- as.numeric((P.left[index.u]-P.left[index.l])/(x.u-x.l))

  ## Complete carpentry at the right boundary.

  index <- which(xnorm==max(x))
  index.l <- index-1
  index.u <- index-5
  x.l <- xnorm[index.l]
  x.u <- xnorm[index.u]
  slope.poly.right <- as.numeric((Pnorm[index.u,ncol(Pnorm)]-Pnorm[index.l,ncol(Pnorm)])/(x.u-x.l))
  index.l <- index-1
  index.u <- index-5
  x.l <- xnorm[index.l]
  x.u <- xnorm[index.u]
  slope.linear.right <- as.numeric((P.right[index.u]-P.right[index.l])/(x.u-x.l))

  P.left <- as.matrix(P.left-1)*slope.poly.left/slope.linear.left+1
  P.right <- as.matrix(P.right-1)*slope.poly.right/slope.linear.right+1

  P.left[xnorm>=min(x),1] <- 0
  P.right[xnorm<=max(x),1] <- 0

  Pnorm[,1] <- Pnorm[,1]+P.left
  Pnorm[,ncol(Pnorm)] <- Pnorm[,ncol(Pnorm)]+P.right

  return(Pnorm)

}

density.deriv.basis <- function(x=NULL,
                                xeval=NULL,
                                xnorm=xnorm,
                                degree=NULL,
                                segments=NULL,
                                basis="tensor",
                                knots="quantiles",
                                monotone=TRUE,
                                deriv.index=1,
                                deriv=1) {

  suppressWarnings(Pnorm.deriv <- prod.spline(x=x,
                                              xeval=xnorm,
                                              K=cbind(degree,segments+if(monotone){2}else{0}),
                                              knots=knots,
                                              basis=basis,
                                              deriv.index=deriv.index,
                                              deriv=deriv))

  if(monotone) Pnorm.deriv <- Pnorm.deriv[,-c(2,degree+segments+1)]

  suppressWarnings(P.lin <- prod.spline(x=x,
                                        xeval=xnorm,
                                        K=cbind(1,1),
                                        knots=knots,
                                        basis=basis,
                                        deriv.index=deriv.index,
                                        deriv=deriv))

  ## For the derivative bases on the extended range `xnorm', above
  ## and below max(x)/min(x) we assign the bases to constants
  ## (zero). We append the linear basis to the left and right of the
  ## bases. The left basis takes on linear values to the left of
  ## min(x), zero elsewhere, the right zero to the left of max(x),
  ## linear elsewhere.

  Pnorm.deriv[xnorm<min(x),] <- 0
  Pnorm.deriv[xnorm>max(x),] <- 0
  P.left <- as.matrix(P.lin[,1])
  P.left[xnorm>=min(x),1] <- 0
  P.right <- as.matrix(P.lin[,2])
  P.right[xnorm<=max(x),1] <- 0
  Pnorm.deriv[,1] <- Pnorm.deriv[,1]+P.left
  Pnorm.deriv[,ncol(Pnorm.deriv)] <- Pnorm.deriv[,ncol(Pnorm.deriv)]+P.right

  return(Pnorm.deriv)

}


clsd <- function(x=NULL,
                 beta=NULL,
                 xeval=NULL,
                 degree=NULL,
                 segments=NULL,
                 degree.min=2,
                 degree.max=25,
                 segments.min=1,
                 segments.max=100,
                 lbound=NULL,
                 ubound=NULL,
                 basis="tensor",
                 knots="quantiles",
                 penalty=NULL,
                 deriv.index=1,
                 deriv=1,
                 elastic.max=TRUE,
                 elastic.diff=3,
                 do.gradient=TRUE,
                 er=NULL,
                 monotone=TRUE,
                 monotone.lb=-250,
                 n.integrate=500,
                 nmulti=1,
                 method = c("L-BFGS-B", "Nelder-Mead", "BFGS", "CG", "SANN"),
                 verbose=FALSE,
                 quantile.seq=seq(.01,.99,by=.01),
                 random.seed=42,
                 maxit=10^5,
                 max.attempts=25,
                 NOMAD=FALSE) {

  if(elastic.max && !NOMAD) {
    degree.max <- 3
    segments.max <- 3
  }

  ptm <- system.time("")

  if(is.null(x)) stop(" You must provide data")

  ## If no er is provided use the following ad-hoc rule which attempts
  ## to ensure we cover the support of the variable for distributions
  ## with moments. This gets the chi-square, t, and Gaussian for n >=
  ## 100 with all degrees of freedom and df=1 is perhaps the worst
  ## case scenario. This rule delivers er = 0.43429448, 0.21714724,
  ## 0.14476483, 0.10857362, 0.08685890, and 0.0723824110, for n = 10,
  ## 10^2, 10^3, 10^4, 10^5, and 10^6. It is probably too aggressive
  ## for the larger samples but one can override - the code traps for
  ## non-finite integration and issues a message when this occurs
  ## along with a suggestion.

  if(!is.null(er) && er < 0) stop(" er must be non-negative")
  if(is.null(er)) er <- 1/log(length(x))

  if(is.null(penalty)) penalty <- log(length(x))/2
  method <- match.arg(method)

  fv <- NULL

  gen.xnorm.out <- gen.xnorm(x=x,
                             lbound=lbound,
                             ubound=ubound,
                             er=er,
                             n.integrate=n.integrate)

  xnorm <- gen.xnorm.out$xnorm
  rank.xnorm <- gen.xnorm.out$rank.xnorm

  if(is.null(beta)) {

    ## If no parameters are provided presume intention is to run
    ## maximum likelihood estimation to obtain the parameter
    ## estimates.

    ptm <- ptm + system.time(ls.ml.out <- ls.ml(x=x,
                                                xnorm=xnorm,
                                                rank.xnorm=rank.xnorm,
                                                degree.min=degree.min,
                                                segments.min=segments.min,
                                                degree.max=degree.max,
                                                segments.max=segments.max,
                                                lbound=lbound,
                                                ubound=ubound,
                                                elastic.max=elastic.max,
                                                elastic.diff=elastic.diff,
                                                do.gradient=do.gradient,
                                                maxit=maxit,
                                                nmulti=nmulti,
                                                er=er,
                                                method=method,
                                                n.integrate=n.integrate,
                                                basis=basis,
                                                knots=knots,
                                                penalty=penalty,
                                                monotone=monotone,
                                                monotone.lb=monotone.lb,
                                                verbose=verbose,
                                                max.attempts=max.attempts,
                                                random.seed=random.seed,
                                                NOMAD=NOMAD))

    beta <- ls.ml.out$beta
    degree <- ls.ml.out$degree
    segments <- ls.ml.out$segments
    fv <- ls.ml.out$fv

  }

  if(!is.null(xeval)) {

    gen.xnorm.out <- gen.xnorm(x=x,
                               xeval=xeval,
                               lbound=lbound,
                               ubound=ubound,
                               er=er,
                               n.integrate=n.integrate)

    xnorm <- gen.xnorm.out$xnorm
    rank.xnorm <- gen.xnorm.out$rank.xnorm

  }

  if(is.null(degree)) stop(" You must provide spline degree")
  if(is.null(segments)) stop(" You must provide number of segments")

  ptm <- ptm + system.time(Pnorm <- density.basis(x=x,
                                                  xeval=xeval,
                                                  xnorm=xnorm,
                                                  degree=degree,
                                                  segments=segments,
                                                  basis=basis,
                                                  knots=knots,
                                                  monotone=monotone))

  if(ncol(Pnorm)!=length(beta)) stop(paste(" Incompatible arguments: beta must be of dimension ",ncol(Pnorm),sep=""))

  Pnorm.beta <- as.numeric(Pnorm%*%as.matrix(beta))

  ## Compute the constant of integration to normalize the density
  ## estimate so that it integrates to one.

  norm.constant <- integrate.trapezoidal.sum(xnorm,exp(Pnorm.beta))
  log.norm.constant <- log(norm.constant)

  if(!is.finite(log.norm.constant))
    stop(paste(" integration not finite - perhaps try reducing `er' (current value = ",round(er,3),")",sep=""))

  ## For the distribution, compute the density over the extended
  ## range, then return values corresponding to either the sample x or
  ## evaluation x (xeval) based on integration over the extended range
  ## for the xnorm points (xnorm contains x and xeval - this ought to
  ## ensure integration to one).

  ## f.norm is the density evaluated on the extended range (including
  ## sample observations and evaluation points if the latter exist),
  ## F.norm the distribution evaluated on the extended range.

  f.norm <- exp(Pnorm.beta-log.norm.constant)
  F.norm <- integrate.trapezoidal(xnorm,f.norm)

  if(deriv > 0) {

    ptm <- ptm + system.time(Pnorm.deriv <- density.deriv.basis(x=x,
                                                                xeval=xeval,
                                                                xnorm=xnorm,
                                                                degree=degree,
                                                                segments=segments,
                                                                basis=basis,
                                                                knots=knots,
                                                                monotone=monotone,
                                                                deriv.index=deriv.index,
                                                                deriv=deriv))

    f.norm.deriv <- as.numeric(f.norm*Pnorm.deriv%*%beta)

  } else {

    f.deriv <- NULL
    f.norm.deriv <- NULL

  }

  ## Compute quantiles using the the quasi-inverse (Definition 2.3.6,
  ## Nelson (2006))

  quantile.vec <- numeric(length(quantile.seq))
  for(i in 1:length(quantile.seq)) {
    if(quantile.seq[i]>=0.5) {
      quantile.vec[i] <- max(xnorm[F.norm<=quantile.seq[i]])
    } else {
      quantile.vec[i] <-  min(xnorm[F.norm>=quantile.seq[i]])
    }
  }

  ## Next, strip off the values of the distribution corresponding to
  ## either sample x or evaluation xeval

  if(is.null(xeval)) {
    f <-   f.norm[rank.xnorm][1:length(x)]
    F <-   F.norm[rank.xnorm][1:length(x)]
    f.norm <- f.norm[rank.xnorm][(length(x)+1):length(f.norm)]
    F.norm <- F.norm[rank.xnorm][(length(x)+1):length(F.norm)]
    xnorm <- xnorm[rank.xnorm][(length(x)+1):length(xnorm)]
    if(deriv>0) {
      f.deriv <- f.norm.deriv[rank.xnorm][1:length(x)]
      f.norm.deriv <- f.norm.deriv[rank.xnorm][(length(x)+1):length(f.norm.deriv)]
    }
    P <-   Pnorm[rank.xnorm,][1:length(x),]
    P.beta <- Pnorm.beta[rank.xnorm][1:length(x)]
  } else {
    f <-   f.norm[rank.xnorm][1:length(xeval)]
    F <-   F.norm[rank.xnorm][1:length(xeval)]
    f.norm <- f.norm[rank.xnorm][(length(x)+length(xeval)+1):length(f.norm)]
    F.norm <- F.norm[rank.xnorm][(length(x)+length(xeval)+1):length(F.norm)]
    xnorm <- xnorm[rank.xnorm][(length(x)+length(xeval)+1):length(xnorm)]
    if(deriv>0) {
      f.deriv <- f.norm.deriv[rank.xnorm][1:length(xeval)]
      f.norm.deriv <- f.norm.deriv[rank.xnorm][(length(x)+length(xeval)+1):length(f.norm.deriv)]
    }
    P <-   Pnorm[rank.xnorm,][1:length(xeval),]
    P.beta <- Pnorm.beta[rank.xnorm][1:length(xeval)]
  }

  clsd.return <- list(density=f,
                      density.deriv=f.deriv,
                      distribution=F,
                      density.er=f.norm,
                      density.deriv.er=f.norm.deriv,
                      distribution.er=F.norm,
                      xer=xnorm,
                      Basis.beta=P.beta,
                      Basis.beta.er=Pnorm.beta,
                      P=P,
                      Per=Pnorm,
                      logl=sum(P.beta-log.norm.constant),
                      constant=norm.constant,
                      degree=degree,
                      segments=segments,
                      knots=knots,
                      basis=basis,
                      nobs=length(x),
                      beta=beta,
                      fv=fv,
                      er=er,
                      penalty=penalty,
                      nmulti=nmulti,
                      x=x,
                      xq=quantile.vec,
                      tau=quantile.seq,
                      ptm=ptm)

  class(clsd.return) <- "clsd"
  return(clsd.return)

}

sum.log.density <- function(beta=NULL,
                            P=NULL,
                            Pint=NULL,
                            xint=NULL,
                            length.x=NULL,
                            penalty=NULL,
                            complexity=NULL,
                            ...) {

  return(2*sum(P%*%beta)-2*length.x*log(integrate.trapezoidal.sum(xint,exp(Pint%*%beta)))-penalty*complexity)

}

sum.log.density.gradient <- function(beta=NULL,
                                     colSumsP=NULL,
                                     Pint=NULL,
                                     xint=NULL,
                                     length.x=NULL,
                                     penalty=NULL,
                                     complexity=NULL,
                                     ...) {


  exp.Pint.beta <- as.numeric(exp(Pint%*%beta))
  exp.Pint.beta.Pint <- exp.Pint.beta*Pint
  int.exp.Pint.beta.Pint <- numeric()
  for(i in 1:complexity) int.exp.Pint.beta.Pint[i] <- integrate.trapezoidal.sum(xint,exp.Pint.beta.Pint[,i])
  return(2*(colSumsP-length.x*int.exp.Pint.beta.Pint/integrate.trapezoidal.sum(xint,exp.Pint.beta)))

}

ls.ml <- function(x=NULL,
                  xnorm=NULL,
                  rank.xnorm=NULL,
                  degree.min=NULL,
                  segments.min=NULL,
                  degree.max=NULL,
                  segments.max=NULL,
                  lbound=NULL,
                  ubound=NULL,
                  elastic.max=FALSE,
                  elastic.diff=NULL,
                  do.gradient=TRUE,
                  maxit=NULL,
                  nmulti=NULL,
                  er=NULL,
                  method=NULL,
                  n.integrate=NULL,
                  basis=NULL,
                  knots=NULL,
                  penalty=NULL,
                  monotone=TRUE,
                  monotone.lb=NULL,
                  verbose=NULL,
                  max.attempts=NULL,
                  random.seed=NULL,
                  NOMAD=FALSE) {

  ## This function conducts log spline maximum
  ## likelihood. Multistarting is supported as is breaking out to
  ## potentially avoid wasted computation (be careful when using this,
  ## however, as it is prone to stopping early).

  ## Save seed prior to setting

  if(exists(".Random.seed", .GlobalEnv)) {
    save.seed <- get(".Random.seed", .GlobalEnv)
    exists.seed = TRUE
  } else {
    exists.seed = FALSE
  }

  set.seed(random.seed)

  if(missing(x)) stop(" You must provide data")

  if(!NOMAD) {

    ## We set some initial parameters that are placeholders to get
    ## things rolling.
    
    d.opt <- Inf
    s.opt <- Inf
    par.opt <- Inf
    value.opt <- -Inf
    length.x <- length(x)
    length.xnorm <- length(xnorm)
    
    ## Loop through all degrees for every segment starting at
    ## segments.min.
    
    d <- degree.min
    
    while(d <= degree.max) {
      
      ## For smooth densities one can simply restrict degree to at least
      ## 2 (or 3 to be consistent with cubic splines)
      
      s <- segments.min
      
      while(s <= segments.max) {
        
        if(verbose) cat("\n")
        cat("\r                                                                                                  ")
        cat("\rOptimizing, degree = ",d,", segments = ",s,", degree.opt = ",d.opt, ", segments.opt = ",s.opt," ",sep="")
        
        ## Generate objects that need not be recomputed for a given d
        ## and s
        
        Pnorm <- density.basis(x=x,
                               xnorm=xnorm,
                               degree=d,
                               segments=s,
                               basis=basis,
                               knots=knots,
                               monotone=monotone)
        
        P <- Pnorm[rank.xnorm,][1:length.x,]
        colSumsP <- colSums(P)
        Pint <- Pnorm[rank.xnorm,][(length.x+1):nrow(Pnorm),]
        xint <- xnorm[rank.xnorm][(length.x+1):nrow(Pnorm)]
        
        complexity <- d+s
        
        ## Multistart if desired.
        
        for(n in 1:nmulti) {
          
          ## Can restart to see if we can improve on min... note initial
          ## values totally ad-hoc...
          
          par.init.out <- par.init(d,s,monotone,monotone.lb)
          par.init <- par.init.out$par.init
          par.upper <- par.init.out$par.upper
          par.lower <- par.init.out$par.lower
          
          ## Trap non-convergence, restart from different initial
          ## points, display message if needed (trace>0 up to 6 provides
          ## ever more detailed information for L-BFGS-B)
          
          optim.out <- list()
          optim.out[[4]] <- 9999
          optim.out$value <- -Inf
          
          m.attempts <- 0
          
          while(tryCatch(suppressWarnings(optim.out <- optim(par=par.init,
                                                             fn=sum.log.density,
                                                             gr=if(do.gradient){sum.log.density.gradient}else{NULL},
                                                             upper=par.upper,
                                                             lower=par.lower,
                                                             method=method,
                                                             penalty=penalty,
                                                             P=P,
                                                             colSumsP=colSumsP,
                                                             Pint=Pint,
                                                             length.x=length.x,
                                                             xint=xint,
                                                             complexity=complexity,
                                                             control=list(fnscale=-1,maxit=maxit,if(verbose){trace=1}else{trace=0}))),
                         error = function(e){return(optim.out)})[[4]]!=0 && m.attempts < max.attempts){
            
            ## If optim fails to converge, reset initial parameters and
            ## try again.
            
            if(verbose && optim.out[[4]]!=0) {
              if(!is.null(optim.out$message)) cat("\n optim message = ",optim.out$message,sep="")
              cat("\n optim failed (degree = ",d,", segments = ",s,", convergence = ", optim.out[[4]],") re-running with new initial values",sep="")
            }
            
            par.init.out <- par.init(d,s,monotone,monotone.lb)
            par.init <- par.init.out$par.init
            par.upper <- par.init.out$par.upper
            par.lower <- par.init.out$par.lower
            
            m.attempts <- m.attempts+1
            
          }
          
          ## Check for a new optimum, overwrite existing values with
          ## new values.
          
          if(optim.out$value > value.opt) {
            if(verbose && n==1) cat("\n optim improved: d = ",d,", s = ",s,", old = ",formatC(value.opt,format="g",digits=6),", new = ",formatC(optim.out$value,format="g",digits=6),", diff = ",formatC(optim.out$value-value.opt,format="g",digits=6),sep="")
            
            if(verbose && n>1) cat("\n optim improved (ms ",n,"/",nmulti,"): d = ",d,", s = ",s,", old = ",formatC(value.opt,format="g",digits=6),", new = ",formatC(optim.out$value,format="g",digits=6),", diff = ",formatC(optim.out$value-value.opt,format="g",digits=6),sep="")
            
            par.opt <- optim.out$par
            d.opt <- d
            s.opt <- s
            value.opt <- optim.out$value
          }
          
        }
        
        if(!(segments.min==segments.max) && elastic.max && s.opt == segments.max) segments.max <- segments.max+elastic.diff
        if(!(segments.min==segments.max) && elastic.max && s.opt < segments.max+elastic.diff) segments.max <- s.opt+elastic.diff
        
        s <- s+1
        
      }
      
      d <- d+1
      
      if(!(degree.min==degree.max) && elastic.max && d.opt == degree.max) degree.max <- degree.max+elastic.diff
      if(!(degree.min==degree.max) && elastic.max && d.opt < degree.max+elastic.diff) degree.max <- d.opt+elastic.diff
      
    }
    
  } else {

    eval.f <- function(input, params) {

      sum.log.density <- params$sum.log.density
      sum.log.density.gradient <- params$sum.log.density.gradient
      method <- params$method
      penalty <- params$penalty
      x <- params$x
      xnorm <- params$xnorm
      knots <- params$knots
      basis <- params$basis
      monotone <- params$monotone
      monotone.lb <- params$monotone.lb
      rank.xnorm <- params$rank.xnorm
      do.gradient <- params$do.gradient
      maxit <- params$maxit
      max.attempts <- params$max.attempts
      verbose <- params$verbose

      length.x <- length(x)
      length.xnorm <- length(xnorm)

      d <- input[1]
      s <- input[2]
      
      complexity <- d+s
        
      Pnorm <- density.basis(x=x,
                             xnorm=xnorm,
                             degree=d,
                             segments=s,
                             basis=basis,
                             knots=knots,
                             monotone=monotone)
      
      P <- Pnorm[rank.xnorm,][1:length.x,]
      colSumsP <- colSums(P)
      Pint <- Pnorm[rank.xnorm,][(length.x+1):nrow(Pnorm),]
      xint <- xnorm[rank.xnorm][(length.x+1):nrow(Pnorm)]
      
      ## NOMAD minimizes only

      optim.out <- list()
      optim.out[[4]] <- 9999
      optim.out$value <- -Inf
      
      m.attempts <- 0
          
      while(tryCatch(suppressWarnings(optim.out <- optim(par=par.init,
                                                         fn=sum.log.density,
                                                         gr=if(do.gradient){sum.log.density.gradient}else{NULL},
                                                         upper=par.upper,
                                                         lower=par.lower,
                                                         method=method,
                                                         penalty=penalty,
                                                         P=P,
                                                         colSumsP=colSumsP,
                                                         Pint=Pint,
                                                         length.x=length.x,
                                                         xint=xint,
                                                         complexity=complexity,
                                                         control=list(fnscale=-1,maxit=maxit,if(verbose){trace=1}else{trace=0}))),
                     error = function(e){return(optim.out)})[[4]]!=0 && m.attempts < max.attempts){
        
        ## If optim fails to converge, reset initial parameters and
        ## try again.
        
        if(verbose && optim.out[[4]]!=0) {
          if(!is.null(optim.out$message)) cat("\n optim message = ",optim.out$message,sep="")
          cat("\n optim failed (degree = ",d,", segments = ",s,", convergence = ", optim.out[[4]],") re-running with new initial values",sep="")
        }
        
        par.init.out <- par.init(d,s,monotone,monotone.lb)
        par.init <- par.init.out$par.init
        par.upper <- par.init.out$par.upper
        par.lower <- par.init.out$par.lower
        
        m.attempts <- m.attempts+1
        
      }

      if(verbose) cat("\n")
      cat("\r                                                                                                  ")
      cat("\rOptimizing, degree = ",d,", segments = ",s,", log likelihood = ",optim.out$value,sep="")

      fv <- -optim.out$value
      
    }

    ## Initial values
    x0 <- c(degree.min,segments.min)
    ## Types of variables
    bbin <-c(1, 1)
    ## Bounds
    lb <- c(degree.min,segments.min)
    ub <- c(degree.max,segments.max)
    ## Type of output
    bbout <- c(0, 2, 1)
    ## Options
    opts <-list("MAX_BB_EVAL"=10000)

    ## Generate params

    params <- list()

    params$sum.log.density <- sum.log.density
    params$sum.log.density.gradient <- sum.log.density.gradient
    params$method <- method
    params$penalty <- penalty
    params$x <- x
    params$xnorm <- xnorm
    params$knots <- knots
    params$basis <- basis    
    params$monotone <- monotone
    params$monotone.lb <- monotone.lb
    params$rank.xnorm <- rank.xnorm
    params$do.gradient <- do.gradient
    params$maxit <- maxit
    params$max.attempts <- max.attempts
    params$verbose <- verbose            

    solution <- snomadr(eval.f=eval.f,
                        n=2,## number of variables
                        x0=x0,
                        bbin=bbin,
                        bbout=bbout,
                        lb=lb,
                        ub=ub,
                        nmulti=nmulti,
                        print.output=FALSE,
                        opts=opts,
                        params=params)

    value.opt <- solution$objective

    d <- solution$solution[1]
    s <- solution$solution[2]

    ## Final call to optim to retrieve beta
    
    length.x <- length(x)
    length.xnorm <- length(xnorm)
    
    complexity <- d+s
    
    par.init.out <- par.init(d,s,monotone,monotone.lb)
    par.init <- par.init.out$par.init
    par.upper <- par.init.out$par.upper
    par.lower <- par.init.out$par.lower
    
    Pnorm <- density.basis(x=x,
                           xnorm=xnorm,
                           degree=d,
                           segments=s,
                           basis=basis,
                           knots=knots,
                           monotone=monotone)
    
    P <- Pnorm[rank.xnorm,][1:length.x,]
    colSumsP <- colSums(P)
    Pint <- Pnorm[rank.xnorm,][(length.x+1):nrow(Pnorm),]
    xint <- xnorm[rank.xnorm][(length.x+1):nrow(Pnorm)]
    
    optim.out <- list()
    optim.out[[4]] <- 9999
    optim.out$value <- -Inf
    
    m.attempts <- 0
    
    while(tryCatch(suppressWarnings(optim.out <- optim(par=par.init,
                                                       fn=sum.log.density,
                                                       gr=if(do.gradient){sum.log.density.gradient}else{NULL},
                                                       upper=par.upper,
                                                       lower=par.lower,
                                                       method=method,
                                                       penalty=penalty,
                                                       P=P,
                                                       colSumsP=colSumsP,
                                                       Pint=Pint,
                                                       length.x=length.x,
                                                       xint=xint,
                                                       complexity=complexity,
                                                       control=list(fnscale=-1,maxit=maxit,if(verbose){trace=1}else{trace=0}))),
                   error = function(e){return(optim.out)})[[4]]!=0 && m.attempts < max.attempts){
      
      ## If optim fails to converge, reset initial parameters and
      ## try again.
      
      if(verbose && optim.out[[4]]!=0) {
        if(!is.null(optim.out$message)) cat("\n optim message = ",optim.out$message,sep="")
        cat("\n optim failed (degree = ",d,", segments = ",s,", convergence = ", optim.out[[4]],") re-running with new initial values",sep="")
      }
      
      par.init.out <- par.init(d,s,monotone,monotone.lb)
      par.init <- par.init.out$par.init
      par.upper <- par.init.out$par.upper
      par.lower <- par.init.out$par.lower
      
      m.attempts <- m.attempts+1
      
    }
    
    d.opt <- d
    s.opt <- s
    par.opt <- optim.out$par
    value.opt <- optim.out$value
    
  }
  
  cat("\r                                                                            ")
  if(!(degree.min==degree.max) && (d.opt==degree.max)) warning(paste(" optimal degree equals search maximum (", d.opt,"): rerun with larger degree.max",sep=""))
  if(!(segments.min==segments.max) && (s.opt==segments.max)) warning(paste(" optimal segment equals search maximum (", s.opt,"): rerun with larger segments.max",sep=""))
  if(par.opt[1]>0|par.opt[length(par.opt)]>0) warning(" optim() delivered a positive weight for linear segment (supposed to be negative)")
  if(!monotone&&par.opt[1]<=monotone.lb) warning(paste(" optimal weight for left nonmonotone basis equals search minimum (",par.opt[1],"): rerun with smaller monotone.lb",sep=""))
  if(!monotone&&par.opt[length(par.opt)]<=monotone.lb) warning(paste(" optimal weight for right nonmonotone basis equals search minimum (",par.opt[length(par.opt)],"): rerun with smaller monotone.lb",sep=""))
  
  ## Restore seed

  if(exists.seed) assign(".Random.seed", save.seed, .GlobalEnv)

  return(list(degree=d.opt,segments=s.opt,beta=par.opt,fv=value.opt))

}

summary.clsd <- function(object,
                         ...) {

  cat("\nCategorical Logspline Density\n",sep="")
  cat(paste("\nModel penalty: ", format(object$penalty), sep=""))
  cat(paste("\nModel degree/segments: ", format(object$degree),"/",format(object$segments), sep=""))
  cat(paste("\nKnot type: ", format(object$knots), sep=""))
  cat(paste("\nBasis type: ",format(object$basis),sep=""))
  cat(paste("\nTraining observations: ", format(object$nobs), sep=""))
  cat(paste("\nLog-likelihood: ", format(object$logl), sep=""))
  cat(paste("\nNumber of multistarts: ", format(object$nmulti), sep=""))
  cat(paste("\nEstimation time: ", formatC(object$ptm[1],digits=1,format="f"), " seconds",sep=""))

  cat("\n\n")

}

print.clsd <- function(x,...)
{
        summary.clsd(x)
}

plot.clsd <- function(x,
                      er=TRUE,
                      distribution=FALSE,
                      derivative=FALSE,
                      ylim,
                      ylab,
                      xlab,
                      type,
                      ...) {

  if(missing(xlab)) xlab <- "Data"
  if(missing(type)) type <- "l"

  if(!er) {
    order.x <- order(x$x)
    if(distribution){
      y <- x$distribution[order.x]
      if(missing(ylab)) ylab <- "Distribution"
      if(missing(ylim)) ylim <- c(0,1)
    }
    if(!distribution&&!derivative) {
      y <- x$density[order.x]
      if(missing(ylab)) ylab <- "Density"
      if(missing(ylim)) ylim <- c(0,max(y))
    }
    if(derivative) {
      y <- x$density.deriv[order.x]
      if(missing(ylab)) ylab <- "Density Derivative"
      if(missing(ylim)) ylim <- c(min(y),max(y))
    }
    x <- x$x[order.x]
  } else {
    order.xer <- order(x$xer)
    if(distribution){
      y <- x$distribution.er[order.xer]
      if(missing(ylab)) ylab <- "Distribution"
      if(missing(ylim)) ylim <- c(0,1)      
    }
    if(!distribution&&!derivative) {
      y <- x$density.er[order.xer]
      if(missing(ylab)) ylab <- "Density"
      if(missing(ylim)) ylim <- c(0,max(y))      
    }
    if(derivative){
      y <- x$density.deriv.er[order.xer]
      if(missing(ylab)) ylab <- "Density Derivative"
      if(missing(ylim)) ylim <- c(min(y),max(y))      
    }
    x <- x$xer[order.xer]
  }

  x <- plot(x,
            y,
            ylim=ylim,
            ylab=ylab,
            xlab=xlab,
            type=type,
            ...)

}

coef.clsd <- function(object, ...) {
  tc <- object$beta
  return(tc)
}

fitted.clsd <- function(object, ...){
 object$density
}

