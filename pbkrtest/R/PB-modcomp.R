##########################################################
###
### Bartlett corrected LRT
###
##########################################################

PBmodcomp <- function(largeModel, smallModel, nsim=1000, ref=NULL, seed=NULL, cl=NULL, details=0){
  UseMethod("PBmodcomp")
}

PBmodcomp.merMod <-
PBmodcomp.mer <-
    function(largeModel, smallModel, nsim=1000, ref=NULL, seed=NULL, cl=NULL, details=0){
        
        ##cat("PBmodcomp.lmerMod\n")
        f.large <- formula(largeModel)
        attributes(f.large) <- NULL
        
        if (inherits(smallModel, c("Matrix", "matrix"))){
            f.small <- smallModel
            smallModel <- restrictionMatrix2model(largeModel, smallModel)
        } else {
            f.small <- formula(smallModel)
            attributes(f.small) <- NULL
        }
        
        if (is.null(ref)){
            ref <- PBrefdist(largeModel, smallModel, nsim=nsim, seed=seed, cl=cl, details=details)
        }
        
        #' samples <- attr(ref, "samples")
        #' if (!is.null(samples)){
        #'     nsim <- samples['nsim']
        #'     npos <- samples['npos']
        #' } else {
        #'     nsim <- length(ref)
        #'     npos <- sum(ref>0)
        #' }
        
        LRTstat     <- getLRT(largeModel, smallModel)
        ans         <- .finalizePB(LRTstat, ref)
        .padPB( ans, LRTstat, ref, f.large, f.small)
    }

.padPB <- function(ans, LRTstat, ref, f.large, f.small){
    ans$LRTstat <- LRTstat
    ans$ref     <- ref
    ans$f.large <- f.large
    ans$f.small <- f.small
    ans
}

PBmodcomp.lm <- function(largeModel, smallModel, nsim=1000, ref=NULL, seed=NULL, cl=NULL, details=0){

  ok.fam <- c("binomial","gaussian","Gamma","inverse.gaussian","poisson")
  f.large <- formula(largeModel)
  attributes(f.large) <- NULL

  if (inherits(smallModel, c("Matrix", "matrix"))){
    f.small <- smallModel
    smallModel <- restrictionMatrix2model(largeModel, smallModel)
  } else {
    f.small <- formula(smallModel)
    attributes(f.small) <- NULL
  }

  if (!all.equal((fam.l <- family(largeModel)), (fam.s <- family(smallModel))))
    stop("Models do not have identical identical family\n")
  if (!(fam.l$family %in% ok.fam)){
    stop(sprintf("family must be of type %s", toString(ok.fam)))
  }

  if (is.null(ref)){
    ref <- PBrefdist(largeModel, smallModel, nsim=nsim, seed=seed, cl=cl, details=details)
  }

  LRTstat     <- getLRT(largeModel, smallModel)
  ans         <- .finalizePB(LRTstat, ref)
    .padPB( ans, LRTstat, ref, f.large, f.small)    
}

.finalizePB <- function(LRTstat, ref){

  tobs <- unname(LRTstat[1])
  ndf  <- unname(LRTstat[2])

  refpos <- ref[ref>0]
  nsim <- length(ref)
  npos <- length(refpos)

  ##cat(sprintf("EE=%f VV=%f\n", EE, VV))
  p.chi <- 1-pchisq(tobs, df=ndf)
  ## Direct computation of tail probability

  n.extreme <- sum(tobs < refpos)
  p.PB  <- (1+n.extreme) / (1+npos)

  test = list(
    LRT      = c(stat=tobs,    df=ndf,    p.value=p.chi),
    PBtest   = c(stat=tobs,    df=NA,     p.value=p.PB))

  test  <- as.data.frame(do.call(rbind, test))
  ans   <- list(test=test, type="X2test", samples=c(nsim=nsim, npos=npos), n.extreme=n.extreme,
                ctime=attr(ref,"ctime"))
  class(ans) <- c("PBmodcomp")
  ans
}

.summarizePB <- function(LRTstat, ref){

  tobs <- unname(LRTstat[1])
  ndf  <- unname(LRTstat[2])

  refpos   <- ref[ref>0]
  nsim <- length(ref)
  npos <- length(refpos)


  EE      <- mean(refpos)
  VV      <- var(refpos)

  ##cat(sprintf("EE=%f VV=%f\n", EE, VV))
  p.chi <- 1-pchisq(tobs, df=ndf)

  ## Direct computation of tail probability
  n.extreme <- sum(tobs < refpos)
  ##p.PB  <- n.extreme / npos
  p.PB  <- (1+n.extreme) / (1+npos)

  p.PB.all  <- (1+n.extreme) / (1+nsim)


  se <- round(sqrt(p.PB*(1-p.PB)/npos),4)
  ci <- round(c(-1.96, 1.96)*se + p.PB,4)

  ## Kernel density estimate
  ##dd <- density(ref)
  ##p.KD <- sum(dd$y[dd$x>=tobs])/sum(dd$y)

  ## Bartlett correction - X2 distribution
  BCstat  <- ndf * tobs/EE
  ##cat(sprintf("BCval=%f\n", ndf/EE))
  p.BC    <- 1-pchisq(BCstat,df=ndf)

  ## Fit to gamma distribution
  scale   <- VV/EE
  shape   <- EE^2/VV
  p.Ga    <- 1-pgamma(tobs, shape=shape, scale=scale)

  ## Fit T/d to F-distribution (1. moment)

  ## FIXME: Think the formula is 2*EE/(EE-1)
  ##ddf  <- 2*EE/(EE-ndf)
  ddf  <- 2*EE/(EE-1)
  Fobs <- tobs/ndf
  if (ddf>0)
      p.FF <- 1-pf(Fobs, df1=ndf, df2=ddf)
  else
      p.FF <- NA

  ## Fit T/d to F-distribution (1. AND 2. moment)
  #' EE2   <- EE/ndf
  #' VV2   <- VV/ndf^2

  #' rho   <- VV2/(2*EE2^2)
  #' ddf2  <- 4 + (ndf+2)/(rho*ndf-1)
  #' lam2  <- (ddf/EE2*(ddf-2))
  #' Fobs2 <- lam2 * tobs/ndf
  #' if (ddf2>0)
  #'   p.FF2 <- 1-pf(Fobs2, df1=ndf, df2=ddf2)
  #' else
  #'   p.FF2 <- NA

  #' cat(sprintf("PB: EE=%f, ndf=%f VV=%f, ddf=%f\n", EE, ndf, VV, ddf))




  test = list(
    PBtest   = c(stat=tobs,    df=NA,  ddf=NA,   p.value=p.PB),
    Gamma    = c(stat=tobs,    df=NA,  ddf=NA,   p.value=p.Ga),
    Bartlett = c(stat=BCstat,  df=ndf, ddf=NA,   p.value=p.BC),
    F        = c(stat=Fobs,    df=ndf, ddf=ddf,  p.value=p.FF),
    LRT      = c(stat=tobs,    df=ndf, ddf=NA,   p.value=p.chi)
    )
    ##          PBkd     = c(stat=tobs,    df=NA,  ddf=NA,   p.value=p.KD),

    ##F2       = c(stat=Fobs2,   df=ndf, ddf=ddf2, p.value=p.FF2),
                                        #,
    #PBtest.all   = c(stat=tobs,    df=NA,  ddf=NA,   p.value=p.PB.all),
    #Bartlett.all = c(stat=BCstat.all,  df=ndf, ddf=NA,   p.value=p.BC.all)
    ##F2       = c(stat=Fobs2,   df=ndf,  p.value=p.FF2, ddf=ddf2)


  test <- as.data.frame(do.call(rbind, test))
  ans <- list(test=test, type="X2test",
              moment = c(mean=EE, var=VV),
              samples= c(nsim=nsim, npos=npos),
              gamma  = c(scale=scale, shape=shape),
              ref    = ref,
              ci     = ci,
              se     = se,
              n.extreme = n.extreme,
              ctime  = attr(ref, "ctime")
              )
  class(ans) <- c("PBmodcomp")
  ans
}



##   rho   <- VV/(2*EE^2)
##   ddf2  <- (ndf*(4*rho+1) - 2)/(rho*ndf-1)
##   lam2  <- (ddf/(ddf-2))/(EE/ndf)
##   cat(sprintf("EE=%f, VV=%f, rho=%f, lam2=%f\n",
##               EE, VV, rho, lam2))

##   ddf2 <- 4 + (ndf+2)/(rho*ndf-1)

##   Fobs2 <- lam2 * tobs/ndf
##   if (ddf2>0)
##     p.FF2 <- 1-pf(Fobs2, df1=ndf, df2=ddf2)
##   else
##     p.FF2 <- NA


### ###########################################################
###
### Utilities
###
### ###########################################################

.PBcommon <- function(x){

    cat(sprintf("Parametric bootstrap test; "))
    if (!is.null((zz<- x$ctime))){
        cat(sprintf("time: %.2f sec; ", round(zz,2)))
    }
    if (!is.null((sam <- x$samples))){
        cat(sprintf("samples: %d extremes: %d;", sam[1], x$n.extreme))
    }
    cat("\n")
    
    
    if (!is.null((sam <- x$samples))){
        if (sam[2]<sam[1]){
            cat(sprintf("Requested samples: %d Used samples: %d Extremes: %d\n",
                        sam[1], sam[2], x$n.extreme))
        }
    }
    if(!is.null(x$f.large)){
        cat("large : "); print(x$f.large)
        cat("small : "); print(x$f.small)
    }
}



print.PBmodcomp <- function(x, ...){
  .PBcommon(x)
  tab <- x$test
  printCoefmat(tab, tst.ind=1, na.print='', has.Pvalue=TRUE)
  return(invisible(x))
}


summary.PBmodcomp <- function(object,...){
  ans <- .summarizePB(object$LRTstat, object$ref)
  ans$f.large <- object$f.large
  ans$f.small <- object$f.small
  class(ans) <- "summaryPB"
  ans
}

print.summaryPB <- function(x,...){
  .PBcommon(x)
  ans <- x$test
  printCoefmat(ans, tst.ind=1, na.print='', has.Pvalue=TRUE)
  cat("\n")

##   ci <- x$ci
##   cat(sprintf("95 pct CI for PBtest   : [%s]\n", toString(ci)))
##   mo <- x$moment
##   cat(sprintf("Reference distribution : mean=%f var=%f\n", mo[1], mo[2]))
##   ga <- x$gamma
##   cat(sprintf("Gamma approximation    : scale=%f shape=%f\n", ga[1], ga[2]))

  return(invisible(x))
}


plot.PBmodcomp <- function(x, ...){

  ref <-x$ref

  ndf  <- x$test$df[1]
  u    <-summary(x)
  ddf  <-u$test['F','ddf']

  EE   <- mean(ref)
  VV   <- var(ref)
  sc   <- var(ref)/mean(ref)
  sh   <- mean(ref)^2/var(ref)
  sc   <- VV/EE
  sh   <- EE^2/VV
  B    <- ndf/EE # if ref is the null distr, so should A*ref follow a chisq(df=ndf) distribution

  upper <- 0.20
  #tail.prob <- c(0.0001, 0.001, 0.01, 0.05, 0.10, 0.20, 0.5)
  tail.prob <-seq(0.001, upper, length.out = 1111)
  PBquant   <- quantile(ref,1-tail.prob) ## tail prob for PB dist

  pLR       <- pchisq(PBquant,df=ndf,           lower.tail=FALSE)
  pF        <- pf(PBquant/ndf,df1=ndf,df2=ddf,  lower.tail=FALSE)
  pGamma    <- pgamma(PBquant,scale=sc,shape=sh,lower.tail=FALSE)
  pBart     <- pchisq(B*PBquant,df=ndf,         lower.tail=FALSE)

  sym.vec <- c(2,4,5,6)
  lwd     <- 2
  plot(pLR~tail.prob,type='l', lwd=lwd, #log="xy",
       xlab='Nominal p-value',ylab='True p-value',
       xlim=c(1e-3, upper),ylim=c(1e-3, upper),
       col=sym.vec[1], lty=sym.vec[1])
  lines(pF~tail.prob,lwd=lwd,     col=sym.vec[2], lty=sym.vec[2])
  lines(pBart~tail.prob,lwd=lwd,  col=sym.vec[3], lty=sym.vec[3])
  lines(pGamma~tail.prob,lwd=lwd, col=sym.vec[4], lty=sym.vec[4])
  abline(c(0,1))

  ZF     <-bquote(paste("F(1,",.(round(ddf,2)),")"))
  Zgamma <-bquote(paste("gamma(scale=",.(round(sc,2)),", shape=", .(round(sh,2)),")" ))
  ZLRT   <-bquote(paste(chi[.(ndf)]^2))
  ZBart  <-bquote(paste("Bartlett scaled ", chi[.(ndf)]^2))

  legend(0.001,upper,legend=as.expression(c(ZLRT,ZF,ZBart,Zgamma)),
         lty=sym.vec,col=sym.vec,lwd=lwd)
}
as.data.frame.XXmodcomp <- function(x, row.names = NULL, optional = FALSE, ...){
    as.data.frame(do.call(rbind, x[-c(1:3)]))
}




## plot.XXmodcomp <- function(x, ...){

##   test <- x$test
##   tobs <- test$LRT['stat']
##   ref <- attr(x,"ref")
##   rr  <- range(ref)
##   xx  <- seq(rr[1],rr[2],0.1)
##   dd  <- density(ref)
##   sc  <- var(ref)/mean(ref)
##   sh  <- mean(ref)^2/var(ref)

##   hist(ref, prob=TRUE,nclass=20, main="Reference distribution")
##   abline(v=tobs)
##   lines(dd, lty=2, col=2, lwd=2)
##   lines(xx,dchisq(xx,df=test$LRT['df']), lty=3, col=3, lwd=2)
##   lines(xx,dgamma(xx,scale=sc, shape=sh), lty=4, col=4, lwd=2)
##   lines(xx,df(xx,df1=test$F['df'], df2=test$F['ddf']), lty=5, col=5, lwd=2)

##   smartlegend(x = 'right', y = 'top',
##               legend = c("kernel density", "chi-square", "gamma","F"),
##               col = 2:5, lty = 2:5)

## }




##   rho   <- VV/(2*EE^2)
##   ddf2  <- (ndf*(4*rho+1) - 2)/(rho*ndf-1)
##   lam2  <- (ddf/(ddf-2))/(EE/ndf)
##   cat(sprintf("EE=%f, VV=%f, rho=%f, lam2=%f\n",
##               EE, VV, rho, lam2))

##   ddf2 <- 4 + (ndf+2)/(rho*ndf-1)

##   Fobs2 <- lam2 * tobs/ndf
##   if (ddf2>0)
##     p.FF2 <- 1-pf(Fobs2, df1=ndf, df2=ddf2)
##   else
##     p.FF2 <- NA
