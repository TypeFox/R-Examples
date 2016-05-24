KRmodcomp <- function(largeModel, smallModel,betaH=0, details=0){
    UseMethod("KRmodcomp")
}

KRmodcomp.lmerMod<-
    KRmodcomp.mer<-
        function(largeModel, smallModel, betaH=0, details=0) {
  ## 'smallModel' can either be an lmerMod (linear mixed) model or a restriction matrix L.
  w <- KRmodcomp_init(largeModel, smallModel, matrixOK = TRUE)
  if (w==-1) {
    stop('Models have either equal fixed mean stucture or are not nested')
  } else {
    if (w==0){
      ##stop('First given model is submodel of second; exchange the models\n')
      tmp <- largeModel
      largeModel <- smallModel
      smallModel <- tmp
    }
  }

  ## Refit large model with REML if necessary
  if (!(getME(largeModel, "is_REML"))){
    largeModel <- update(largeModel,.~.,REML=TRUE)
  }



  ## All computations are based on 'largeModel' and the restriction matrix 'L'
  ## -------------------------------------------------------------------------
  t0    <- proc.time()
  L     <- .model2restrictionMatrix(largeModel, smallModel)

  PhiA  <- vcovAdj(largeModel, details)
  stats <- .KR_adjust(PhiA, Phi=vcov(largeModel), L, beta=fixef(largeModel), betaH)
  stats <- lapply(stats, c) ## To get rid of all sorts of attributes
  ans   <- .finalizeKR(stats)

  f.small <-
    if (.is.lmm(smallModel)){
      .zzz <- formula(smallModel)
      attributes(.zzz) <- NULL
      .zzz
    } else {
      list(L=L, betaH=betaH)
    }
  f.large <- formula(largeModel)
  attributes(f.large) <- NULL

  ans$f.large <- f.large
  ans$f.small <- f.small
  ans$ctime   <- (proc.time()-t0)[1]
  ans$L       <- L
  ans
}

.finalizeKR <- function(stats){

  test = list(
    Ftest      = c(stat=stats$Fstat,     ndf=stats$ndf,  ddf=stats$ddf,  F.scaling=stats$F.scaling,  p.value=stats$p.value),
    FtestU     = c(stat=stats$FstatU,    ndf=stats$ndf,  ddf=stats$ddf,  F.scaling=NA,               p.value=stats$p.valueU))
  test  <- as.data.frame(do.call(rbind, test))
  test$ndf <- as.integer(test$ndf)
  ans   <- list(test=test, type="F", aux=stats$aux, stats=stats)
  ## Notice: stats are carried to the output. They are used for get getKR function...
  class(ans)<-c("KRmodcomp")
  ans
}

KRmodcomp_internal <- function(largeModel, LL, betaH=0, details=0){

  PhiA  <- vcovAdj(largeModel, details)
  stats <- .KR_adjust(PhiA, Phi=vcov(largeModel), LL, beta=fixef(largeModel), betaH)
  stats <- lapply(stats, c) ## To get rid of all sorts of attributes
  ans   <- .finalizeKR(stats)
  ans
}


## --------------------------------------------------------------------
## This is the function that calculates the Kenward-Roger approximation
## --------------------------------------------------------------------
.KR_adjust <- function(PhiA, Phi, L, beta, betaH){
  Theta  <-  t(L) %*% solve( L %*% Phi %*% t(L), L)
  P <- attr( PhiA, "P" )
  W <- attr( PhiA, "W" )

  A1 <- A2 <- 0
  ThetaPhi <- Theta %*% Phi
  n.ggamma <- length(P)
  for (ii in 1:n.ggamma) {
    for (jj in c(ii:n.ggamma)) {
      e  <- ifelse(ii==jj, 1, 2)
      ui <- ThetaPhi %*% P[[ii]] %*% Phi
      uj <- ThetaPhi %*% P[[jj]] %*% Phi
      A1 <- A1 + e* W[ii,jj] * (.spur(ui) * .spur(uj))
      A2 <- A2 + e* W[ii,jj] * sum(ui * t(uj))
    }
  }

  q <- rankMatrix(L)
  B <- (1/(2*q)) * (A1+6*A2)
  g <- ( (q+1)*A1 - (q+4)*A2 )  / ((q+2)*A2)
  c1<- g/(3*q+ 2*(1-g))
  c2<- (q-g) / (3*q + 2*(1-g))
  c3<- (q+2-g) / ( 3*q+2*(1-g))
  ##  cat(sprintf("q=%i B=%f A1=%f A2=%f\n", q, B, A1, A2))
  ##  cat(sprintf("g=%f, c1=%f, c2=%f, c3=%f\n", g, c1, c2, c3))
###orgDef: E<-1/(1-A2/q)
###orgDef: V<- 2/q * (1+c1*B) /  ( (1-c2*B)^2 * (1-c3*B) )

  ##EE     <- 1/(1-A2/q)
  ##VV     <- (2/q) * (1+c1*B) /  ( (1-c2*B)^2 * (1-c3*B) )
  EE     <- 1 + (A2/q)
  VV     <- (2/q)*(1+B)
  EEstar <- 1/(1-A2/q)
  VVstar <- (2/q)*((1+c1*B)/((1-c2*B)^2 * (1-c3*B)))
  ##  cat(sprintf("EE=%f VV=%f EEstar=%f VVstar=%f\n", EE, VV, EEstar, VVstar))
  V0<-1+c1*B
  V1<-1-c2*B
  V2<-1-c3*B
  V0<-ifelse(abs(V0)<1e-10,0,V0)
  ##  cat(sprintf("V0=%f V1=%f V2=%f\n", V0, V1, V2))

###orgDef: V<- 2/q* V0 /(V1^2*V2)
###orgDef: rho <-  V/(2*E^2)

  rho <- 1/q * (.divZero(1-A2/q,V1))^2 * V0/V2
  df2 <- 4 + (q+2)/ (q*rho-1)          ## Here are the adjusted degrees of freedom.

###orgDef: F.scaling <-  df2 /(E*(df2-2))
###altCalc F.scaling<- df2 * .divZero(1-A2/q,df2-2,tol=1e-12)
  ## this does not work because df2-2 can be about 0.1
  F.scaling<-ifelse( abs(df2-2)<1e-2, 1 , df2*(1-A2/q)/(df2-2))
  ##cat(sprintf("KR: rho=%f, df2=%f F.scaling=%f\n", rho, df2, F.scaling))

  ## Vector of auxilary values; just for checking etc...
  aux <- c(A1=A1, A2=A2, V0=V0, V1=V1, V2=V2, rho=rho, F.scaling=F.scaling)

### The F-statistic; scaled and unscaled
  betaDiff <- cbind( beta - betaH )
  Wald     <- as.numeric(t(betaDiff) %*% t(L) %*% solve(L%*%PhiA%*%t(L), L%*%betaDiff))
  WaldU    <- as.numeric(t(betaDiff) %*% t(L) %*% solve(L%*%Phi%*%t(L), L%*%betaDiff))

  FstatU <- Wald/q
  pvalU  <- pf(FstatU,df1=q,df2=df2,lower.tail=FALSE)

  Fstat  <- F.scaling * FstatU
  pval   <- pf(Fstat,df1=q,df2=df2,lower.tail=FALSE)

  stats<-list(ndf=q, ddf=df2,
              Fstat  = Fstat,  p.value=pval, F.scaling=F.scaling,
              FstatU = FstatU, p.valueU = pvalU,
              aux = aux)
  stats

}















.KRcommon <- function(x){
  cat(sprintf("F-test with Kenward-Roger approximation; computing time: %.2f sec.\n",
              x$ctime))
  cat("large : ")
  print(x$f.large)

  if (inherits(x$f.small,"call")){
    cat("small : ")
    print(x$f.small)
  } else {
    formSmall <- x$f.small
    cat("small : L beta = L betaH \n")
    cat('L=\n')
    print(formSmall$L)
    cat('betaH=\n')
    print(formSmall$betaH)
  }

}


print.KRmodcomp <- function(x,...){

  .KRcommon(x)
  FF.thresh <- 0.2
  F.scale <- x$aux['F.scaling']
  tab <- x$test

  if (max(F.scale)>FF.thresh){
    printCoefmat(tab[1,,drop=FALSE], tst.ind=c(1,2,3), na.print='', has.Pvalue=TRUE)
  } else {
    printCoefmat(tab[2,,drop=FALSE], tst.ind=c(1,2,3), na.print='', has.Pvalue=TRUE)
  }
  return(invisible(x))
}


summary.KRmodcomp <- function(object,...){

  .KRcommon(object)
  FF.thresh <- 0.2
  F.scale <- object$aux['F.scaling']
  tab <- object$test

  printCoefmat(tab, tst.ind=c(1,2,3), na.print='', has.Pvalue=TRUE)

  if (F.scale<0.2 & F.scale>0) {
    cat('Note: The scaling factor for the F-statistic is smaller than 0.2 \n')
    cat('The Unscaled statistic might be more reliable \n ')
  } else {
    if (F.scale<=0){
      cat('Note: The scaling factor for the F-statistic is negative \n')
      cat('Use the Unscaled statistic instead. \n ')
    }
  }
}






















  #stats <- .KRmodcompPrimitive(largeModel, L, betaH, details)


## .KRmodcompPrimitive<-function(largeModel, L, betaH, details) {
##   PhiA<-vcovAdj(largeModel, details)
##   .KR_adjust(PhiA, Phi=vcov(largeModel), L, beta=fixef(largeModel), betaH )
## }


### SHD addition: calculate bartlett correction and gamma approximation
###
##   ## Bartlett correction - X2 distribution
##   BCval   <- 1 / EE
##   BCstat  <- BCval * Wald
##   p.BC    <- 1-pchisq(BCstat,df=q)
## #  cat(sprintf("Wald=%f BCval=%f BC.stat=%f p.BC=%f\n", Wald, BCval, BCstat, p.BC))
##   ## Gamma distribution
##   scale   <- q*VV/EE
##   shape   <- EE^2/VV
##   p.Ga    <- 1-pgamma(Wald, shape=shape, scale=scale)
## #  cat(sprintf("shape=%f scale=%f p.Ga=%f\n", shape, scale, p.Ga))
