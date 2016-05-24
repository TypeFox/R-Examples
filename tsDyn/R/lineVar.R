#'Multivariate linear models: VAR and VECM
#'
#'Estimate either a VAR or a VECM.
#'
#'This function provides basic functionalities for VAR and VECM models. More
#'comprehensive functions are in package \pkg{vars}. A few differences appear
#'in the VECM estimation: \describe{ 
#' \item{Engle-Granger estimator}{The Engle-Granger estimator is available}
#' \item{Presentation}{Results are printed in a different ways, using a matrix form}
#' \item{lateX export}{The matrix of coefficients can be exported to latex, 
#'        with or without standard-values and significance stars}
#' } 
#'
#'Two estimators are available: the Engle-Granger two step
#'approach (\code{2OLS}) or the Johansen (\code{ML}). For the 2OLS,
#'deterministic regressors (or external variables if \code{LRinclude} is of
#'class numeric) can be added for the estimation of the cointegrating value and
#'for the ECT. This is only working when the beta value is not pre-specified.
#'
#' The argument \code{beta} is only for \code{\link{VECM}}, look at the specific help page for more details. 
#'
#'@param data multivariate time series (first row being first=oldest value)
#'@param lag Number of lags to include in each regime
#'@param r Number of cointegrating relationships
#'@param include Type of deterministic regressors to include
#'@param model Model to estimate. Either a VAR or a VECM
#'@param I For VAR only: whether in the VAR the variables are to be taken in
#'levels (original series) or in difference, or similarly to the univariate ADF
#'case.
#'@param beta for VECM only: imposed cointegrating value. If null, will be estimated
#'@param LRinclude Possibility to include in the long-run relationship and the
#'ECT trend, constant... Can also be a matrix with exogeneous regressors
#'@param estim Type of estimator for the VECM: '2OLS' for the two-step approach
#'or 'ML' for Johansen MLE
#'@param exogen Inclusion of exogenous variables (first row being first=oldest
#'value). Is either of same size than data (then automatically cut) or than
#'end-sample.
#'@return Fitted model data
#'@author Matthieu Stigler
#'@seealso \code{\link{VECM}} which is just a wrapper for
#'\code{lineVar(...,model="VECM")}
#'
#'\code{\link{TVAR}} and \code{\link{TVECM}} for the corresponding threshold
#'models. \code{\link{linear}} for the univariate AR model.
#'@keywords ts
#'@export
#'@examples
#'
#'data(zeroyld)
#'data<-zeroyld
#'
#'#Fit a VAR
#'VAR<-lineVar(data, lag=1)
#'VAR
#'summary(VAR)
#'
#'#compare results with package vars:
#'if(require(vars)) {
#'	a<-VAR(data, p=1)
#'	coef_vars <- t(sapply(coef(a), function(x) x[c(3,1,2),1]))
#'	all.equal(coef(VAR),coef_vars, check.attributes=FALSE)
#'}
#'
#'###VECM
#'VECM.EG<-lineVar(data, lag=2, model="VECM")
#'VECM.EG
#'summary(VECM.EG)
#'
#'VECM.ML<-lineVar(data, lag=2, model="VECM", estim="ML")
#'VECM.ML
#'summary(VECM.ML)
#'
#'
#'###Check Johansen MLE
#'myVECM<-lineVar(data, lag=1, include="const", model="VECM", estim="ML")
#'summary(myVECM, digits=7) 
#'#comparing with vars package
#'if(require(vars)){
#'	a<-ca.jo(data, spec="trans")
#'	summary(a)
#'	#same answer also!
#'}
#'
#'##export to Latex
#'toLatex(VECM.EG)
#'toLatex(summary(VECM.EG))
#'options("show.signif.stars"=FALSE)
#'toLatex(summary(VECM.EG), parenthese="Pvalue")
#'options("show.signif.stars"=TRUE)
#'
#'
#'
lineVar<-function(data, lag, r=1,include = c( "const", "trend","none", "both"), model=c("VAR", "VECM"), 
		  I=c("level", "diff", "ADF"),beta=NULL, estim=c("2OLS", "ML"),
		  LRinclude=c("none", "const", "trend","both"), exogen=NULL){

  y <- as.matrix(data)
  Torigin <- nrow(y) 	#Size of original sample
  T <- nrow(y) 		#Size of start sample

  if(length(lag)==1){
    p <- lag
    notAllLags<-FALSE
    Lags<-1:p
  } else {
    notAllLags<-TRUE
    p<-max(lag)
    Lags<-lag
  }

  t <- T-p 		#Size of end sample
  k <- ncol(y) 		#Number of variables
  t<-T-p			#Size of end sample

  if(is.null(colnames(data)))
	  colnames(data)<-paste("Var", c(1:k), sep="")

###Check args
  include<-match.arg(include)
  LRinclude<-match.arg(LRinclude)
  if(lag==0) {
    warning("Lag=0 not fully implemented, methods not expected to work: fevd, predict, irf,...")
  }
  if(LRinclude%in%c("const", "both"))  include<-"none"
  ninclude<-switch(include, "const"=1, "trend"=1,"none"=0, "both"=2)
  model<-match.arg(model)
  estim<-match.arg(estim)
  I<-match.arg(I)
  if(!r%in%1:(k-1)) stop("Arg r, the number of cointegrating relationships, should be between 1 and K-1\n")
  if(model=="VECM"&estim=="2OLS"&r>1){
    warning("Estimation of more than 1 coint relationship is not possible with estim '2OLS'. Switched to Johansen 'ML'\n")
    estim<-"ML"
  }

  minPara<-p*k+ninclude
  if(!t>minPara) stop("Not enough observations. Try reducing lag number\n")

###Construct variables
  Y <- y[(p+1):T,] #
  X <- embed(y, p+1)[, -seq_len(k)]	#Lags matrix

  #Set up of dependant and independant variables matrices
  if(notAllLags)
    X<-X[,sort(rep((Lags*k-k+1), k))+0:(k-1)]

  DeltaY<-diff(y)[(p+1):(T-1),]
  Xminus1<-embed(y,p+2)[,(k+1):(k+k)]
  DeltaX<-embed(diff(y),p+1)[,-(1:k)]

  if(model=="VAR"){
    if(I=="level"){
      Z<-X
      Y<-Y
    } else if(I=="diff"){
      Z<-DeltaX
      Y<-DeltaY
      t<-t-1
    } else if(I=="ADF"){
      Z<-cbind(Xminus1, DeltaX)
      Y<-DeltaY
      t<-t-1
    }
  } else if(model=="VECM"){
    Z<-DeltaX
    Y<-DeltaY
    t<-t-1
  }

###Regressors matrix
  incReg <- switch(include,
                   "const" = matrix(1, nrow=t, ncol=1), 
                   "trend" = matrix(seq_len(t), ncol=1),
                   "both"  = cbind(rep(1,t),seq_len(t)),
                   "none"=NULL)
  Z <-if(lag==0) incReg else cbind(incReg, Z)

  if(!is.null(exogen)){
    if(is.data.frame(exogen)|is.vector(exogen)) exogen <- as.matrix(exogen)
    n_exo <- NROW(exogen)
    if(n_exo!=t){
      if(n_exo!=T)  warning("exogen is of size ", n_exo, "while full/end-sample size is of size", T,"/", nrow(Z), "series shortened")
      exogen <- myTail(exogen, t, addrownums=FALSE)
    }
    Z <- if(lag==0 & include=="none") exogen else cbind(Z, exogen)
  }


##VECM: Long-run relationship OLS estimation
  if(model=="VECM"&estim=="2OLS"){
	  #beta has to be estimated
    beta.estimated <- is.null(beta)
    if(is.null(beta) ){

    ## build LRplus: deterministic/exogeneous regressor in coint
      if(class(LRinclude)=="character"){
        LRplus <-switch(LRinclude, "none"=NULL,"const"=rep(1,T),"trend"=seq_len(T),"both"=cbind(rep(1,T),seq_len(T)))
        LRinc_name <- switch(LRinclude, "const"="const", "trend"="trend", "both"=c("const", "trend"), "none"=NULL)
        LRinc_dim <- switch(LRinclude, "const"=1, "trend"=1, "both"=2, "none"=0)
      } else if(class(LRinclude)%in%c("matrix", "numeric")) {
        LRplus<-LRinclude
      } else{
        stop("Argument LRinclude badly given")
      }
    ## run coint regression
      if(LRinclude=="none"){
        cointLM<-lm(y[,1] ~  y[,-1]-1)
      } else {
        cointLM<-lm(y[,1] ~  y[,-1]-1+ LRplus)
        Xminus1 <- cbind(Xminus1, tail(LRplus,nrow(Xminus1)))
      }
      
      betaLT<-coint<-c(1,-cointLM$coef)
      betaLT_std <- c(1,summary(cointLM)$coef[,2])
      names(betaLT_std)<-c(colnames(data), LRinc_name)

## case beta pre-estimated
    } else {
      if(length(beta)!=k-1) stop("Arg 'beta' should be of length k-1")
      if(LRinclude!="none")
        warning("Arg LRinclude not taken into account when beta is given by user")
      LRinc_name <- NULL
      LRinc_dim <- 0
      coint<-c(1, -beta)
      betaLT<-c(1,-beta)
    }

    coint_export<-matrix(coint, nrow=k+LRinc_dim , dimnames=list(c(colnames(data),LRinc_name), "r1"))
    betaLT<-matrix(betaLT, nrow=k+LRinc_dim , dimnames=list(c(colnames(data),LRinc_name),"r1"))
    ECTminus1<-Xminus1%*%betaLT
    Z<-cbind(ECTminus1,Z)
  }

##VECM: ML (Johansen ) estimation of cointegrating vector
  else if(model=="VECM"&estim=="ML"){
    beta.estimated<- is.null(beta)
    if(lag==0 & include=="none" & is.null(exogen)){
     u <- Y
     v <- Xminus1
    } else {
      #Auxiliary regression 1
      reg_res1<-lm.fit(Z,Y)
      u<-residuals(reg_res1)
      #Auxiliary regression 2
      reg_res2<-lm.fit(Z,Xminus1)
      v<-residuals(reg_res2)
    }
    #Auxiliary regression 3
    if(LRinclude!="none"){
      add <- switch(LRinclude, "const"=matrix(1, nrow=t), 
                    "trend"=matrix(1:t, nrow=t), 
                    "both"=cbind(1,1:t))
      reg_res3 <- if(lag==0) add else residuals(lm.fit(Z,add))
      v<-cbind(v,reg_res3) # equ 20.2.46 in Hamilton 
    }
    #Moment matrices
    S00<-crossprod(u)
    S11<-crossprod(v)
    S01<-crossprod(u,v)
    if(beta.estimated){
      SSSS<-solve(S11)%*%t(S01)%*%solve(S00)%*%S01
      eig<-eigen(SSSS)
      ve<-Re(eig$vectors)
      va<-Re(eig$values)
      #normalize eigenvectors
      ve_no<-apply(ve,2, function(x) x/sqrt(t(x)%*%S11%*%x))
      ve_2<-t(t(ve_no)/diag(ve_no)) 
      ve_3<-ve_2[,1:r, drop=FALSE]
      C2 <- matrix(0, nrow = nrow(ve_2) - r, ncol = r)
      C <- rbind(diag(r), C2)
      ve_4 <- ve_3 %*% solve(t(C) %*% ve_3)

      #compute A (speed adjustment)
      z0<-t(u)%*%v%*%tcrossprod(ve_no[,1:r])

	###Slope parameters
      if(LRinclude!="none"){
        ECTminus1<-cbind(Xminus1,add)%*%ve_4
      }else{
        ECTminus1<-Xminus1%*%ve_4
      }
      coin_ve_names <- switch(LRinclude, "const"="const", "trend"="trend", "both"=c("const", "trend"), "none"=NULL)
      dimnames(ve_4)<-list(c(colnames(data), coin_ve_names), paste("r", 1:r, sep=""))
      betaLT<-ve_4
#     if beta restricted:
    } else {
      betaLT <- as.matrix(beta)
      if(ncol(betaLT)!=r) stop("Argument 'beta' should have as many columns as 'r'")
      if(ncol(betaLT)==1&&nrow(betaLT)==k-1) betaLT <- rbind(1,betaLT)
      if(nrow(betaLT)!=k) stop("Argument 'beta' should have as many rows as 'k'")
      ECTminus1 <- Xminus1%*%betaLT 
      ## restricted ML:
      S11_r <- t(betaLT)%*%S11%*%betaLT ## equa 20.3.11 Hamilton 1994, p. 649
      S01_r <- S01 %*%betaLT ## equa 20.3.12 Hamilton 1994, p. 649
      SSSS_r <-solve(S11_r)%*%t(S01_r)%*%solve(S00)%*%S01_r
      eig_r<-eigen(SSSS_r)
      va<-Re(eig_r$values)
     }
    Z<-cbind(ECTminus1,Z)
  }#end model=="VECM"&estim=="ML"


###Slope parameters, residuals and fitted
#   B<-t(Y)%*%Z%*%solve(t(Z)%*%Z)		#B: OLS parameters, dim 2 x npar
  B<- t(qr.coef(qr(Z),Y))
  fitted<-Z%*%t(B)
  res<-Y-fitted

###naming of variables and parameters
  npar<-ncol(B)*nrow(B)
  rownames(B)<-paste("Equation",colnames(data))
  if(p>0){
    LagNames<-c(paste(rep(colnames(data),length(Lags)), -rep(Lags, each=k)))
    if(I=="ADF") LagNames <- paste("D", LagNames,sep="_")
  } else {
    LagNames <- NULL
  }
  ECT<- if(model=="VECM") paste("ECT", if(r>1) 1:r else NULL, sep="") else NULL
  Xminus1Names<- if(I=="ADF") paste(colnames(data),"-1",sep="") else NULL
  BnamesInter<-switch(include,"const"="Intercept","none"=NULL,"trend"="Trend","both"=c("Intercept","Trend"))
  if(!is.null(exogen)){
    exo_names <- colnames(as.matrix(exogen))
    exp_res <- "-[0-9]+|\\.l[0-9]+|ECT|Intercept|Trend"
    if(any(grepl(exp_res, exo_names))) {
      warning("Exogen contains reserved names (", exp_res, ". Changed to exo_x")
      exo_names[grepl(exp_res, exo_names)] <- paste("exo", 1:sum(grepl(exp_res, exo_names)), sep="_")
    }
    if(any(is.null(exo_names))) exo_names <- paste("exo", 1:NCOL(exogen), sep="_")
  } else {
    exo_names <- NULL
  }
  Bnames<-c(ECT,BnamesInter,Xminus1Names, LagNames, exo_names)
  colnames(B)<-Bnames

###Y and regressors matrix to be returned
  naX<-rbind(matrix(NA, ncol=ncol(Z), nrow=T-t), Z)
  rownames(naX)<-rownames(data)
  YnaX<-cbind(data, naX)
  colnames(YnaX)<-c(colnames(data),Bnames)

###Return outputs
  model.specific<-list()
  model.specific$nthresh<-0

  if(model=="VECM"){
    model.specific$beta<- betaLT
    model.specific$coint <- betaLT
    model.specific$r<-r
    model.specific$estim<-estim
    model.specific$LRinclude<-LRinclude
    model.specific$beta.estimated<-beta.estimated
    model.specific$estim <- estim
    
    if(estim=="ML"){
      model.specific$S00<-S00
      model.specific$lambda<-va
    } 
  }


  z<-list(residuals=res,  
          coefficients=B,  k=k, t=t,T=T, npar=npar, nparB=ncol(B), type="linear", 
          fitted.values=fitted, 
          model.x=Z, 
          include=include,
          lag=lag, 
          model=YnaX, 
          df.residual=t-npar/k, 
          exogen = !is.null(exogen),
          num_exogen = if(!is.null(exogen)) NCOL(exogen) else 0,
          model.specific=model.specific)
  if(model=="VAR"){
    class(z)<-c("VAR","nlVar")
  } else {
    class(z)<-c("VECM","VAR", "nlVar")
    I<-"diff"
  }

  attr(z, "varsLevel")<-I
  attr(z, "model")<-model
  return(z)
}


#### VECM function: wrapper to lineVar


#' Estimation of Vector error correction model (VECM)
#'
#' Estimate either a VECM by Engle-Granger or Johansen (MLE) method.
#'
#' This function is just a wrapper for the \code{\link{lineVar}}, with
#' model="VECM".
#'
#' More comprehensive functions for VECM are in package \pkg{vars}. A few
#' differences appear in the VECM estimation: \describe{ \item{Engle-Granger
#' estimator}{The Engle-Granger estimator is available}
#' \item{Presentation}{Results are printed in a different ways, using a matrix
#' form} \item{lateX export}{The matrix of coefficients can be exported to
#' latex, with or without standard-values and significance stars}
#' \item{Prediction}{The \code{predict} method contains a \code{newdata}
#' argument allowing to compute rolling forecasts.} }
#'
#' Two estimators are available: the Engle-Granger two step approach
#' (\code{2OLS}) or the Johansen (\code{ML}). For the 2OLS, deterministics
#' regressors (or external variables if LRinclude is of class numeric) can be
#' added for the estimation of the cointegrating value and for the ECT. This is
#' only working when the beta value is not pre-specified.
#'
#' The arg beta is the cointegrating value, the cointegrating vector will be
#' taken as: (1, -beta).
#'
#' Note that the lag specification corresponds to the lags in the VECM
#' representation, not in the VAR (as is done in package vars or software
#' GRETL). Basically, a VAR with 2 lags corresponds here to a VECM with 1 lag.
#' Lag 0 in the VECM is not allowed.
#' 
#' #'The arg \code{beta} allows to specify constrained cointegrating values, leading to
#' \eqn{ECT= \beta^{'}X_{t-1}}. It should be specified as a \eqn{K \times r} matrix. In case of
#' \eqn{r=1}, can also be specified as a vector. Note that the vector should be normalised, 
#' with the first value to 1, and the next values showing the opposite sign in the long-run relationship \eqn{- \beta}. 
#' In case the vector has \eqn{K-1} values, this is what \code{lineVar} is doing, setting \eqn{(1, - \beta)}. 
#' Note finally one should provide values for all
#' the coefficients (eventually except for special case of r=1 and k-1), if you want to provide only part of the 
#' parameters, and let the others be estimated, look at the functions in package urca. 
#'
#' @param data multivariate time series (first row being first=oldest value)
#' @param lag Number of lags (in the VECM representation, see Details)
#' @param r Number of cointegrating relationships
#' @param include Type of deterministic regressors to include
#' @param beta for VECM only: imposed cointegrating value. If null, will be estimated
#'so values will be estimated
#'@param LRinclude Type of deterministic regressors to include in the long-term
#'relationship. Can also be a matrix with exogeneous regressors (2OLS only).
#'@param estim Type of estimator: \code{2OLS} for the two-step approach or
#'\code{ML} for Johansen MLE
#'@param exogen Inclusion of exogenous variables (first row being first=oldest
#'value). Is either of same size than data (then automatically cut) or than
#'end-sample.
#'@return An object of class \code{VECM} (and higher classes \code{VAR} and
#'\code{nlVar}) with methods: \describe{ \item{Usual methods}{Print, summary,
#'plot, residuals, fitted, vcov} \item{Fit criteria}{AIC, BIC,
#'\code{\link{MAPE}}, \code{\link{mse}}, \code{\link[=logLik.VECM]{logLik}}
#'(latter only for models estimated with MLE)} \item{Prediction}{Predict and
#'\code{\link{predict_rolling}}} \item{VAR/VECM methods}{Impulse response
#'function (\code{\link[=irf.nlVar]{irf}}) and forecast error variance
#'decomposition (\code{\link[=fevd.nlVar]{fevd}})} \item{LaTeX}{toLatex} }
#'@author Matthieu Stigler
#'@seealso \code{\link{lineVar}} \code{\link{TVAR}} and \code{\link{TVECM}} for
#'the correspoding threshold models. \code{\link{linear}} for the univariate AR
#'model.
#'@keywords ts
#'@export
#'@examples
#'
#'data(zeroyld)
#'data<-zeroyld
#'
#'#Fit a VECM with Engle-Granger 2OLS estimator:
#'vecm.eg<-VECM(zeroyld, lag=2)
#'
#'#Fit a VECM with Johansen MLE estimator:
#'vecm.jo<-VECM(zeroyld, lag=2, estim="ML")
#'
#'#compare results with package vars:
#'if(require(vars)) {
#'  data(finland)
#'  #check long coint values
#'    all.equal(VECM(finland, lag=2, estim="ML", r=2)$model.specific$beta, 
#'              cajorls(ca.jo(finland, K=3, spec="transitory"), r=2)  $beta, check.attributes=FALSE)
#' # check OLS parameters
#'   all.equal(t(coefficients(VECM(finland, lag=2, estim="ML", r=2))), 
#'     coefficients(cajorls(ca.jo(finland, K=3, spec="transitory"), r=2)$rlm), check.attributes=FALSE)
#'
#'}
#'
#'
#'##export to Latex
#'toLatex(vecm.eg)
#'toLatex(summary(vecm.eg))
#'options("show.signif.stars"=FALSE)
#'toLatex(summary(vecm.eg), parenthese="Pvalue")
#'options("show.signif.stars"=TRUE)
#'
#'
#'
VECM<-function(data, lag,r=1, include = c( "const", "trend","none", "both"), beta=NULL, estim=c("2OLS", "ML"),LRinclude=c("none", "const", "trend","both"), exogen=NULL)
  lineVar(data, lag, r=r,include = include, model="VECM" ,beta=beta, estim=estim,LRinclude=LRinclude,exogen=exogen)



###Testing
if(FALSE) { #usage example
###Hansen Seo data
library(tsDyn)
environment(lineVar)<-environment(star)
environment(summary.VAR)<-environment(star)
environment(toLatex.VAR)<-environment(star)
#data(zeroyld)
dat<-zeroyld

#tests
aVAR<-lineVar(dat[1:100,], lag=c(1,2), include="both", model="VAR")
#lag2, 2 thresh, trim00.05: 561.46
class(aVAR)
aVAR
print(aVAR)
logLik(aVAR)
AIC(aVAR)
BIC(aVAR)
deviance(aVAR)
coef(aVAR)
summary(aVAR)
toLatex(aVAR)
toLatex(summary(aVAR))
}



#' @S3method print VAR
print.VAR <- function(x,...){
	print(coef(x))
}

#' @S3method summary VAR
summary.VAR<-function(object, digits=4,...){
  x<-object
  r<-4
  t<-x$t
  k<-x$k
  Sigma<-matrix(1/(object$df.residual)*crossprod(x$residuals),ncol=k)
  XX <- crossprod(x$model.x)
  cov.unscaled <- try(solve(XX), silent=TRUE)
  if(inherits(cov.unscaled, "try-error")) {
    Qr <- qr(x$model.x)
    p1 <- 1:Qr$rank
    cov.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
    warning("Potential numerical unstability, beware of standard errors\n")
  }
  VarCovB<-cov.unscaled%x%Sigma
  StDevB<-matrix(sqrt(diag(VarCovB)), nrow=k)
  Tvalue<-x$coefficients/StDevB
  
  Pval<-pt(abs(Tvalue), df=(nrow(x$model.x)-ncol(x$model.x)), lower.tail=FALSE)+pt(-abs(Tvalue), df=(nrow(x$model.x)-ncol(x$model.x)), lower.tail=TRUE)
	#Pval<-round(Pval,4)
  symp <- symnum(Pval, corr=FALSE,cutpoints = c(0,  .001,.01,.05, .1, 1), symbols = c("***","** ","*  ",".  ","    "))
  stars<-matrix(symp, nrow=nrow(Pval))
  ab<-matrix(paste(myformat(x$coefficients,digits),"(", myformat(StDevB,digits),")",stars,sep=""), nrow=nrow(Pval))
  dimnames(ab)<-dimnames(x$coefficients)		

  x$bigcoefficients<-ab
  x$cov.unscaled<-cov.unscaled
  x$sigma<-Sigma
  x$StDev<-StDevB
  x$Pvalues<-Pval
  x$stars<-stars
  x$starslegend<-symp
  x$aic<-AIC.nlVar(x)
  x$bic<-BIC.nlVar(x)
  x$SSR<-deviance.nlVar(x)
  class(x)<-c("summary.VAR", "VAR")
  return(x)
}


#' @S3method print summary.VAR
print.summary.VAR<-function(x,...){
  cat("#############\n###Model", attr(x,"model"),"\n#############")
  cat("\nFull sample size:",x$T, "\tEnd sample size:", x$t) 
  cat("\nNumber of variables:", x$k,"\tNumber of estimated slope parameters", x$npar)
  cat("\nAIC",x$aic , "\tBIC", x$bic, "\tSSR", x$SSR)
  if(attr(x,"model")=="VECM"){
    cat("\nCointegrating vector ", if(x$model.specific$beta.estimated) paste("(estimated by ", x$model.specific$estim, "):\n", sep="") else "(fixed by user):\n", sep="")
    print(t(x$model.specific$beta))
  }
  cat("\n\n")
  print(noquote(x$bigcoefficients))

}

#' @S3method vcov VAR
vcov.VAR<-function(object, ...){
  sum<-summary.VAR(object)
  so<-sum$cov.unscaled%x%sum$sigma
  co.names<-gsub(" ", "", colnames(coef(object)))
  eq.names<-gsub("Equation ", "",rownames(coef(object)))
  together.names<-paste(rep(eq.names,each= length(co.names)), co.names, sep=":")
  dimnames(so)<-list(together.names, together.names)
  so
}

#' @S3method toLatex VAR
toLatex.VAR<-function(object,..., digits=4, parenthese=c("StDev","Pvalue"), label){
  x<-object
  if(attr(x,"model")=="VECM"&&x$model.specific$LRinclude!="none") stop("toLatex not implemented now for models with arg 'LRinclude' different from 'none'") 
  parenthese<-match.arg(parenthese)
  if(inherits(x,"summary.VAR")){
    a<-myformat(x$coefficients,digits, toLatex=TRUE)
    inp<-switch(parenthese, "StDev"= x$StDev, "Pvalue"= x$Pvalues )
    b<-myformat(inp ,digits,toLatex=TRUE)
    if(getOption("show.signif.stars"))    
      stars<-paste("^{",x$stars,"}", sep="")
    else
      stars<-NULL
    coeftoprint<-matrix(paste(a,"(",b,")",stars, sep=""),ncol=ncol(a), nrow=nrow(a))
  }#end if x is of class summary
  
  else{
    coeftoprint <-myformat(x$coefficients, digits, toLatex=TRUE)}
  varNames<-rownames(x$coefficients)
  res<-character()
  res[1]<-"%insert in the preamble and uncomment the line you want for usual /medium /small matrix"
  res[2]<-"%\\usepackage{amsmath} \\newenvironment{smatrix}{\\begin{pmatrix}}{\\end{pmatrix}} %USUAL"
  res[3]<-"%\\usepackage{amsmath} \\newenvironment{smatrix}{\\left(\\begin{smallmatrix}}{\\end{smallmatrix}\\right)} %SMALL"
  res[4]<-"%\\usepackage{nccmath} \\newenvironment{smatrix}{\\left(\\begin{mmatrix}}{\\end{mmatrix}\\right)} %MEDIUM"
  res[5]<-"\\begin{equation}"
  if(!missing(label)) res[5]<- paste(res[5], "\\label{", label, "}", sep="")
  res[6]<- "\\begin{smatrix} %explained vector"

###explained vector
  if(attr(x, "varsLevel")=="diff")
    res[7]<-TeXVec(paste("slashDelta X_{t}^{",seq(1, x$k),"}", sep=""))
  else
    res[7]<-TeXVec(paste("X_{t}^{",seq(1, x$k),"}", sep=""))
  res[8]<- "\\end{smatrix}="
###ECT
  ninc<-switch(x$include, "const"=1, "trend"=1,"none"=0, "both"=2)
  if(attr(x,"model")=="VECM"){
    r<-x$model.specific$r
    len<-length(res)
    if(r==1){
      res[len+1]<-"+\\begin{smatrix}  %ECT"
      res[len+2]<-TeXVec(coeftoprint[,1]) #see nlVar-methods.R
      res[len+3]<-"\\end{smatrix}ECT_{-1}"
    }else{
      res[len+1]<-"+\\begin{smatrix}  %ECT"
      res[(len+2):(len+x$k+1)]<-TeXMat(coeftoprint[,1:r]) #see nlVar-methods.R
      res[len+x$k+2]<-"\\end{smatrix}ECT_{-1}"
    }
  }
###Const
  a<-if(attr(x,"model")=="VECM") r else 0
  res<-include(x, res, coeftoprint, skip=a)	#nlVar-methods.R
###Lags
  res<-LagTeX(res, x, coeftoprint, ninc+a)	#nlVar-methods.R
  res[length(res)+1]<-"\\end{equation}"
  res<-gsub("slash", "\\", res, fixed=TRUE)
  res<-res[res!="blank"]
  
  return(structure(res, class="Latex"))
}


if(FALSE){
###TODO
#check if const/trend/both in LR rel and VECM makes sense!
#check for standaard deviation of coint vector whith ML estim!
#consistency between ML and OLS coint estimator?
}


if(FALSE) { #usage example
###Hansen Seo data
library(tsDyn)
#data(zeroyld)
dat<-zeroyld
environment(lineVar)<-environment(star)
environment(summary.VAR)<-environment(star)

aVAR<-lineVar(dat, lag=1, include="both", model="VAR")
aVAR<-lineVar(dat, lag=1, include="const", model="VECM", estim="ML", beta=0.98)
#lag2, 2 thresh, trim00.05: 561.46
aVAR
summary(aVAR)
sqrt(diag(summary(aVAR, cov=0)$sigma))
vcov.VAR(aVAR)
vcovHC.VAR(aVAR)
logLik(aVAR)
AIC(aVAR)
BIC(aVAR)
deviance(aVAR)
coef(aVAR)
environment(toLatex.VAR)<-environment(star)
toLatex(aVAR)
toLatex(summary(aVAR))

###Check VAR: comparing with vars
myVAR<-lineVar(dat, lag=1)

library(vars)
var<-VAR(dat, lag=1)

vaco1<-coef(var)$short.run[c(3,1,2),1]
vaco2<-coef(var)$long.run[c(3,1,2),1]
round(coef(myVAR),8)==round(rbind(vaco1, vaco2),8)

###Check Johansen MLE
myVECM<-lineVar(dat, lag=1, include="const", model="VECM", estim="ML")
summary(myVECM, digits=7) 
#comparing with Hansen paper:reported in Gauss procedure is:
#coint vector: 1.02206: ok!
#coeff: 
#comparing with vars package
a<-ca.jo(dat, spec="trans")
summary(a)
#same answer also!

set.seed(123)
rn <- rnorm(2*nrow(dat))
ex_1_n <- ex_1 <- rn[1:nrow(dat)]
ex_2_n <- ex_2 <- matrix(rn, ncol=2)
names(ex_1_n) <- "exoVar"
colnames(ex_2_n) <- c("exoVar1", "exoVar2")

lineVar(dat, lag=1, include="both", model="VAR", exogen=ex_1)
lineVar(dat, lag=1, include="both", model="VAR", exogen=ex_2)
lineVar(dat, lag=1, include="both", model="VAR", exogen=ex_1_n)
lineVar(dat, lag=1, include="both", model="VAR", exogen=ex_2_n)

colnames(ex_2_n) <- c("exoVar1", "a -1")
lineVar(dat, lag=1, include="both", model="VAR", exogen=ex_2_n)

colnames(ex_2_n) <- c("a.l2", "Trend")
lineVar(dat, lag=1, include="both", model="VAR", exogen=ex_2_n)

}
