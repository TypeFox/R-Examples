#'Selection of the cointegrating rank with Information criterion.
#'
#'Selection of the cointegrating rank and the lags with Information criterion
#'(AIC, BIC).
#'
#'This function estimates the AIC, BIC and Hannan-Quinn for each rank (up to
#'\code{lags.max}) and lags (up to \code{lags.max}). This method has been shown
#'to be useful to select simultaneously the rank and the lags, see references.
#'
#'@aliases rank.select print.rank.select summary.rank.select
#'@param data multivariate time series.
#'@param lag.max Maximum number of lags to investigate.
#'@param r.max Maximum rank to investigate.
#'@param include Type of deterministic regressors to innclude. See
#'\code{\link{VECM}} or \code{\link{lineVar}}.
#'@param fitMeasure Whether the AIC/BIC should be based on the full likelihood,
#'or just the SSR. See explanations in \code{\link{logLik.VECM}}.
#'@param sameSample Logical. Whether the data should be shortened so that the
#'AIC/BIC are estimated on the same sample. Default to TRUE.
#'@param returnModels Logical, default to FALSE. Whether the output should also
#'contain the list of each model computed.
#'@param x The output from \code{rank.select} for the print method.
#'@param object The output from \code{rank.select} for the summary method.
#'@param ... Unused.
#'@return An object of class \sQuote{rank.select}, with \sQuote{print} and
#'\sQuote{summary methods}, containing among other the matrices of AIC/BIC/HQ,
#'the Likelihood, and best ranks according to each criterion.
#'@author Matthieu Stigler
#'@export
#'@seealso \code{\link{VECM}} for estimating a VECM. \code{\link{rank.test}}
#'(or \code{\link[urca]{ca.jo}} in package \pkg{urca}) for the classical
#'Johansen cointegration test.
#'@references - Aznar A and Salvador M (2002). Selecting The Rank Of The
#'Cointegration Space And The Form Of The Intercept Using An Information
#'Criterion. Econometric Theory, *18*(04), pp. 926-947. <URL:
#'http://ideas.repec.org/a/cup/etheor/v18y2002i04p926-947_18.html>.
#'
#'-Cheng X and Phillips PCB (2009). Semiparametric cointegrating rank
#'selection. Econometrics Journal , *12*(s1), pp. S83-S104. <URL:
#'http://ideas.repec.org/a/ect/emjrnl/v12y2009is1ps83-s104.html>.
#'
#'- Gonzalo J and Pitarakis J (1998). Specification via model selection in
#'vector error correction models. Economics Letters, *60*(3), pp. 321 - 328.
#'ISSN 0165-1765, <URL: http://dx.doi.org/DOI: 10.1016/S0165-1765(98)00129-3>.
#'
#'- Kapetanios G (2004). The Asymptotic Distribution Of The Cointegration Rank
#'Estimator Under The Akaike Information Criterion. Econometric Theory,
#'*20*(04), pp. 735-742. <URL:
#'http://ideas.repec.org/a/cup/etheor/v20y2004i04p735-742_20.html>.
#'
#'- Wang Z and Bessler DA (2005). A Monte Carlo Study On The Selection Of
#'Cointegrating Rank Using Information Criteria. Econometric Theory, *21*(03),
#'pp. 593-620. <URL:
#'http://ideas.repec.org/a/cup/etheor/v21y2005i03p593-620_05.html>.
#'@keywords ts VECM cointegration
#'@examples
#'
#'
#'data(barry)
#'
#'# 
#'rk_sel <- rank.select(barry)
#'rk_sel
#'summary(rk_sel)
#'
#'
rank.select <- function(data, lag.max=10, r.max=ncol(data)-1, include = c( "const", "trend","none", "both"), fitMeasure=c("SSR", "LL"), sameSample=TRUE, returnModels=FALSE) {

  fitMeasure <- match.arg(fitMeasure)
  include <- match.arg(include)

  models_list <- list()

  VARtype <- if(r.max==0) "level" else "diff"

## Loop to get each model:
  for(i in 1:lag.max){
    models_list[[i]] <- list()
    if(sameSample){
      data_cut <- if(i==lag.max) data else data[-c(1:(lag.max-i)),]
    } else {
      data_cut <- data
    }
  ## VAR: rank 0 (on diffs)
    models_list[[i]][[1]] <- try(lineVar(data_cut, lag=i, I=VARtype, include=include ), silent=TRUE)
    if(inherits(models_list[[i]][[1]], "try-error")) models_list[[i]][[1]] <- NA

  ## VECM: rank 1 to k-1
    if(r.max>0){
      for(j in 1:r.max){
	models_list[[i]][[j+1]] <- try(VECM(data_cut, lag=i, r=j, estim="ML", include=include), silent=TRUE)
	if(inherits(models_list[[i]][[j+1]], "try-error")) models_list[[i]][[j+1]] <- NA
      }
      ## VAR: full rank
      models_list[[i]][[r.max+2]] <- try(lineVar(data_cut, lag=i+1, I="level", include=include), silent=TRUE)
    }
  }

  logLik.logical <- function(x) NA
  if_nona <- function(x, fun,...) if(inherits(x, "VAR")) fun(x,...) else NA

  AICs <- matrix(sapply(models_list, function(x) sapply(x,function(x) if(inherits(x, "VAR")) AIC(x,fitMeasure=fitMeasure) else NA)), ncol=lag.max)
  AICs <- matrix(sapply(models_list, function(x) sapply(x,function(x) if_nona(x,AIC, fitMeasure=fitMeasure))), ncol=lag.max)
  HQs <- matrix(sapply(models_list, function(x) sapply(x,function(x) if_nona(x,AIC, fitMeasure=fitMeasure,k=2*log(log(x$t))))), ncol=lag.max)
  LLs <- matrix(sapply(models_list, function(x) sapply(x,function(x) if_nona(x,logLik))), ncol=lag.max)
  BICs <- matrix(sapply(models_list, function(x) sapply(x,function(x) if_nona(x,BIC, fitMeasure=fitMeasure))), ncol=lag.max)

  names_mats <- list(paste("r", 0:(r.max+min(r.max,1)), sep="="), paste("lag", 1:lag.max, sep="="))
  dimnames(AICs) <- dimnames(BICs) <- dimnames(LLs) <- dimnames(HQs)<- names_mats

## Best values:
  AIC_min <- which(AICs==min(AICs, na.rm=TRUE),arr.ind=TRUE)
  BIC_min <- which(BICs==min(BICs, na.rm=TRUE),arr.ind=TRUE)
  HQ_min <- which(HQs==min(HQs, na.rm=TRUE),arr.ind=TRUE)
  AIC_min[1,1] <- AIC_min[1,1]-1
  BIC_min[1,1] <- BIC_min[1,1]-1
  HQ_min[1,1] <- HQ_min[1,1]-1
  colnames(AIC_min) <- colnames(BIC_min) <- colnames(HQ_min) <- c("rank", "lag")

## Best rank:
  best_ranks1 <- sapply(list(AIC_min, BIC_min, HQ_min), rownames)
  best_ranks <- as.numeric(gsub("r=", "", best_ranks1))
  names(best_ranks) <- c("AIC", "BIC", "HQ")

## return result
  res <- list(AICs=AICs, BICs=BICs, HQs=HQs, AIC_min=AIC_min, HQ_min=HQ_min, BIC_min=BIC_min, LLs=LLs, best_ranks=best_ranks)
  res$models_list <-  if(returnModels) models_list else NULL
  class(res ) <- "rank.select"
  return(res)

}




#'Selection of the lag with Information criterion.
#'
#'Selection of the cointegrating rank and the lags with Information criterion
#'(AIC, BIC).
#'
#'This function selects the lag according to AIC, BIC and Hannan-Quinn.
#'
#'@param data multivariate time series.
#'@param lag.max Maximum number of lags to investigate.
#'@param include Type of deterministic regressors to include.
#'@param fitMeasure Whether the AIC/BIC should be based on the full likelihood,
#'or just the SSR. See explanations in \code{\link{logLik.VECM}}.
#'@param sameSample Logical. Whether the data should be shortened so that the
#'AIC/BIC are estimated on the same sample. Default to TRUE.
#'@return An object of class \code{\link{rank.select}}, with \sQuote{print} and
#'\sQuote{summary methods}, containing among other the matrices of AIC/BIC/HQ.
#'@author Matthieu Stigler
#'@seealso \code{\link{rank.select}}, the underlying function, to estimate the
#'rank also.
#'
#'\code{\link[vars]{VARselect}} in package \pkg{vars}, does basically the same.
#'@keywords ts
#'@export
#'@examples
#'
#'
#'data(barry)
#'
#'# 
#'rk_sel <- lags.select(barry)
#'rk_sel
#'summary(rk_sel)
#'
#'
lags.select <- function(data, lag.max=10, include = c( "const", "trend","none", "both"), fitMeasure=c("SSR", "LL"), sameSample=TRUE) {
  rank.select(data=data, lag.max=lag.max, r.max=0, include=include, fitMeasure=fitMeasure, sameSample=sameSample) 
}

#' @rdname rank.select
#' @method print rank.select
#' @S3method print rank.select
print.rank.select <- function(x,...){

  AIC_rank_info <- if(nrow(x$AICs)>1) paste("rank=",x$AIC_min[1,"rank"])
  BIC_rank_info <- if(nrow(x$AICs)>1) paste("rank=",x$BIC_min[1,"rank"])
  HQ_rank_info <- if(nrow(x$AICs)>1) paste("rank=",x$HQ_min[1,"rank"])

  cat("Best AIC:", AIC_rank_info,  " lag=", x$AIC_min[1,"lag"],  "\n")
  cat("Best BIC:", BIC_rank_info, " lag=", x$BIC_min[1,"lag"],  "\n")
  cat("Best HQ :",  HQ_rank_info, " lag=", x$HQ_min[1,"lag"],  "\n")
}

#' @rdname rank.select
#' @method summary rank.select
#' @S3method summary rank.select
summary.rank.select <- function(object,...){

  print.rank.select(object)

  cat("\nBest number of lags:\n")
  AIC_minlag<-apply(object$AICs, 1, which.min)
  BIC_minlag<-apply(object$BICs, 1, which.min)
  HQ_minlag<-apply(object$BICs, 1, which.min)

  mat <- rbind(AIC_minlag, BIC_minlag,HQ_minlag)
  dimnames(mat) <- list(c("AIC", "BIC", "HQ"), paste("r", if(nrow(object$BICs)==1) 0 else 0:(nrow(object$BICs)-1), sep="="))
  print(mat)
}




####################################################
####### TESTS
####################################################
if(FALSE){
library(tsDyn)
library(vars)
#data(Canada)


resu <- rank.select(Canada, sameSample=TRUE)
resu
resu$LLs
resu$AICs
summary(resu)

resu_LL <- rank.select(Canada, sameSample=TRUE, fitMeasure="LL")
resu_LL
summary(resu_LL)
resu_LL$AICs

resvar_SSR <- lags.select(Canada, fitMeasure="SSR")
resvar_SSR
summary(resvar_SSR)
resvar_LL <- lags.select(Canada, fitMeasure="LL")
resvar_LL


resvar2 <- VARselect(Canada)
resvar2$criteria[3,]
resvar_SSR$BICs/(nrow(Canada)-10)
resvar_LL$BICs/(nrow(Canada)-10)

AIC(lineVar(Canada, lag=1), fitMeasure="LL")
AIC(lineVar(Canada, lag=1), fitMeasure="SSR")


v<- lineVar(Canada, lag=1)
# AIC(v, fitMeasure="LL")
AIC(v)





### Get AIC from one model?
ve <- VECM(Canada, lag=2, estim="ML")
var_r0 <- lineVar(Canada, lag=2, I="diff")
var_rk_a <- lineVar(Canada, lag=3, I="level")


logLiks <- Vectorize(logLik.VECM, vectorize.args="r")
AICs <- Vectorize(tsDyn:::AIC.VECM, vectorize.args="r")

logLiks(ve, r=0:3)
AICs(ve, r=0:3, fitMeasure="LL")

### rank 0
logLik(ve, r=0)
logLik(var_r0)
AIC(ve, r=0, fitMeasure="LL")
AIC(var_r0, fitMeasure="LL")

## rank 1:3
logLik(ve, r=1)
logLik(ve, r=2)
logLik(ve, r=3)

## rank 4
logLik(ve, r=4)
logLik(var_rk_a)
AIC(ve, r=4, fitMeasure="LL")
AIC(var_rk_a, fitMeasure="LL")

AIC(ve, r=4, fitMeasure="SSR")
AIC(var_rk_a, fitMeasure="SSR")

(AIC(var_rk_a, fitMeasure="LL")-AIC(ve, r=3, fitMeasure="LL"))/AIC(ve, r=3, fitMeasure="LL")
(AIC(var_rk_a, fitMeasure="SSR")-AIC(ve, r=3, fitMeasure="SSR"))/AIC(ve, r=3, fitMeasure="SSR")

AIC(VECM(Canada, lag=1, r=1, estim="ML"), r=2)
AIC(VECM(Canada, lag=1, r=2, estim="ML"))

AIC(VECM(Canada, lag=1, r=1, estim="ML"), r=2, fitMeasure="LL")
AIC(VECM(Canada, lag=1, r=2, estim="ML"), fitMeasure="LL")

####
ve_r1 <- VECM(Canada, lag=1, r=1, estim="ML")
ve_r2 <- VECM(Canada, lag=1, r=2, estim="ML")
ve_r3 <- VECM(Canada, lag=1, r=3, estim="ML")
var_r4 <- lineVar(Canada, lag=2, I="level")

AIC(ve_r1)-AIC(ve_r2)
AIC(ve_r1, fitMeasure="LL")-AIC(ve_r2, fitMeasure="LL")
AIC(ve_r1, fitMeasure="LL")-AIC(ve_r1,r=2, fitMeasure="LL")

AIC(var_r4)-AIC(ve_r3)
AIC(var_r4, fitMeasure="LL")-AIC(ve_r3, fitMeasure="LL")
AIC(ve_r3,r=4, fitMeasure="LL")-AIC(ve_r3, fitMeasure="LL")

rs <- rank.select(Canada, lag.max=6)$AICs
rs_LL <- rank.select(Canada, lag.max=6, fitMeasure="LL")$AICs
all.equal(rs[1,]-rs[2,], rs_LL[1,]-rs_LL[2,])
all.equal(rs[2,]-rs[3,], rs_LL[2,]-rs_LL[3,])
all.equal(rs[3,]-rs[4,], rs_LL[3,]-rs_LL[4,])
all.equal(rs[4,]-rs[5,], rs_LL[4,]-rs_LL[5,])


rank.select(Canada, lag.max=6, sameSample=FALSE)$AICs
AICs(VECM(Canada, lag=2, estim="ML"), r=0:4)

rank.select(Canada, lag.max=6, sameSample=FALSE, fitMeasure="LL")$AICs
AICs(VECM(Canada, lag=2, estim="ML"), r=0:4, fitMeasure="LL")






tsDyn:::npar.VECM(ve, r=4)
tsDyn:::npar.nlVar(var_rk_a)

-2*logLik(ve, r=3) +2*tsDyn:::npar.VECM(ve, r=3)
-2*logLik(var_rk_a) +2*(tsDyn:::npar.nlVar(var_rk_a)-1)

ve_1 <- VECM(Canada, lag=1, r=1, estim="ML")
round(ve_1$model.specific$lambda,5)



logLik(lineVar(Canada, lag=2, I="diff"))

logLik(VECM(Canada, lag=2, r=1, estim="ML"), r=1)

logLik(VECM(Canada, lag=2, r=1, estim="ML"),r=2)


AIC(VECM(Canada, lag=2, r=1, estim="ML"),r=0, fitMeasure="LL")
AIC(lineVar(Canada, lag=2, I="diff"), fitMeasure="LL")


## remember: ADF(1) = diff(2) !!
deviance(linear(x=sunspots, m=2, type="level"))
deviance(linear(x=sunspots, m=1, type="ADF"))

deviance(lineVar(Canada, lag=2, I="level"))
deviance(lineVar(Canada, lag=1, I="ADF"))



}
