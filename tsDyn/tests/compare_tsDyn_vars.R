
library(tsDyn)
library(vars)

data(Canada)


#########################
##### VECM #####
#########################

### unrestricted cons
vecm_l1_co_tsD <-VECM(Canada, lag=1, include="const", estim="ML")
vecm_l3_co_tsD <-VECM(Canada, lag=3, include="const", estim="ML")

vecm_l1_co_var <- ca.jo(Canada, K=2, ecdet="none", spec="transitory")
vecm_l3_co_var <- ca.jo(Canada, K=4, ecdet="none", spec="transitory")

### restricted cons
vecm_l1_LRco_tsD <-VECM(Canada, lag=1, LRinclude="const", estim="ML")
vecm_l1_LRco_var <- ca.jo(Canada, K=2, ecdet="const", spec="transitory")

vecm_l3_LRco_tsD <-VECM(Canada, lag=3, LRinclude="const", estim="ML")
vecm_l3_LRco_var <- ca.jo(Canada, K=4, ecdet="const", spec="transitory")

### restricted trend
vecm_l1_LRtr_tsD <-VECM(Canada, lag=1, LRinclude="trend", estim="ML")
vecm_l1_LRtr_var <- ca.jo(Canada, K=2, ecdet="trend", spec="transitory")

vecm_l3_LRtr_tsD <-VECM(Canada, lag=3, LRinclude="trend", estim="ML")
vecm_l3_LRtr_var <- ca.jo(Canada, K=4, ecdet="trend", spec="transitory")

all_models <- list(
		    list(vecm_l1_co_var, vecm_l1_co_tsD), 
		    list(vecm_l3_co_var, vecm_l3_co_tsD), 
		    list(vecm_l1_LRco_var, vecm_l1_LRco_tsD),
		    list(vecm_l3_LRco_var, vecm_l3_LRco_tsD),
		    list(vecm_l1_LRtr_var, vecm_l1_LRtr_tsD),
		    list(vecm_l3_LRtr_var, vecm_l3_LRtr_tsD))

deftol <- .Machine$double.eps ^ 0.5
lowtol <- 1.e-07
lowlowtol <- 1.e-06
comp_teststat <- function(x, tol=deftol) all.equal(x[[1]]@teststat, rev(rank.test(x[[2]])$res_df[,"eigen"]), check.attributes=FALSE, tolerance=tol)
comp_betas <- function(x, tol=deftol) all.equal(cajorls(x[[1]])$beta, x[[2]]$model.specific$beta, check.attributes=FALSE, tolerance=tol)
comp_coefs <- function(x, tol=deftol) all.equal(coefficients(cajorls(x[[1]])$rlm), t(coefficients(x[[2]])), check.attributes=FALSE, tolerance=tol)
comp_LL <- function(x) all.equal(as.numeric(logLik(vec2var(x[[1]]))), logLik(x[[2]]), check.attributes=FALSE)
comp_IRF <- function(x) all.equal(irf(vec2var(x[[1]]), boot=FALSE)$irf, irf(x[[2]], boot=FALSE)$irf, check.attributes=FALSE)
comp_IRF_rand <- function(x) all.equal(irf(vec2var(x[[1]]), runs=2, seed=1234)$irf, irf(x[[2]], runs=2, seed=1234)$irf, check.attributes=FALSE)
comp_FEVD <- function(x) all.equal(fevd(vec2var(x[[1]])), fevd(x[[2]]), check.attributes=FALSE)
comp_resid <- function(x, tol=deftol) all.equal(residuals(vec2var(x[[1]])), residuals(x[[2]]), check.attributes=FALSE, tolerance=tol)
comp_fitted <- function(x) all.equal(fitted(vec2var(x[[1]])), fitted(x[[2]], level="original"), check.attributes=FALSE)
comp_predictOld <- function(x) all.equal(predict(vec2var(x[[1]]))$fcst, tsDyn:::predictOld.VECM(x[[2]])$fcst, check.attributes=FALSE)
comp_predict <- function(x) all.equal(sapply(predict(vec2var(x[[1]]), n.ahead=5)$fcst,function(x) x[,"fcst"]), predict(x[[2]]), check.attributes=FALSE)
comp_VECM_coefA <- function(x,...) all.equal(coefA(x[[1]]), coefA(x[[2]]), check.attributes=FALSE,...)
comp_VECM_coefB <- function(x,...) all.equal(coefB(x[[1]]), coefB(x[[2]]), check.attributes=FALSE,...)
comp_VECM_coefPI <- function(x) all.equal(coefPI(x[[1]]), coefPI(x[[2]]), check.attributes=FALSE)
comp_VECM_logLik <- function(x, r=1) all.equal(logLik(x[[1]], r=r), logLik(x[[2]], r=r), check.attributes=FALSE)

### Small function to print nicely output of all.equal, rounding the number:

roundAll.Equal <- function(x, round=8){
  isFALSE <- x!="TRUE"
  xFalse <- x[isFALSE]
  # extract the number (i.e remove all the rest)
  xf<- gsub("(Component ([0-9]+)?([[:punct:]][[:alnum:]]+[[:punct:]])?: )?Mean relative difference: ", 
            "", xFalse)
  xf2<- round(as.numeric(xf),round)
  x[isFALSE] <- paste("Mean relative difference: ", xf2, sep="")
  x
}


### Compare VECM methods:
print(sapply(all_models, comp_teststat ))
print(sapply(all_models, comp_betas, tol=lowtol)) # 2 and 6
roundAll.Equal(sapply(all_models, comp_coefs), round=7) # 5 and 6
print(sapply(all_models, comp_LL)) # 2 and 6
print(sapply(all_models, comp_IRF))
print(sapply(all_models, comp_IRF_rand))
print(sapply(all_models, comp_FEVD))
print(sapply(all_models, comp_resid, tol=lowtol)) # 5 and 6
print(sapply(all_models, comp_fitted)) 
roundAll.Equal(sapply(all_models, comp_predict)) # 5 and 6
lapply(sapply(all_models, comp_predictOld),roundAll.Equal, round=7) # 5 and 6
sapply(all_models, comp_VECM_coefA, tolerance=1e-07)
sapply(all_models, comp_VECM_coefB, tolerance=1e-07)
sapply(all_models, comp_VECM_coefPI)
mapply(function(r) sapply(all_models, comp_VECM_logLik, r=r), r=0:4)

## restricted betas:
all.equal(coefA(vecm_l1_co_tsD),coefA(vecm_l1_co_var), check.attributes=FALSE)
all.equal(coefB(vecm_l1_co_tsD),coefB(vecm_l1_co_var), check.attributes=FALSE)

R <- matrix(c(1,0.1, -0.24, 3.6), ncol=1)

vecm_Rrest_tsD <-VECM(Canada, lag=1, include="const", estim="ML", beta=R)
vecm_Rrest_var <-blrtest(vecm_l1_co_var, H=R, r=1)

all.equal(coefA(vecm_Rrest_tsD),coefA(vecm_Rrest_var), check.attributes=FALSE)
all.equal(coefB(vecm_Rrest_tsD),coefB(vecm_Rrest_var), check.attributes=FALSE)
all.equal(logLik(vecm_Rrest_tsD, r=1),logLik(vecm_Rrest_var, r=1), check.attributes=FALSE)

all.equal(deviance(vecm_Rrest_tsD),sum(deviance(cajorls(vecm_Rrest_var)[[1]])), check.attributes=FALSE)
all.equal(-2*(logLik(vecm_Rrest_tsD)-logLik(vecm_l1_co_tsD)), vecm_Rrest_var@teststat)

#########################
##### VAR #####
#########################

var_l1_co_tsD <-lineVar(Canada, lag=1, include="const")
var_l1_tr_tsD <-lineVar(Canada, lag=1, include="trend")
var_l1_bo_tsD <-lineVar(Canada, lag=1, include="both")
var_l1_no_tsD <-lineVar(Canada, lag=1, include="none")

var_l3_co_tsD <-lineVar(Canada, lag=3, include="const")
var_l3_tr_tsD <-lineVar(Canada, lag=3, include="trend")
var_l3_bo_tsD <-lineVar(Canada, lag=3, include="both")
var_l3_no_tsD <-lineVar(Canada, lag=3, include="none")

var_l1_co_var <- VAR(Canada, p=1, type="const")
var_l1_tr_var <- VAR(Canada, p=1, type="trend")
var_l1_bo_var <- VAR(Canada, p=1, type="both")
var_l1_no_var <- VAR(Canada, p=1, type="none")

var_l3_co_var <- VAR(Canada, p=3, type="const")
var_l3_tr_var <- VAR(Canada, p=3, type="trend")
var_l3_bo_var <- VAR(Canada, p=3, type="both")
var_l3_no_var <- VAR(Canada, p=3, type="none")



all_var_models <- list(
		    list(var_l1_co_tsD, var_l1_co_var), 
		    list(var_l1_tr_tsD, var_l1_tr_var), 
		    list(var_l1_bo_tsD, var_l1_bo_var), 
		    list(var_l1_no_tsD, var_l1_no_var), 
		    list(var_l3_co_tsD, var_l3_co_var), 
		    list(var_l3_tr_tsD, var_l3_tr_var), 
		    list(var_l3_bo_tsD, var_l3_bo_var), 
		    list(var_l3_no_tsD, var_l3_no_var))

all_var_models_noNoBo <- all_var_models <- list(
		    list(var_l1_co_tsD, var_l1_co_var), 
		    list(var_l1_tr_tsD, var_l1_tr_var), 
		    list(var_l3_co_tsD, var_l3_co_var), 
		    list(var_l3_tr_tsD, var_l3_tr_var))

## test functions
coef_to_vars <- function(x){
  c <- coef(x)
  if(any(grepl("Intercept|Trend", colnames(c)))){
    wh.deter <- grep("Intercept|Trend", colnames(c))
    res <- cbind(c[,-wh.deter], c[,wh.deter,drop=FALSE])
  } else {
    res <- c
  }
  colnames(res) <-gsub(" -", "\\.l", colnames(res))
  rownames(res) <-gsub("Equation ", "", rownames(res))
  return(res)
}



comp_var_coefs <- function(x) all.equal(coef_to_vars (x[[1]]), t(sapply(coef(x[[2]]), function(x) x[,"Estimate"])), check.attributes=FALSE)
comp_var_logLik <- function(x, tol=deftol) all.equal(logLik(x[[1]]), as.numeric(logLik(x[[2]])), check.attributes=FALSE, tolerance=tol)

comp_var_pred <- function(x) all.equal(predict(x[[1]]), sapply(predict(x[[2]], n.ahead=5)$fcst, function(x) x[,"fcst"]),check.attributes=FALSE)
comp_var_predOld <- function(x) all.equal(sapply(tsDyn:::predictOld.VAR(x[[1]], n.ahead=5)$fcst, function(x) x[,"fcst"]), sapply(predict(x[[2]], n.ahead=5)$fcst, function(x) x[,"fcst"]),check.attributes=FALSE)
comp_var_fevd <- function(x) all.equal(sapply(fevd(x[[1]]), head,2), sapply(fevd(x[[2]]), head,2), check.attributes=FALSE)
comp_var_IRF <- function(x) isTRUE(all.equal(irf(x[[1]], boot=FALSE)$irf, irf(x[[2]], boot=FALSE)$irf, check.attributes=FALSE))

### Compare VECM methods:
sapply(all_var_models, comp_var_coefs)
sapply(all_var_models, comp_var_logLik, tol=lowlowtol)
roundAll.Equal(sapply(all_var_models_noNoBo, comp_var_pred),7)
roundAll.Equal(sapply(all_var_models_noNoBo, comp_var_predOld),7)
sapply(all_var_models_noNoBo, comp_var_IRF)

