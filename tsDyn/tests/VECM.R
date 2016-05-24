library(tsDyn)


data(zeroyld)
data(barry)




## Test a few VECM models
vecm_OLS_l1_co <-VECM(barry, lag=1)
vecm_OLS_l3_co <-VECM(barry, lag=3, include="const")
vecm_OLS_l3_co_betaGiven<-VECM(barry, lag=3, include="const", beta=c(0.1, -0.05))
vecm_OLS_l1_tr <-VECM(barry, lag=1, include="trend")
vecm_OLS_l1_bo <-VECM(barry, lag=1, include="both")
vecm_OLS_l1_no <-VECM(barry, lag=1, include="none")
vecm_OLS_l1_coAsExo <-VECM(barry, lag=1, include="none", exogen=rep(1, nrow(barry)))
vecm_OLS_l3_coAsExo <-VECM(barry, lag=3, include="none", exogen=rep(1, nrow(barry)))

vecm_OLS_l0_co <-VECM(barry, lag=0)
vecm_OLS_l0_tr <-VECM(barry, lag=0, include="trend")
vecm_ML_l0_co <-VECM(barry, lag=0, include="const", estim="ML")
vecm_ML_l0_tr <-VECM(barry, lag=0, include="trend", estim="ML")
vecm_ML_l0_bo <-VECM(barry, lag=0, include="both", estim="ML")
vecm_ML_l0_no <-VECM(barry, lag=0, include="none", estim="ML")
vecm_ML_l0_LRco <-VECM(barry, lag=0, LRinclude="const", estim="ML")
vecm_ML_l0_LRtr <-VECM(barry, lag=0, LRinclude="trend", estim="ML")

vecm_OLS_l1_LRco <-VECM(barry, lag=1, LRinclude="const")
vecm_OLS_l1_LRtr <-VECM(barry, lag=1, LRinclude="trend")
vecm_OLS_l1_LRtr_noCo <-VECM(barry, lag=1, LRinclude="trend", include="none")
vecm_OLS_l1_LRbo <-VECM(barry, lag=1, LRinclude="both")

vecm_ML_l1_co <-VECM(barry, lag=1, estim="ML")
vecm_ML_l3_co <-VECM(barry, lag=3, include="const", estim="ML")
vecm_ML_l3_co_betaGiven<-VECM(barry, lag=3, include="const", beta=c(-0.035,0.04), estim="ML")
vecm_ML_l1_tr <-VECM(barry, lag=1, include="trend", estim="ML")
vecm_ML_l1_bo <-VECM(barry, lag=1, include="both", estim="ML")
vecm_ML_l1_no <-VECM(barry, lag=1, include="none", estim="ML")
vecm_ML_l1_coAsExo <-VECM(barry, lag=1, include="none", exogen=rep(1, nrow(barry)), estim="ML")
vecm_ML_l3_coAsExo <-VECM(barry, lag=3, include="none", exogen=rep(1, nrow(barry)), estim="ML")

set.seed(1234)
exoVar <- rnorm(n=nrow(barry))
vecm_ML_l1_LRco <-VECM(barry, lag=1, LRinclude="const", estim="ML")
vecm_ML_l1_LRc_exo <-VECM(barry, lag=1, LRinclude="const", estim="ML", exogen=exoVar)
vecm_ML_l1_LRtr <-VECM(barry, lag=1, LRinclude="trend", estim="ML")
vecm_ML_l1_LRtr_exo <-VECM(barry, lag=1, LRinclude="trend", estim="ML", exogen=exoVar)
vecm_ML_l1_LRtr_noCo <-VECM(barry, lag=1, LRinclude="trend", include="none", estim="ML")
vecm_ML_l1_LRbo <-VECM(barry, lag=1, LRinclude="both", estim="ML")
vecm_ML_l1_LRbo_exo <-VECM(barry, lag=1, LRinclude="both", estim="ML", exogen=exoVar)

vecm_all <- list(
		vecm_OLS_l1_co=vecm_OLS_l1_co,
		vecm_OLS_l3_co=vecm_OLS_l3_co,
		vecm_OLS_l3_co_betaGiven=vecm_OLS_l3_co_betaGiven,
		vecm_OLS_l1_tr=vecm_OLS_l1_tr, 
		vecm_OLS_l1_bo=vecm_OLS_l1_bo,
		vecm_OLS_l1_no=vecm_OLS_l1_no,
		vecm_OLS_l1_coAsExo=vecm_OLS_l1_coAsExo,
		vecm_OLS_l3_coAsExo=vecm_OLS_l3_coAsExo, 
		vecm_OLS_l1_LRco=vecm_OLS_l1_LRco,
		vecm_OLS_l1_LRtr=vecm_OLS_l1_LRtr, 
		vecm_OLS_l1_LRtr_noCo=vecm_OLS_l1_LRtr_noCo, 
		vecm_OLS_l1_LRbo=vecm_OLS_l1_LRbo, 

		vecm_OLS_l0_co=vecm_OLS_l0_co,
		vecm_OLS_l0_tr=vecm_OLS_l0_tr,
		vecm_ML_l0_co=vecm_ML_l0_co,
		vecm_ML_l0_tr=vecm_ML_l0_tr,
		vecm_ML_l0_bo=vecm_ML_l0_bo,
		vecm_ML_l0_no=vecm_ML_l0_no,
		vecm_ML_l0_LRco=vecm_ML_l0_LRco,
		vecm_ML_l0_LRtr=vecm_ML_l0_LRtr,
		
    
		vecm_ML_l1_co=vecm_ML_l1_co,
		vecm_ML_l3_co=vecm_ML_l3_co,
		vecm_ML_l1_tr=vecm_ML_l1_tr, 
		vecm_ML_l1_bo=vecm_ML_l1_bo,
		vecm_ML_l1_no=vecm_ML_l1_no,
		vecm_ML_l1_coAsExo=vecm_ML_l1_coAsExo,
		vecm_ML_l3_coAsExo=vecm_ML_l3_coAsExo, 
		vecm_ML_l1_LRco=vecm_ML_l1_LRco,
		vecm_ML_l1_LRc_exo=vecm_ML_l1_LRc_exo,
		vecm_ML_l1_LRtr=vecm_ML_l1_LRtr,
		vecm_ML_l1_LRtr_exo=vecm_ML_l1_LRtr_exo,
		vecm_ML_l1_LRtr_noCo=vecm_ML_l1_LRtr_noCo,
		vecm_ML_l1_LRbo=vecm_ML_l1_LRbo,
		vecm_ML_l1_LRbo_exo=vecm_ML_l1_LRbo_exo)


vecm_ML <- vecm_all[grep("ML", names(vecm_all))]
vecm_no_l0 <- vecm_all[!names(vecm_all)%in%grep("l0", names(vecm_all), value=TRUE)]

lapply(vecm_all, print)
lapply(vecm_all, summary)

lapply(vecm_all, function(x) head(residuals(x), 3))
lapply(vecm_all, function(x) head(fitted(x), 3))
lapply(vecm_no_l0, function(x) head(fitted(x, level="original"), 3))
sapply(vecm_all, deviance)



## logLik
sapply(vecm_all, logLik)
sapply(vecm_ML, logLik, r=0)
sapply(vecm_ML, logLik, r=1)
sapply(vecm_ML, logLik, r=2)

## AIC/BIC
sapply(vecm_all, AIC)
sapply(vecm_ML, AIC, r=0)
sapply(vecm_ML, AIC, r=1)
sapply(vecm_ML, AIC, r=2)
sapply(vecm_ML, AIC, r=0, fitMeasure="LL")
sapply(vecm_ML, AIC, r=1, fitMeasure="LL")
sapply(vecm_ML, AIC, r=2, fitMeasure="LL")

sapply(vecm_all, BIC)
sapply(vecm_ML, BIC, r=0)
sapply(vecm_ML, BIC, r=0, fitMeasure="LL")


## coint
sapply(vecm_all, function(x) x$model.specific$beta)

### VARrep
lapply(vecm_all, function(x) round(VARrep(x),9))

### fevd
lapply(vecm_all, function(x) sapply(fevd(x, n.ahead=2), head))

### irf
vecm_irf <- vecm_all[-grep("l1_no|bo|exo|Exo|l0", names(vecm_all))] ## does not work for these models
lapply(vecm_irf, function(x) sapply(irf(x, runs=1)$irf,head,2))

## predict
vecm_all_pred <- vecm_all[-grep("_bo|_no|_noCo|LRbo|coAsExo|exo", names(vecm_all))]
lapply(vecm_all_pred, predict,  n.ahead=2)
lapply(vecm_all_pred, function(x) sapply(tsDyn:::predictOld.VECM(x, n.ahead=2)$fcst, function(y) y[,"fcst"]))

lapply(vecm_all, function(x) predict_rolling(x,nroll=2)$pred)
lapply(vecm_all, function(x) predict_rolling(x,nroll=2, refit.every=1)$pred)

## VECM boot
vecm_all_noExo_noLRinc <- vecm_all[-grep("LR|exo|Exo|l0", names(vecm_all))]
options(warn=-1)
sapply(vecm_all_noExo_noLRinc, function(x) try(tsDyn:::check.VECM.boot(x), silent=TRUE))
options(warn=0)

### CoefA, coefB, coefPI
lapply(vecm_all, coefB)
lapply(vecm_all, coefA)
options(digits=6)
lapply(vecm_all, coefPI)
options(digits=7)

### rank test
vecm_ML_rtest <- vecm_ML[-grep("vecm_ML_l1_LRtr_noCo|vecm_ML_l1_LRbo", names(vecm_ML))] ## does not work for these models

rank.tests <- lapply(vecm_ML_rtest , rank.test)
rank.tests_rnull1 <- lapply(vecm_ML_rtest , rank.test, r_null=1)
rank.tests_tr <- lapply(vecm_ML_rtest , rank.test, type="trace")
rank.tests_tr_rnull1 <- lapply(vecm_ML_rtest , rank.test, r_null=1, type="trace")

rank.tests.all <- c(rank.tests , rank.tests_rnull1, rank.tests_tr,rank.tests_tr_rnull1 )

lapply(rank.tests.all, print)
lapply(rank.tests.all, summary)


### rank select
data(barry)
r_sel <- rank.select(barry)
r_sel_tre <- rank.select(barry, include="trend")
r_sel_none <- rank.select(barry, include="none")
r_sel_both <- rank.select(barry, include="both")

r_sel$LLs
r_sel$AICs

r_sel_tre$LLs
r_sel_tre$AICs

r_sel_none$LLs
r_sel_none$AICs

r_sel_both$LLs
r_sel_both$AICs

#### exogen: check equalities
check.same <- function(x1, x2) {
  co_x2 <- coef(x2)
  t1 <- isTRUE(all.equal(coef(x1), co_x2[,c(1:x2$model.specific$r, ncol(co_x2),(x2$model.specific$r+1):(ncol(co_x2)-1))], check.attributes=FALSE))
  t2 <- isTRUE(all.equal(AIC(x1), AIC(x2), check.attributes=FALSE))
  t3 <- isTRUE(all.equal(BIC(x1), BIC(x2), check.attributes=FALSE))
  if(x1$model.specific$estim=="ML"){
    t5 <- isTRUE(all.equal(BIC(x1,fitMeasure="LL"), BIC(x2,fitMeasure="LL"), check.attributes=FALSE))
    t6 <- isTRUE(all.equal(BIC(x1,fitMeasure="LL", r=2), BIC(x2,fitMeasure="LL",r=2), check.attributes=FALSE))
    t7 <- isTRUE(all.equal(logLik(x1,fitMeasure="LL", r=2), logLik(x2,fitMeasure="LL",r=2), check.attributes=FALSE))
    t8 <- isTRUE(all.equal(rank.test(x1)$res_df[,c("trace", "eigen")], rank.test(x2)$res_df[,c("trace", "eigen")], check.attributes=FALSE))
  } else {
    t5 <- t6 <- t7 <- t8 <-  NULL
  }
  t4 <- isTRUE(all.equal(residuals(x1), residuals(x2), check.attributes=FALSE))
  c(t1, t1, t3,t4, t5, t6, t7, t8)
}

check.same(x1=vecm_OLS_l1_co, x2=vecm_OLS_l1_coAsExo)
check.same(x1=vecm_OLS_l3_co, x2=vecm_OLS_l3_coAsExo)
check.same(x1=vecm_ML_l1_co, x2=vecm_ML_l1_coAsExo)
check.same(x1=vecm_ML_l3_co, x2=vecm_ML_l3_coAsExo)



###################################
####### predict_rolling check
###################################
n_ca<- nrow(barry)
as.M <- function(x) as.matrix(x)

#### VECM No refit lag=1
mod_vecm_l1_full <- VECM(barry, lag=1)
mod_vecm_l1_sub <- VECM(tsDyn:::myHead(barry,n_ca-10), lag=1)
mod_vecm_l1_sub_ref5 <- VECM(tsDyn:::myHead(barry,n_ca-5), lag=1)


pred_l3_vecm_roll_12<-predict_rolling(object=mod_vecm_l1_full, nroll=10, n.ahead=1:2)
pred_l3_vecm_0_12_nd <- predict(mod_vecm_l1_sub, n.ahead=2, newdata=barry[n_ca-c(12:11),,drop=FALSE])
pred_l3_vecm_1_12_nd_ref5 <- predict(mod_vecm_l1_sub_ref5, n.ahead=2, newdata=barry[n_ca-c(6:5),,drop=FALSE])
pred_l3_vecm_0_12_nd_ref5 <- predict(mod_vecm_l1_sub_ref5, n.ahead=2, newdata=barry[n_ca-c(7:6),,drop=FALSE])
pred_l3_vecm_1_12_nd <- predict(mod_vecm_l1_sub, n.ahead=2, newdata=barry[n_ca-c(11:10),,drop=FALSE])
pred_l3_vecm_2_12_nd <- predict(mod_vecm_l1_sub, n.ahead=2, newdata=barry[n_ca-c(10:9),,drop=FALSE])
pred_l3_vecm_1_12 <- predict(mod_vecm_l1_sub, n.ahead=2)
all.equal(pred_l3_vecm_1_12_nd, pred_l3_vecm_1_12) ## minor: consistency in predict with/withotut newdata=dataset


pred_l3_vecm_nd <- rbind(pred_l3_vecm_0_12_nd, pred_l3_vecm_1_12_nd, pred_l3_vecm_2_12_nd)
all.equal(pred_l3_vecm_nd[c(3,5),], as.M(pred_l3_vecm_roll_12$pred[1:2,1:3]), check.attributes=FALSE) ## check 1-ahead
all.equal(pred_l3_vecm_nd[c(2,4),], as.M(pred_l3_vecm_roll_12$pred[11:12,1:3]), check.attributes=FALSE) ## check 2 ahead



#### VECM No refit lag=3
mod_var_l3_vecm_full <- VECM(barry, lag=3)
mod_var_l3_vecm_sub <- VECM(tsDyn:::myHead(barry,n_ca-10), lag=3)
mod_var_l3_vecm_sub_ref5 <- VECM(tsDyn:::myHead(barry,n_ca-6), lag=3)

pred_l3_vecm_roll_12<-predict_rolling(object=mod_var_l3_vecm_full, nroll=10, n.ahead=1:2)
pred_l3_vecm_0_12_nd <- predict(mod_var_l3_vecm_sub, n.ahead=2, newdata=barry[n_ca-c(14:11),,drop=FALSE])
pred_l3_vecm_1_12_nd_ref5 <- predict(mod_var_l3_vecm_sub_ref5, n.ahead=2, newdata=barry[n_ca-c(6:9),,drop=FALSE])
pred_l3_vecm_0_12_nd_ref5 <- predict(mod_var_l3_vecm_sub_ref5, n.ahead=2, newdata=barry[n_ca-c(7:10),,drop=FALSE])
pred_l3_vecm_1_12_nd <- predict(mod_var_l3_vecm_sub, n.ahead=2, newdata=barry[n_ca-c(13:10),,drop=FALSE])
pred_l3_vecm_2_12_nd <- predict(mod_var_l3_vecm_sub, n.ahead=2, newdata=barry[n_ca-c(12:9),,drop=FALSE])
pred_l3_vecm_1_12 <- predict(mod_var_l3_vecm_sub, n.ahead=2)
all.equal(pred_l3_vecm_1_12_nd, pred_l3_vecm_1_12) ## minor: consistency in predict with/withotut newdata=dataset


pred_l3_vecm_nd <- rbind(pred_l3_vecm_0_12_nd, pred_l3_vecm_1_12_nd, pred_l3_vecm_2_12_nd)
all.equal(pred_l3_vecm_nd[c(3,5),], as.M(pred_l3_vecm_roll_12$pred[1:2,1:3]), check.attributes=FALSE) ## check 1-ahead
all.equal(pred_l3_vecm_nd[c(2,4),], as.M(pred_l3_vecm_roll_12$pred[11:12,1:3]), check.attributes=FALSE) ## check 2 ahead


### VECM Refit lag=3
pred_l3_vecm_roll_12_ref<-predict_rolling(object=mod_var_l3_vecm_full, nroll=10, n.ahead=1:2, refit=5)
pred_l3_vecm_roll_12_ref_b<-predict_rolling(object=mod_var_l3_vecm_full, nroll=10, n.ahead=2, refit=5)


pred_l3_vecm_0_12_nd_ref5 <- predict(mod_var_l3_vecm_sub_ref5, n.ahead=2, newdata=barry[n_ca-c(9:6),,drop=FALSE])
pred_l3_vecm_1_12_nd_ref5 <- predict(mod_var_l3_vecm_sub_ref5, n.ahead=2, newdata=barry[n_ca-c(8:5),,drop=FALSE])

pred_l3_vecm_nd_ref5 <- rbind(pred_l3_vecm_0_12_nd_ref5, pred_l3_vecm_1_12_nd_ref5)
all.equal(pred_l3_vecm_nd_ref5[c(2,4),], as.M(pred_l3_vecm_roll_12_ref$pred[16:17,1:3]), check.attributes=FALSE) ## check refited2 ahead


