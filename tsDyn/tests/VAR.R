library(tsDyn)


data(zeroyld)
data(barry)




## Test a few VAR models
var_l1_co <-lineVar(barry, lag=1, include="const")
var_l1_tr <-lineVar(barry, lag=1, include="trend")
var_l1_bo <-lineVar(barry, lag=1, include="both")
var_l1_no <-lineVar(barry, lag=1, include="none")
var_l1_coAsExo <-lineVar(barry, lag=1, include="none", exogen=rep(1, nrow(barry)))

var_l0_co <-lineVar(barry, lag=0, include="const")
var_l0_tr <-lineVar(barry, lag=0, include="trend")
var_l0_bo <-lineVar(barry, lag=0, include="both")
var_l0_coAsExo <-lineVar(barry, lag=0, include="trend", exogen=rep(1, nrow(barry)))

var_l3_co <-lineVar(barry, lag=3, include="const")
var_l3_tr <-lineVar(barry, lag=3, include="trend")
var_l3_bo <-lineVar(barry, lag=3, include="both")
var_l3_no <-lineVar(barry, lag=3, include="none")
var_l3_coAsExo <-lineVar(barry, lag=3, include="none", exogen=rep(1, nrow(barry)))

var_l2_diff_co <-lineVar(barry, lag=2, include="const", I="diff")
var_l2_diff_tr <-lineVar(barry, lag=2, include="trend", I="diff")
var_l2_diff_bo <-lineVar(barry, lag=2, include="both", I="diff")
var_l2_diff_no <-lineVar(barry, lag=2, include="none", I="diff")
var_l2_diff_coAsExo <-lineVar(barry, lag=2, include="none", I="diff", exogen=rep(1, nrow(barry)))

var_l2_adf_co <-lineVar(barry, lag=2, include="const", I="ADF")
var_l2_adf_tr <-lineVar(barry, lag=2, include="trend", I="ADF")
var_l2_adf_bo <-lineVar(barry, lag=2, include="both", I="ADF")
var_l2_adf_no <-lineVar(barry, lag=2, include="none", I="ADF")
var_l2_adf_coAsExo <-lineVar(barry, lag=2, include="none", I="ADF", exogen=rep(1, nrow(barry)))


var_all <- list(
		var_l1_co, var_l1_tr, var_l1_bo, var_l1_no, var_l1_coAsExo,
		var_l0_co, var_l0_tr, var_l0_bo, var_l0_coAsExo,
		var_l3_co, var_l3_tr, var_l3_bo, var_l3_no, var_l3_coAsExo,
		var_l2_diff_co, var_l2_diff_tr, var_l2_diff_bo, var_l2_diff_no, var_l2_diff_coAsExo,
		var_l2_adf_co, var_l2_adf_tr, var_l2_adf_bo, var_l2_adf_no, var_l2_adf_coAsExo)


names(var_all) <-c(
		"var_l1_co", "var_l1_tr", "var_l1_bo", "var_l1_no", "var_l1_coAsExo",
		"var_l0_co", "var_l0_tr", "var_l0_bo",  "var_l0_coAsExo",
		"var_l3_co", "var_l3_tr", "var_l3_bo", "var_l3_no", "var_l3_coAsExo",
		"var_l2_diff_co", "var_l2_diff_tr", "var_l2_diff_bo", "var_l2_diff_no", "var_l2_diff_coAsExo",
		"var_l2_adf_co", "var_l2_adf_tr", "var_l2_adf_bo", "var_l2_adf_no","var_l2_adf_coAsExo")


## Check methods:
lapply(var_all, print)
lapply(var_all, summary)

lapply(var_all, function(x) head(residuals(x), 3))
lapply(var_all, function(x) head(fitted(x), 3))
sapply(var_all, deviance)


## logLik/AIC/BIC
sapply(var_all, logLik)
sapply(var_all, AIC)
sapply(var_all, AIC, fitMeasure="LL")
sapply(var_all, BIC)
sapply(var_all, BIC, fitMeasure="LL")

### VARrep
var_all_noADF <- var_all[-grep("adf", names(var_all))]
lapply(var_all_noADF, VARrep)


### fevd
var_all_level <- var_all[-grep("diff|adf|Exo|l0", names(var_all))]
lapply(var_all_level , function(x) sapply(fevd(x, n.ahead=2), head))


## predict
var_all_pred <- var_all[-grep("bo|no|adf|diff|Exo|l0", names(var_all))]
var_all_pred2 <- var_all[-grep("adf|diff|Exo", names(var_all))]
lapply(var_all_pred2, predict, n.ahead=2)
lapply(var_all, function(x) try(predict(x, n.ahead=2), silent=TRUE))
lapply(var_all_pred, function(x) sapply(tsDyn:::predictOld.VAR(x, n.ahead=2)$fcst, function(y) y[,"fcst"]))
lapply(var_all, function(x) try(sapply(tsDyn:::predictOld.VAR(x, n.ahead=2)$fcst, function(y) y[,"fcst"]), silent=TRUE))

all.equal(lapply(var_all_pred, predict, n.ahead=2), lapply(var_all_pred, function(x) sapply(tsDyn:::predictOld.VAR(x, n.ahead=2)$fcst, function(y) y[,"fcst"])), check.attributes=FALSE)

lapply(var_all_level , function(x) predict_rolling(x,nroll=2)$pred)
lapply(var_all_level , function(x) predict_rolling(x,nroll=2, refit.every=1)$pred)

## check "retro" predictions against fitted
check.pred <- function(x){
  true <- tail(fitted(x),1)
  if(x$lag>0){
    newD <- barry[nrow(barry)-(x$lag:1),,drop=FALSE] 
    check <- predict(x, newdata=newD, newdataTrendStart=x$t, n.ahead=1)
  } else {
    check <- predict(x, n.ahead=1, newdataTrendStart=x$t)
  }
  
  isTRUE(all.equal(true, check, check.attributes=FALSE))
}
sapply(var_all_pred2, check.pred)


## boot
var_all_boot <- var_all[-grep("adf|diff|Exo|l0", names(var_all))]
lapply(var_all_boot, function(x) tail(VAR.boot(x, seed=1234),2))
checkBoot <- function(x){
  check <- VAR.boot(x, boot.scheme="check")
  all.equal(check, as.matrix(as.data.frame(barry)), check.attributes = FALSE)
}
sapply(var_all_boot, checkBoot)

## sim 
comp_tvar_sim <- function(mod, serie){
  ns <- nrow(serie)
  sim_mod <- TVAR.sim(B=coef(mod), lag=mod$lag, include=mod$include,nthresh=0, n=ns-mod$lag, innov=residuals(mod), starting=serie[1:mod$lag,,drop=FALSE])
  all.equal(sim_mod, as.matrix(serie)[-c(1:mod$lag),], check.attributes=FALSE)
}

lapply(var_all_level, comp_tvar_sim, serie=barry)



#### exogen: check equalities
check.same <- function(x1, x2) {
  co_x2 <- coef(x2)
  t1 <- isTRUE(all.equal(coef(x1), co_x2[,c(ncol(co_x2),1:(ncol(co_x2)-1))], check.attributes=FALSE))
  t2 <- isTRUE(all.equal(AIC(x1), AIC(x2), check.attributes=FALSE))
  t3 <- isTRUE(all.equal(BIC(x1), BIC(x2), check.attributes=FALSE))
  t4 <- isTRUE(all.equal(BIC(x1,fitMeasure="LL"), BIC(x2,fitMeasure="LL"), check.attributes=FALSE))
  t5 <- isTRUE(all.equal(residuals(x1), residuals(x2), check.attributes=FALSE))
  if(attr(x1, "varsLevel")!="ADF"){
    va_x2 <- VARrep(x2)
    t6 <- isTRUE(all.equal(VARrep(x1), va_x2[,c(ncol(va_x2),1:(ncol(va_x2)-1))], check.attributes=FALSE))
  } else {
    t6 <- NULL
  }
  c(t1, t1, t3,t4, t5, t6)
}

check.same(x1=var_l1_co, x2=var_l1_coAsExo)
check.same(x1=var_l3_co, x2=var_l3_coAsExo)
check.same(x1=var_l2_diff_co, x2=var_l2_diff_coAsExo)
check.same(x1=var_l2_adf_co, x2=var_l2_adf_coAsExo)


###################################
####### predict_rolling check
###################################
n_ca<- nrow(barry)
as.M <- function(x) as.matrix(x)

mod_var_l1_full <- lineVar(barry, lag=1)
mod_var_l1_sub <- lineVar(tsDyn:::myHead(barry,n_ca-10), lag=1)
mod_var_l1_sub_ref5 <- lineVar(tsDyn:::myHead(barry,n_ca-5), lag=1)

pred_l1_roll_12<-predict_rolling(object=mod_var_l1_full, nroll=10, n.ahead=1:2)
pred_l1_0_12_nd <- predict(mod_var_l1_sub, n.ahead=2, newdata=barry[n_ca-11,,drop=FALSE])
pred_l1_1_12_nd_ref5 <- predict(mod_var_l1_sub_ref5, n.ahead=2, newdata=barry[n_ca-5,,drop=FALSE])
pred_l1_0_12_nd_ref5 <- predict(mod_var_l1_sub_ref5, n.ahead=2, newdata=barry[n_ca-6,,drop=FALSE])
pred_l1_1_12_nd <- predict(mod_var_l1_sub, n.ahead=2, newdata=barry[n_ca-10,,drop=FALSE])
pred_l1_2_12_nd <- predict(mod_var_l1_sub, n.ahead=2, newdata=barry[n_ca-9,,drop=FALSE])
pred_l1_1_12 <- predict(mod_var_l1_sub, n.ahead=2)
all.equal(pred_l1_1_12_nd, pred_l1_1_12) ## minor: consistency in predict with/withotut newdata=dataset


pred_l1_nd <- rbind(pred_l1_0_12_nd, pred_l1_1_12_nd, pred_l1_2_12_nd)
all.equal(pred_l1_nd[c(3,5),], as.M(pred_l1_roll_12$pred[1:2,1:3]), check.attributes=FALSE) ## check 1-ahead
all.equal(pred_l1_nd[c(2,4),], as.M(pred_l1_roll_12$pred[11:12,1:3]), check.attributes=FALSE) ## check 2 ahead


#### No refit lag=3
mod_var_l3_full <- lineVar(barry, lag=3)
mod_var_l3_sub <- lineVar(tsDyn:::myHead(barry,n_ca-10), lag=3)

pred_l3_0_12_nd <- predict(mod_var_l3_sub, n.ahead=2, newdata=barry[n_ca-c(13:11),,drop=FALSE])
pred_l3_1_12_nd <- predict(mod_var_l3_sub, n.ahead=2, newdata=barry[n_ca-(12:10),,drop=FALSE])
pred_l3_2_12_nd <- predict(mod_var_l3_sub, n.ahead=2, newdata=barry[n_ca-c(11:9),,drop=FALSE])
pred_l3_1_12 <- predict(mod_var_l3_sub, n.ahead=2)
all.equal(pred_l3_1_12_nd, pred_l3_1_12) ## minor: consistency in predict with/withotut newdata=dataset

  
pred_l3_roll_12<-predict_rolling(object=mod_var_l3_full, nroll=10, n.ahead=1:2)
pred_l3_nd <- rbind(pred_l3_0_12_nd, pred_l3_1_12_nd, pred_l3_2_12_nd)
all.equal(pred_l3_nd[c(3,5),], as.M(pred_l3_roll_12$pred[1:2,1:3]), check.attributes=FALSE) ## check 1-ahead
all.equal(pred_l3_nd[c(2,4),], as.M(pred_l3_roll_12$pred[11:12,1:3]), check.attributes=FALSE) ## check 2 ahead


### Refit lag=1
pred_l1_ref_roll_12<-predict_rolling(object=mod_var_l1_full, nroll=10, n.ahead=1:2, refit=5)

mod_var_l1_sub_ref5 <- lineVar(tsDyn:::myHead(barry,n_ca-6), lag=1)
pred_l1_1_12_nd_ref5 <- predict(mod_var_l1_sub_ref5, n.ahead=2, newdata=barry[n_ca-5,,drop=FALSE])
pred_l1_0_12_nd_ref5 <- predict(mod_var_l1_sub_ref5, n.ahead=2, newdata=barry[n_ca-6,,drop=FALSE])

pred_l1_nd_ref5 <- rbind(pred_l1_0_12_nd_ref5, pred_l1_1_12_nd_ref5)
all.equal(pred_l1_nd[c(3,5),], as.M(pred_l1_ref_roll_12$pred[1:2,1:3]), check.attributes=FALSE) ## check 1-ahead
all.equal(pred_l1_nd[c(2,4),], as.M(pred_l1_ref_roll_12$pred[11:12,1:3]), check.attributes=FALSE) ## check 2 ahead
all.equal(pred_l1_nd_ref5[c(2,4),], as.M(pred_l1_ref_roll_12$pred[16:17,1:3]), check.attributes=FALSE) ## check refited2 ahead


### Refit lag=3
mod_var_l3_sub_ref5 <- lineVar(tsDyn:::myHead(barry,n_ca-6), lag=3)

pred_l3_ref_roll_12<-predict_rolling(object=mod_var_l3_full, nroll=10, n.ahead=1:2, refit=5)

all.equal(pred_l3_nd[c(3,5),], as.M(pred_l3_ref_roll_12$pred[1:2,1:3]), check.attributes=FALSE) ## check 1-ahead
all.equal(pred_l3_nd[c(2,4),], as.M(pred_l3_ref_roll_12$pred[11:12,1:3]), check.attributes=FALSE) ## check 2 ahead


pred_l3_0_12_nd_ref5 <- predict(mod_var_l3_sub_ref5, n.ahead=2, newdata=barry[n_ca-c(8:6),,drop=FALSE])
pred_l3_1_12_nd_ref5 <- predict(mod_var_l3_sub_ref5, n.ahead=2, newdata=barry[n_ca-c(7:5),,drop=FALSE])
pred_l3_nd_ref5 <- rbind(pred_l3_0_12_nd_ref5, pred_l3_1_12_nd_ref5)
all.equal(pred_l3_nd_ref5[c(2,4),], as.M(pred_l3_ref_roll_12$pred[16:17,1:3]), check.attributes=FALSE) ## check refited2 ahead


#### No refit: VAR diff,  lag=1
mod_var_l1_diff_full <- lineVar(barry, lag=1, I="diff")
mod_var_l1_diff_sub <- lineVar(tsDyn:::myHead(barry,n_ca-10), lag=1, I="diff")
mod_var_l1_diff_sub_ref5 <- lineVar(tsDyn:::myHead(barry,n_ca-5), lag=1, I="diff")

pred_l1_diff_roll_12<-predict_rolling(object=mod_var_l1_diff_full, nroll=10, n.ahead=1:2)
pred_l1_diff_0_12_nd <- predict(mod_var_l1_diff_sub, n.ahead=2, newdata=barry[n_ca-c(12:11),,drop=FALSE])
pred_l1_diff_1_12_nd_ref5 <- predict(mod_var_l1_diff_sub_ref5, n.ahead=2, newdata=barry[n_ca-c(6:5),,drop=FALSE])
pred_l1_diff_0_12_nd_ref5 <- predict(mod_var_l1_diff_sub_ref5, n.ahead=2, newdata=barry[n_ca-c(7:6),,drop=FALSE])
pred_l1_diff_1_12_nd <- predict(mod_var_l1_diff_sub, n.ahead=2, newdata=barry[n_ca-c(11:10),,drop=FALSE])
pred_l1_diff_2_12_nd <- predict(mod_var_l1_diff_sub, n.ahead=2, newdata=barry[n_ca-c(10:9),,drop=FALSE])
pred_l1_diff_1_12 <- predict(mod_var_l1_diff_sub, n.ahead=2)
all.equal(pred_l1_diff_1_12_nd, pred_l1_diff_1_12) ## minor: consistency in predict with/withotut newdata=dataset


pred_l1_diff_nd <- rbind(pred_l1_diff_0_12_nd, pred_l1_diff_1_12_nd, pred_l1_diff_2_12_nd)
all.equal(pred_l1_diff_nd[c(3,5),], as.M(pred_l1_diff_roll_12$pred[1:2,1:3]), check.attributes=FALSE) ## check 1-ahead
all.equal(pred_l1_diff_nd[c(2,4),], as.M(pred_l1_diff_roll_12$pred[11:12,1:3]), check.attributes=FALSE) ## check 2 ahead


