getIndic = function(lp,lpnew,Surv.rsp,Surv.rsp.new,times.auc=seq(10,1000,10),times.prederr=1:500,train.fit,train.fit.cph,tmax.train=365,tmax.test=365,TR,TE,plot.it=TRUE){
  try(attachNamespace("survival"),silent=TRUE)
  #on.exit(try(unloadNamespace("survival"),silent=TRUE))
  try(attachNamespace("risksetROC"),silent=TRUE)
  on.exit(try(unloadNamespace("risksetROC"),silent=TRUE))
  try(attachNamespace("survcomp"),silent=TRUE)
  on.exit(try(unloadNamespace("survcomp"),silent=TRUE))
  try(attachNamespace("rms"),silent=TRUE)
  on.exit(try(unloadNamespace("rms"),silent=TRUE))
  
  object <- NULL
  
  object$tmax.train <- tmax.train
  object$tmax.test <- tmax.test
  
  object$TR <- TR
  object$TE <- TE
  
  object$lp <- lp
  object$lpnew <- lpnew
  object$Surv.rsp <- Surv.rsp
  object$Surv.rsp.new <- Surv.rsp.new
  object$times.auc <- times.auc
  object$times.prederr <- times.prederr
  object$train.fit <- train.fit
  object$train.fit.cph <- train.fit.cph
  
  object$nulltrain.fit <- coxph(object$Surv.rsp~1)
  object$lp0 <- predict(object$nulltrain.fit)
  object$nulltest.fit <- coxph(object$Surv.rsp.new~1)
  object$lp0new <- predict(object$nulltest.fit)
  object$test.fit <- coxph(object$Surv.rsp.new~object$lpnew, iter.max=0, init=1)
  
  object$Erestrain <- predict(train.fit, type='expected')
  object$Erestest <- predict(train.fit, newdata=TE, type='expected')
  
  #predicted Martingale resid from coxph
  object$mtrainresid <- Surv.rsp[,"status"] - object$Erestrain
  object$mtestresid <- Surv.rsp.new[,"status"] - object$Erestest
  #then var of them
  object$var_mtrainresid <- var(object$mtrainresid)
  object$var_mtestresid <- var(object$mtestresid)
  
  #LRT test stat and/or p-value
  
  #R^2 Neglekerke
  # cox$loglik[1] loglik with initial values of coef
  # cox$loglik[2] loglik with final values of coef
  # cox$n number of observations
  object$logtraintest <- -2 * (object$nulltrain.fit$loglik - object$train.fit$loglik[2])
  object$logtesttest <- -2 * (object$nulltest.fit$loglik - object$test.fit$loglik[2])
  
  object$rval <- NULL
  # Cox and Snell
  object$rval$train$rsq.cs <- c(rsq.cs = 1 - exp(-object$logtraintest/object$train.fit$n), maxrsq.cs = 1 -exp(2 * object$train.fit$loglik[1]/object$train.fit$n))
  object$rval$test$rsq.cs <- c(rsq.cs = 1 - exp(-object$logtesttest/object$test.fit$n), maxrsq.cs = 1 -exp(2 * object$nulltest.fit$loglik/object$test.fit$n))
  # Nagelkerke
  object$rval$train$rsq.nagel <- c(rsq.nagel = unname(object$rval$train$rsq.cs[1]/object$rval$train$rsq.cs[2]), maxrsq.nagel = 1)
  object$rval$test$rsq.nagel <- c(rsq.nagel = unname(object$rval$test$rsq.cs[1]/object$rval$test$rsq.cs[2]), maxrsq.nagel = 1)
  
  
#  library(survAUC)
  # iAUC
  object$AUC_CD <- AUC.cd(object$Surv.rsp, object$Surv.rsp.new, object$lp, object$lpnew, object$times.auc)  #CoxModel
  object$AUC_hc <- AUC.hc(object$Surv.rsp, object$Surv.rsp.new, object$lpnew, object$times.auc)      #No model
  object$AUC_sh <- AUC.sh(object$Surv.rsp, object$Surv.rsp.new, object$lp, object$lpnew, object$times.auc)  #CoxModel
  object$AUC_Uno <- AUC.uno(object$Surv.rsp, object$Surv.rsp.new, object$lpnew, object$times.auc)    #No model
  
  #AUC_CD$iauc
  #AUC_hc$iauc
  #AUC_sh$iauc
  #AUC_Uno$iauc
  
  if(plot.it){
    layout(matrix(1:4,nrow=2,byrow=TRUE))
    plot(object$AUC_CD)
    abline(h = 0.5)
    plot(object$AUC_hc)
    abline(h = 0.5)
    plot(object$AUC_sh)
    abline(h = 0.5)
    plot(object$AUC_Uno)
    abline(h = 0.5)
  }
  
  #iAUC HZ train.set
#  library(risksetROC)
  ## first find the estimated survival probabilities at unique failure times
  object$surv.prob.train <- unique(survival::survfit(object$Surv.rsp~1)$surv)
  object$eta.train <- predict(object$train.fit)
  #model.score <- eta.train
  
  object$utimes.train <- unique( object$Surv.rsp[,"time"][ object$Surv.rsp[,"status"] == 1 ] )
  object$utimes.train <- object$utimes.train[ order(object$utimes.train) ]
  
  ## find AUC at unique failure times
  object$AUC_hz.train <- NULL
  object$AUC_hz.train$auc <- rep( NA, length(object$utimes.train) )
  for( j in 1:length(object$utimes.train) )
  {
    object$out.train <- CoxWeights(object$eta.train, object$Surv.rsp[,"time"], object$Surv.rsp[,"status"], object$utimes.train[j])
    object$AUC_hz.train$auc[j] <- object$out.train$AUC
  }
  ## integrated AUC to get concordance measure
  object$AUC_hz.train$times <- object$utimes.train
  object$AUC_hz.train$iauc <- IntegrateAUC( object$AUC_hz.train$auc, object$utimes.train, object$surv.prob.train, tmax=object$tmax.train)
  class(object$AUC_hz.train) <- "survAUC"
  if(plot.it){
    layout(matrix(1:4,nrow=2,byrow=TRUE))
    plot(object$AUC_hz.train)
    abline(h = 0.5)
  }
  
  #iAUC HZ test.set
#  library(risksetROC)
  ## first find the estimated survival probabilities at unique failure times
  object$surv.prob.test <- unique(survival::survfit(object$Surv.rsp.new~1)$surv)
  object$eta.test <- predict(object$test.fit)
  #model.score <- eta.test
  
  object$utimes.test <- unique( object$Surv.rsp.new[,"time"][ object$Surv.rsp.new[,"status"] == 1 ] )
  object$utimes.test <- object$utimes.test[ order(object$utimes.test) ]
  
  ## find AUC at unique failure times
  object$AUC_hz.test <- NULL
  object$AUC_hz.test$auc <- rep( NA, length(object$utimes.test) )
  for( j in 1:length(object$utimes.test) )
  {
    object$out.test <- CoxWeights( object$eta.test, object$Surv.rsp.new[,"time"], object$Surv.rsp.new[,"status"], object$utimes.test[j])
    object$AUC_hz.test$auc[j] <- object$out.test$AUC
  }
  ## integrated AUC to get concordance measure
  object$AUC_hz.test$times <- object$utimes.test
  object$AUC_hz.test$iauc <- NA
  if(length(object$utimes.test)>1){
    object$AUC_hz.test$iauc <- IntegrateAUC( object$AUC_hz.test$auc, object$utimes.test, object$surv.prob.test, tmax=object$tmax.test )
  } else {if(length(object$utimes.test)==1){object$AUC_hz.test$iauc<-object$AUC_hz.test$auc}}
  class(object$AUC_hz.test) <- "survAUC"
  if(plot.it){
    plot(object$AUC_hz.test)
    abline(h = 0.5)
  }
  
  
  ##time-dependent ROC curves
  ##time-dependent ROC curves
#  library(survcomp)
  ##train
  mytdroc.train <- NULL
  mytdroc.train <- NULL
  object$AUC_survivalROC.train <- NULL
  object$AUC_survivalROC.train$auc <- rep(NA,length(object$utimes.train))
  object$AUC_survivalROC.train$iauc <- NA
  object$AUC_survivalROC.train$times <- object$utimes.train
  class(object$AUC_survivalROC.train) <- "survAUC"
  train.cc.ix <- complete.cases(object$lp, object$Surv.rsp[,"time"], object$Surv.rsp[,"status"], NULL)
  train.surv.event.cc.ix <- object$Surv.rsp.new[,"status"][train.cc.ix]
  if (all(sort(unique(train.surv.event.cc.ix)) == c(0, 1))) {
    for(i in 1:length(object$utimes.train)) {
      rr.train <- tdrocc(x=object$lp, surv.time=object$Surv.rsp[,"time"], surv.event=object$Surv.rsp[,"status"], time=object$utimes.train[i], na.rm=TRUE, verbose=FALSE)
      mytdroc.train <- c(mytdroc.train, list(rr.train))
    }
    object$AUC_survivalROC.train$auc <- unlist(lapply(mytdroc.train, function(x) { return(x$AUC) }))
    cc.ix.train <- complete.cases(object$AUC_survivalROC.train$auc)
    auc.survivalROC.train.cc <- object$AUC_survivalROC.train$auc[cc.ix.train]
    time.train.cc <- object$utimes.train[cc.ix.train]
    if(length(time.train.cc)>0){
      diffs.train.cc <- c(time.train.cc[1], time.train.cc[2:length(time.train.cc)] - time.train.cc[1:(length(time.train.cc) - 1)])
      object$AUC_survivalROC.train$iauc <- sum(diffs.train.cc * auc.survivalROC.train.cc)/max(time.train.cc)
      if(plot.it){
        plot(object$AUC_survivalROC.train)
        abline(h = 0.5)
      }
    }
  }
  
  
#  library(survcomp)
  ##test
  mytdroc.test <- NULL
  object$AUC_survivalROC.test <- NULL
  object$AUC_survivalROC.test$auc <- rep(NA,length(object$utimes.test))
  object$AUC_survivalROC.test$iauc <- NA
  object$AUC_survivalROC.test$times <- object$utimes.test
  class(object$AUC_survivalROC.test) <- "survAUC"
  test.cc.ix <- complete.cases(object$lpnew, object$Surv.rsp.new[,"time"], object$Surv.rsp.new[,"status"], NULL)
  test.surv.event.cc.ix <- object$Surv.rsp.new[,"status"][test.cc.ix]
  if (all(sort(unique(test.surv.event.cc.ix)) == c(0, 1))) {
    for(i in 1:length(object$utimes.test)) {
      rr.test <- tdrocc(x=object$lpnew, surv.time=object$Surv.rsp.new[,"time"], surv.event=object$Surv.rsp.new[,"status"], time=object$utimes.test[i], na.rm=TRUE, verbose=FALSE)
      mytdroc.test <- c(mytdroc.test, list(rr.test))
    }
    object$AUC_survivalROC.test$auc <- unlist(lapply(mytdroc.test, function(x) { return(x$AUC) }))
    cc.ix.test <- complete.cases(object$AUC_survivalROC.test$auc)
    auc.survivalROC.test.cc <- object$AUC_survivalROC.test$auc[cc.ix.test]
    time.test.cc <- object$utimes.test[cc.ix.test]
    if(length(time.test.cc)>0){
      diffs.test.cc <- c(time.test.cc[1], time.test.cc[2:length(time.test.cc)] - time.test.cc[1:(length(time.test.cc) - 1)])
      object$AUC_survivalROC.test$iauc <- sum(diffs.test.cc * auc.survivalROC.test.cc)/max(time.test.cc)
      if(plot.it){
        plot(object$AUC_survivalROC.test)
        abline(h = 0.5)
      }
    }
  }
  
  
  object$HarrelC <- BeggC(object$Surv.rsp, object$Surv.rsp.new, object$lp, object$lpnew)  #CoxModel C-statistic by Begg et al.
  object$GonenHellerCI <- GHCI(object$lpnew)  #CoxModel Gonen and Heller?s Concordance Index for Cox models
  
  object$rval$train$rsq.OXS <- OXS(object$Surv.rsp, object$lp, object$lp0)
  object$rval$train$rsq.Nagelk <- Nagelk(object$Surv.rsp, object$lp, object$lp0)
  object$rval$train$rsq.XO <- XO(object$Surv.rsp, object$lp, object$lp0)
  
  object$rval$test$rsq.OXS <- OXS(object$Surv.rsp.new, object$lpnew, object$lp0new)
  object$rval$test$rsq.Nagelk <- Nagelk(object$Surv.rsp.new, object$lpnew, object$lp0new)
  object$rval$test$rsq.XO <- XO(object$Surv.rsp.new, object$lpnew, object$lp0new)
  
  #Surv.rsp A Surv(.,.) object containing to the outcome of the test data.
  #lp The vector of predictors.
  #lp0 The vector of predictors obtained from the covariate-free null model.
  
  #prederr #ierror for iBrier
  object$prederr <- NULL
  object$times.prederr <- times.prederr
  object$prederr$brier.unw <- predErr(object$Surv.rsp, object$Surv.rsp.new, object$lp, object$lpnew, object$times.prederr, type = "brier", int.type = "unweighted")
  object$prederr$robust.unw <- predErr(object$Surv.rsp, object$Surv.rsp.new, object$lp, object$lpnew, object$times.prederr, type = "robust", int.type = "unweighted")
  object$prederr$brier.w <- predErr(object$Surv.rsp, object$Surv.rsp.new, object$lp, object$lpnew, object$times.prederr, type = "brier", int.type = "weighted")
  object$prederr$robust.w <- predErr(object$Surv.rsp, object$Surv.rsp.new, object$lp, object$lpnew, object$times.prederr, type = "robust", int.type = "weighted")
  object$prederr$brier0.unw <- predErr(object$Surv.rsp, object$Surv.rsp.new, object$lp0, object$lp0new, object$times.prederr, type = "brier", int.type = "unweighted")
  object$prederr$robust0.unw <- predErr(object$Surv.rsp, object$Surv.rsp.new, object$lp0, object$lp0new, object$times.prederr, type = "robust", int.type = "unweighted")
  object$prederr$brier0.w <- predErr(object$Surv.rsp, object$Surv.rsp.new, object$lp0, object$lp0new, object$times.prederr, type = "brier", int.type = "weighted")
  object$prederr$robust0.w <- predErr(object$Surv.rsp, object$Surv.rsp.new, object$lp0, object$lp0new, object$times.prederr, type = "robust", int.type = "weighted")
  
  if(plot.it){
    layout(matrix(1:4,nrow=2))
    plot(object$prederr$brier.unw)
    abline(h = 0.25)
    plot(object$prederr$robust.unw)
    abline(h = 0.25)
    plot(object$prederr$brier.w)
    abline(h = 0.25)
    plot(object$prederr$robust.w)
    abline(h = 0.25)
  }
  
  
  #integrated R2 from BS        max(T_i)=temps (censur? ou non) le plus grand=max(time). Fonction continue lin?aire par morceaux.
  # R2_{BS}(t)=1-BS(t)/BS_0(t)
  # iR2BS=1/max(T_i)int_0^{max(T_i)}R2_{BS}(t)dt.
  object$rval$test$R2.bs.unw <- 1-object$prederr$brier.unw$error/object$prederr$brier0.unw$error
  object$rval$test$R2.bs.unw[object$prederr$brier.unw$error==0] <- 0
  object$rval$test$R2.bs.unw[object$rval$test$R2.bs.unw<=0] <- 0
  object$rval$test$R2.rs.unw <- 1-object$prederr$robust.unw$error/object$prederr$robust0.unw$error
  object$rval$test$R2.rs.unw[object$prederr$robust.unw$error==0] <- 0
  object$rval$test$R2.rs.unw[object$rval$test$R2.rs.unw<=0] <- 0
  object$rval$test$R2.bs.w <- 1-object$prederr$brier.w$error/object$prederr$brier0.w$error
  object$rval$test$R2.bs.w[object$prederr$brier.w$error==0] <- 0
  object$rval$test$R2.bs.w[object$rval$test$R2.bs.w<=0] <- 0
  object$rval$test$R2.rs.w <- 1-object$prederr$robust.w$error/object$prederr$robust0.w$error
  object$rval$test$R2.rs.w[object$prederr$robust.w$error==0] <- 0
  object$rval$test$R2.rs.w[object$rval$test$R2.rs.w<=0] <- 0
  
  object$rval$test$iRbs.unw <- sum(object$rval$test$R2.bs.unw[-1]*diff(object$prederr$brier.unw$time))/max(object$prederr$brier.unw$time)
  object$rval$test$iRrs.unw <- sum(object$rval$test$R2.rs.unw[-1]*diff(object$prederr$robust.unw$time))/max(object$prederr$robust.unw$time)
  object$rval$test$iRbs.w <- sum(object$rval$test$R2.bs.w[-1]*diff(object$prederr$brier.w$time))/max(object$prederr$brier.w$time)
  object$rval$test$iRrs.w <- sum(object$rval$test$R2.rs.w[-1]*diff(object$prederr$robust.w$time))/max(object$prederr$robust.w$time)
  
  #UnoC C-statistic by Uno et al.
  object$Cstat <- UnoC(object$Surv.rsp, object$Surv.rsp.new, object$lpnew) #no model
  
  #schemper Distance-based estimator of survival predictive accuracy proposed by Schemper and Henderson
  #Schemper and Henderson's estimator of the absolute deviation between survival functions
#  library(rms)
  object$Schemper <- schemper(object$train.fit.cph, object$TR, object$TE)
  return(invisible(object))
}











getIndicCV = function(lp,lpnew,Surv.rsp,Surv.rsp.new,times.auc=seq(10,1000,10),times.prederr=1:500,train.fit,plot.it=FALSE,tmax.train=max(Surv.rsp[,"time"][ object$Surv.rsp[,"status"] == 1 ]),tmax.test=max(Surv.rsp.new[,"time"][ object$Surv.rsp.new[,"status"] == 1 ])){
  try(attachNamespace("survival"),silent=TRUE)
  #on.exit(try(unloadNamespace("survival"),silent=TRUE))
  try(attachNamespace("risksetROC"),silent=TRUE)
  on.exit(try(unloadNamespace("risksetROC"),silent=TRUE))
  try(attachNamespace("survcomp"),silent=TRUE)
  on.exit(try(unloadNamespace("survcomp"),silent=TRUE))
  #  library(survAUC)
  object <- NULL
  
  object$lp <- lp
  object$lpnew <- lpnew
  object$Surv.rsp <- Surv.rsp
  object$Surv.rsp.new <- Surv.rsp.new
  object$times.auc <- times.auc
  object$train.fit <- train.fit
  object$test.fit <- coxph(object$Surv.rsp.new~object$lpnew, iter.max=0, init=1)
  
  object$nulltrain.fit <- coxph(object$Surv.rsp~1)
  object$lp0 <- predict(object$nulltrain.fit)
  object$nulltest.fit <- coxph(object$Surv.rsp.new~1)
  object$lp0new <- predict(object$nulltest.fit)
  
  object$tmax.train <- tmax.train
  object$tmax.test <- tmax.test
  
  # iAUC
  object$AUC_CD <- AUC.cd(object$Surv.rsp, object$Surv.rsp.new, object$lp, object$lpnew, object$times.auc)  #CoxModel
  object$AUC_hc <- AUC.hc(object$Surv.rsp, object$Surv.rsp.new, object$lpnew, object$times.auc)      #No model
  object$AUC_sh <- AUC.sh(object$Surv.rsp, object$Surv.rsp.new, object$lp, object$lpnew, object$times.auc)  #CoxModel
  object$AUC_Uno <- list(auc=rep(0,length(object$times.auc)),times=times.auc,iauc=0)
  #class(object$AUC_Uno) <- "survAUC"
  #if(var(object$lpnew)>1e-8){try(
  object$AUC_Uno <- AUC.uno(object$Surv.rsp, object$Surv.rsp.new, object$lpnew, object$times.auc)#  )}   #No model
  
  #AUC_CD$iauc
  #AUC_hc$iauc
  #AUC_sh$iauc
  #AUC_Uno$iauc
  
  if(plot.it){
    layout(matrix(1:4,nrow=2,byrow=TRUE))
    plot(object$AUC_CD)
    abline(h = 0.5)
    plot(object$AUC_hc)
    abline(h = 0.5)
    plot(object$AUC_sh)
    abline(h = 0.5)
    plot(object$AUC_Uno)
    abline(h = 0.5)
  }
  
  #iAUC HZ train.set
  #  library(risksetROC)
  ## first find the estimated survival probabilities at unique failure times
  object$surv.prob.train <- unique(survival::survfit(object$Surv.rsp~1)$surv)
  object$eta.train <- predict(object$train.fit)
  #model.score <- eta.train
  
  object$utimes.train <- unique( object$Surv.rsp[,"time"][ object$Surv.rsp[,"status"] == 1 ] )
  object$utimes.train <- object$utimes.train[ order(object$utimes.train) ]
  
  ## find AUC at unique failure times
  object$AUC_hz.train <- NULL
  object$AUC_hz.train$auc <- rep( NA, length(object$utimes.train) )
  for( j in 1:length(object$utimes.train) )
  {
    object$out.train <- CoxWeights(object$eta.train, object$Surv.rsp[,"time"], object$Surv.rsp[,"status"], object$utimes.train[j])
    object$AUC_hz.train$auc[j] <- object$out.train$AUC
  }
  ## integrated AUC to get concordance measure
  object$AUC_hz.train$times <- object$utimes.train
  object$AUC_hz.train$iauc <- IntegrateAUC( object$AUC_hz.train$auc, object$utimes.train, object$surv.prob.train, tmax=object$tmax.train)
  class(object$AUC_hz.train) <- "survAUC"
  if(plot.it){
    layout(matrix(1:4,nrow=2))
    plot(object$AUC_hz.train)
    abline(h = 0.5)
  }
  
  #iAUC HZ test.set
  #  library(risksetROC)
  ## first find the estimated survival probabilities at unique failure times
  object$surv.prob.test <- unique(survival::survfit(object$Surv.rsp.new~1)$surv)
  object$eta.test <- predict(object$test.fit)
  #model.score <- eta.test
  
  object$utimes.test <- unique( object$Surv.rsp.new[,"time"][ object$Surv.rsp.new[,"status"] == 1 ] )
  object$utimes.test <- object$utimes.test[ order(object$utimes.test) ]
  
  ## find AUC at unique failure times
  object$AUC_hz.test <- NULL
  object$AUC_hz.test$auc <- rep( NA, length(object$utimes.test) )
  for( j in 1:length(object$utimes.test) )
  {
    object$out.test <- CoxWeights( object$eta.test, object$Surv.rsp.new[,"time"], object$Surv.rsp.new[,"status"], object$utimes.test[j])
    object$AUC_hz.test$auc[j] <- object$out.test$AUC
  }
  ## integrated AUC to get concordance measure
  object$AUC_hz.test$times <- object$utimes.test
  object$AUC_hz.test$iauc <- NA
  if(length(object$utimes.test)>1){
    object$AUC_hz.test$iauc <- IntegrateAUC( object$AUC_hz.test$auc, object$utimes.test, object$surv.prob.test, tmax=object$tmax.test )
  } else {if(length(object$utimes.test)==1){object$AUC_hz.test$iauc<-object$AUC_hz.test$auc}}
  class(object$AUC_hz.test) <- "survAUC"
  if(plot.it){
    plot(object$AUC_hz.test)
    abline(h = 0.5)
  }
  
  ##time-dependent ROC curves
  #  library(survcomp)
  ##train
  mytdroc.train <- NULL
  mytdroc.train <- NULL
  object$AUC_survivalROC.train <- NULL
  object$AUC_survivalROC.train$auc <- rep(NA,length(object$utimes.train))
  object$AUC_survivalROC.train$iauc <- NA
  object$AUC_survivalROC.train$times <- object$utimes.train
  class(object$AUC_survivalROC.train) <- "survAUC"
  train.cc.ix <- complete.cases(object$lp, object$Surv.rsp[,"time"], object$Surv.rsp[,"status"], NULL)
  train.surv.event.cc.ix <- object$Surv.rsp.new[,"status"][train.cc.ix]
  if (all(sort(unique(train.surv.event.cc.ix)) == c(0, 1))) {
    for(i in 1:length(object$utimes.train)) {
      rr.train <- tdrocc(x=object$lp, surv.time=object$Surv.rsp[,"time"], surv.event=object$Surv.rsp[,"status"], time=object$utimes.train[i], na.rm=TRUE, verbose=FALSE)
      mytdroc.train <- c(mytdroc.train, list(rr.train))
    }
    object$AUC_survivalROC.train$auc <- unlist(lapply(mytdroc.train, function(x) { return(x$AUC) }))
    cc.ix.train <- complete.cases(object$AUC_survivalROC.train$auc)
    auc.survivalROC.train.cc <- object$AUC_survivalROC.train$auc[cc.ix.train]
    time.train.cc <- object$utimes.train[cc.ix.train]
    if(length(time.train.cc)>0){
      diffs.train.cc <- c(time.train.cc[1], time.train.cc[2:length(time.train.cc)] - time.train.cc[1:(length(time.train.cc) - 1)])
      object$AUC_survivalROC.train$iauc <- sum(diffs.train.cc * auc.survivalROC.train.cc)/max(time.train.cc)
      if(plot.it){
        plot(object$AUC_survivalROC.train)
        abline(h = 0.5)
      }
    }
  }
  
  
  #  library(survcomp)
  ##test
  mytdroc.test <- NULL
  object$AUC_survivalROC.test <- NULL
  object$AUC_survivalROC.test$auc <- rep(NA,length(object$utimes.test))
  object$AUC_survivalROC.test$iauc <- NA
  object$AUC_survivalROC.test$times <- object$utimes.test
  class(object$AUC_survivalROC.test) <- "survAUC"
  test.cc.ix <- complete.cases(object$lpnew, object$Surv.rsp.new[,"time"], object$Surv.rsp.new[,"status"], NULL)
  test.surv.event.cc.ix <- object$Surv.rsp.new[,"status"][test.cc.ix]
  if (all(sort(unique(test.surv.event.cc.ix)) == c(0, 1))) {
    for(i in 1:length(object$utimes.test)) {
      rr.test <- tdrocc(x=object$lpnew, surv.time=object$Surv.rsp.new[,"time"], surv.event=object$Surv.rsp.new[,"status"], time=object$utimes.test[i], na.rm=TRUE, verbose=FALSE)
      mytdroc.test <- c(mytdroc.test, list(rr.test))
    }
    object$AUC_survivalROC.test$auc <- unlist(lapply(mytdroc.test, function(x) { return(x$AUC) }))
    cc.ix.test <- complete.cases(object$AUC_survivalROC.test$auc)
    auc.survivalROC.test.cc <- object$AUC_survivalROC.test$auc[cc.ix.test]
    time.test.cc <- object$utimes.test[cc.ix.test]
    if(length(time.test.cc)>0){
      diffs.test.cc <- c(time.test.cc[1], time.test.cc[2:length(time.test.cc)] - time.test.cc[1:(length(time.test.cc) - 1)])
      object$AUC_survivalROC.test$iauc <- sum(diffs.test.cc * auc.survivalROC.test.cc)/max(time.test.cc)
      if(plot.it){
        plot(object$AUC_survivalROC.test)
        abline(h = 0.5)
      }
    }
  }
  #Surv.rsp A Surv(.,.) object containing to the outcome of the test data.
  #lp The vector of predictors.
  #lp0 The vector of predictors obtained from the covariate-free null model.
  
  #prederr #ierror for iBrier
  object$prederr <- NULL
  object$times.prederr <- times.prederr[times.prederr>1]
  object$prederr$brier.unw <- predErr(object$Surv.rsp, object$Surv.rsp.new, object$lp, object$lpnew, object$times.prederr, type = "brier", int.type = "unweighted")
  object$prederr$robust.unw <- predErr(object$Surv.rsp, object$Surv.rsp.new, object$lp, object$lpnew, object$times.prederr, type = "robust", int.type = "unweighted")
  object$prederr$brier.w <- predErr(object$Surv.rsp, object$Surv.rsp.new, object$lp, object$lpnew, object$times.prederr, type = "brier", int.type = "weighted")
  object$prederr$robust.w <- predErr(object$Surv.rsp, object$Surv.rsp.new, object$lp, object$lpnew, object$times.prederr, type = "robust", int.type = "weighted")
  
  if(plot.it){
    layout(matrix(1:4,nrow=2))
    plot(object$prederr$brier.unw)
    abline(h = 0.25)
    plot(object$prederr$robust.unw)
    abline(h = 0.25)
    plot(object$prederr$brier.w)
    abline(h = 0.25)
    plot(object$prederr$robust.w)
    abline(h = 0.25)
  }
  
  return(object)
}



getIndicCViAUCSH = function(lp,lpnew,Surv.rsp,Surv.rsp.new,times.auc=seq(10,1000,10),times.prederr=1:500,train.fit,plot.it=FALSE,tmax.train=max(Surv.rsp[,"time"][ object$Surv.rsp[,"status"] == 1 ]),tmax.test=max(Surv.rsp.new[,"time"][ object$Surv.rsp.new[,"status"] == 1 ])){
  try(attachNamespace("survival"),silent=TRUE)
  #on.exit(try(unloadNamespace("survival"),silent=TRUE))
  try(attachNamespace("risksetROC"),silent=TRUE)
  on.exit(try(unloadNamespace("risksetROC"),silent=TRUE))
  try(attachNamespace("survcomp"),silent=TRUE)
  on.exit(try(unloadNamespace("survcomp"),silent=TRUE))
  #  library(survAUC)
  object <- NULL
  
  object$lp <- lp
  object$lpnew <- lpnew
  object$Surv.rsp <- Surv.rsp
  object$Surv.rsp.new <- Surv.rsp.new
  object$times.auc <- times.auc
  object$train.fit <- train.fit
  object$test.fit <- coxph(object$Surv.rsp.new~object$lpnew, iter.max=0, init=1)
  
  object$nulltrain.fit <- coxph(object$Surv.rsp~1)
  object$lp0 <- predict(object$nulltrain.fit)
  object$nulltest.fit <- coxph(object$Surv.rsp.new~1)
  object$lp0new <- predict(object$nulltest.fit)
  
  object$tmax.train <- tmax.train
  object$tmax.test <- tmax.test
  
  # iAUC
  object$AUC_sh <- AUC.sh(object$Surv.rsp, object$Surv.rsp.new, object$lp, object$lpnew, object$times.auc)  #CoxModel
  
  if(plot.it){
    plot(object$AUC_sh)
    abline(h = 0.5)
  }
  
  return(object)
}



getIndicCViAUCSurvROCTest = function(lp,lpnew,Surv.rsp,Surv.rsp.new,times.auc=seq(10,1000,10),times.prederr=1:500,train.fit,plot.it=FALSE,tmax.train=max(Surv.rsp[,"time"][ object$Surv.rsp[,"status"] == 1 ]),tmax.test=max(Surv.rsp.new[,"time"][ object$Surv.rsp.new[,"status"] == 1 ])){
  try(attachNamespace("survival"),silent=TRUE)
  #on.exit(try(unloadNamespace("survival"),silent=TRUE))
  try(attachNamespace("survcomp"),silent=TRUE)
  on.exit(try(unloadNamespace("survcomp"),silent=TRUE))
  object <- NULL
  
  object$lp <- lp
  object$lpnew <- lpnew
  object$Surv.rsp <- Surv.rsp
  object$Surv.rsp.new <- Surv.rsp.new
  object$times.auc <- times.auc
  object$train.fit <- train.fit
  object$test.fit <- coxph(object$Surv.rsp.new~object$lpnew, iter.max=0, init=1)
  
  object$nulltrain.fit <- coxph(object$Surv.rsp~1)
  object$lp0 <- predict(object$nulltrain.fit)
  object$nulltest.fit <- coxph(object$Surv.rsp.new~1)
  object$lp0new <- predict(object$nulltest.fit)
  
  object$tmax.train <- tmax.train
  object$tmax.test <- tmax.test
  
  
  object$utimes.test <- unique( object$Surv.rsp.new[,"time"][ object$Surv.rsp.new[,"status"] == 1 ] )
  object$utimes.test <- object$utimes.test[ order(object$utimes.test) ]
      
  #  library(survcomp)
  ##test
  mytdroc.test <- NULL
  object$AUC_survivalROC.test <- NULL
  object$AUC_survivalROC.test$auc <- rep(NA,length(object$utimes.test))
  object$AUC_survivalROC.test$iauc <- NA
  object$AUC_survivalROC.test$times <- object$utimes.test
  class(object$AUC_survivalROC.test) <- "survAUC"
  test.cc.ix <- complete.cases(object$lpnew, object$Surv.rsp.new[,"time"], object$Surv.rsp.new[,"status"], NULL)
  test.surv.event.cc.ix <- object$Surv.rsp.new[,"status"][test.cc.ix]
  if (all(sort(unique(test.surv.event.cc.ix)) == c(0, 1))) {
    for(i in 1:length(object$utimes.test)) {
      rr.test <- tdrocc(x=object$lpnew, surv.time=object$Surv.rsp.new[,"time"], surv.event=object$Surv.rsp.new[,"status"], time=object$utimes.test[i], na.rm=TRUE, verbose=FALSE)
      mytdroc.test <- c(mytdroc.test, list(rr.test))
    }
    object$AUC_survivalROC.test$auc <- unlist(lapply(mytdroc.test, function(x) { return(x$AUC) }))
    cc.ix.test <- complete.cases(object$AUC_survivalROC.test$auc)
    auc.survivalROC.test.cc <- object$AUC_survivalROC.test$auc[cc.ix.test]
    time.test.cc <- object$utimes.test[cc.ix.test]
    if(length(time.test.cc)>0){
      diffs.test.cc <- c(time.test.cc[1], time.test.cc[2:length(time.test.cc)] - time.test.cc[1:(length(time.test.cc) - 1)])
      object$AUC_survivalROC.test$iauc <- sum(diffs.test.cc * auc.survivalROC.test.cc)/max(time.test.cc)
      if(plot.it){
        plot(object$AUC_survivalROC.test)
        abline(h = 0.5)
      }
    }
  }
  
  return(object)
}



logplik = function (x, time, status, b, method = c("breslow", "efron"),
                    return.all = FALSE)
{
  method <- match.arg(method)
  n <- length(time)
  o <- order(status, decreasing = T)
  oo <- o[order(time[o])]
  time <- time[oo]
  status <- status[oo]
  rept <- rep(0, n)
  for (i in 1:n) rept[i] <- sum(time[i:n] == time[i] & status[i:n] ==
                                  1)
  complete <- which(status == 1)
  nnc <- length(complete)
  if (nnc == 0) {
    stop("No complete observation. Failed to compute partial likelihood.")
  }
  dmat <- matrix(0, n, nnc)
  for (i in 1:nnc) {
    dmat[time >= time[complete[i]], i] <- 1
    if (method == "efron") {
      if (rept[complete[i]] > 0) {
        tie <- time == time[complete[i]] & status ==
          1
        di <- max(rept[tie])
        dmat[tie, i] <- dmat[tie, i] - (di - rept[complete[i]])/di
      }
    }
  }
  eta <- x %*% b
  eeta <- exp(eta)
  k <- ncol(eta)
  loglik <- rep(0, k)
  for (i in 1:k) {
    w <- dmat * eeta[oo, i]
    wsum <- apply(w, 2, sum)
    loglik[i] <- sum(eta[oo, i][status == 1]) - sum(log(wsum))
  }
  if (return.all) {
    return(list(loglik = loglik, w = scale(w, F, wsum), eta = eta,
                dmat = dmat, oo = oo))
  }
  else {
    return(loglik)
  }
}

getmin2 = function (lambda, cvm, cvsd)
{
  cvmin = min(cvm, na.rm = TRUE)
  idmin = cvm <= cvmin
  lambda.min = max(lambda[idmin], na.rm = TRUE)
  idminl = match(lambda.min, lambda)
  semin = (cvm + cvsd)[idminl]
  idmin2 = cvm >= semin
  #    if(lambda.min==-Inf){
  #    lambda.1se = -Inf
  #    } else {
  idmin2[idminl:length(idmin2)] = FALSE 
  lambda.1se = min(c(lambda[idmin2],min(lambda)), na.rm = TRUE)
  #    }
  list(lambda.min = lambda.min, lambda.1se = lambda.1se)
}


correctp.cox=function (x, y, eta, K, kappa, select, fit, verbose=FALSE) 
{
  force(K)
  if (min(eta) < 0 | max(eta) >= 1) {
    if (max(eta) == 1) {
      stop("eta should be strictly less than 1!")
    }
    if (length(eta) == 1) {
      stop("eta should be between 0 and 1!")
    }
    else {
      stop("eta should be between 0 and 1! \n  Choose appropriate range of eta!")
    }
  }
  if (max(K) > ncol(x)) {
    stop("K cannot exceed the number of predictors! Pick up smaller K!")
  }
  if (max(K) >= nrow(x)) {
    stop("K cannot exceed the sample size! Pick up smaller K!")
  }
  if (min(K) <= 0 | !all(K%%1 == 0)) {
    if (length(K) == 1) {
      stop("K should be a positive integer!")
    }
    else {
      stop("K should be a positive integer! \n  Choose appropriate range of K!")
    }
  }
  if (kappa > 0.5 | kappa < 0) {
    if(verbose){cat("kappa should be between 0 and 0.5! kappa=0.5 is used. \n\n")}
    kappa <- 0.5
  }
  if (select != "pls2" & select != "simpls") {
    if(verbose){cat("Invalid PLS algorithm for variable selection.\n")}
    if(verbose){cat("pls2 algorithm is used. \n\n")}
    select <- "pls2"
  }
  fits <- c("regression", "canonical", "invariant", "classic")
  if (!any(fit == fits)) {
    if(verbose){cat("Invalid PLS algorithm for model fitting\n")}
    if(verbose){cat("regression algorithm is used. \n\n")}
    fit <- "regression"
  }
  list(K = K, eta = eta, kappa = kappa, select = select, fit = fit)
}

spls.cox=function (x, y, K, eta, kappa = 0.5, select = "pls2", fit = "regression", 
                   scale.x = TRUE, scale.y = FALSE, eps = 1e-04, maxstep = 100, 
                   trace = FALSE) 
{
  force(K)
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  ip <- c(1:p)
  y <- as.matrix(y)
  q <- ncol(y)
  one <- matrix(1, 1, n)
  mu <- one %*% y/n
  y <- scale(y, drop(mu), FALSE)
  meanx <- drop(one %*% x)/n
  x <- scale(x, meanx, FALSE)
  if (scale.x) {
    normx <- sqrt(drop(one %*% (x^2))/(n - 1))
    if (any(normx < .Machine$double.eps)) {
      stop("Some of the columns of the predictor matrix have zero variance.")
    }
    x <- scale(x, FALSE, normx)
  }
  else {
    normx <- rep(1, p)
  }
  if (scale.y) {
    normy <- sqrt(drop(one %*% (y^2))/(n - 1))
    if (any(normy < .Machine$double.eps)) {
      stop("Some of the columns of the response matrix have zero variance.")
    }
    y <- scale(y, FALSE, normy)
  }
  else {
    normy <- rep(1, q)
  }
  betahat <- matrix(0, p, q)
  betamat <- list()
  x1 <- x
  y1 <- y
  type <- correctp.cox(x, y, eta, K, kappa, select, fit)
  eta <- type$eta
  K <- type$K
  kappa <- type$kappa
  select <- type$select
  fit <- type$fit
  if (is.null(colnames(x))) {
    xnames <- c(1:p)
  }
  else {
    xnames <- colnames(x)
  }
  new2As <- list()
  if (trace) {
    cat("The variables that join the set of selected variables at each step:\n")
  }
  for (k in 1:K) {
    Z <- t(x1) %*% y1
    what <- spls.dv(Z, eta, kappa, eps, maxstep)
    A <- unique(ip[what != 0 | betahat[, 1] != 0])
    new2A <- ip[what != 0 & betahat[, 1] == 0]
    xA <- x[, A, drop = FALSE]
    plsfit <- pls.cox(X=xA, Y=y, ncomp = min(k, length(A)), mode = fit, 
                      scale.X = FALSE, scale.Y=FALSE)
    predplsfit <- predict.pls.cox(plsfit,newdata=xA,scale.X = FALSE, scale.Y=FALSE)
    betahat <- matrix(0, p, q)
    betahat[A, ] <- matrix(predplsfit$B.hat[,,plsfit$ncomp], length(A), q)
    betamat[[k]] <- betahat
    #        pj <- plsfit$projection
    if (select == "pls2") {
      y1 <- y - predplsfit$predict[,,plsfit$ncomp]
    }
    #        if (select == "simpls") {
    #            pw <- pj %*% solve(t(pj) %*% pj) %*% t(pj)
    #            x1 <- x
    #            x1[, A] <- x[, A, drop = FALSE] - x[, A, drop = FALSE] %*% 
    #                pw
    #        }
    new2As[[k]] <- new2A
    if (trace) {
      if (length(new2A) <= 10) {
        cat(paste("- ", k, "th step (K=", k, "):\n", 
                  sep = ""))
        cat(xnames[new2A])
        cat("\n")
      }
      else {
        cat(paste("- ", k, "th step (K=", k, "):\n", 
                  sep = ""))
        nlines <- ceiling(length(new2A)/10)
        for (i in 0:(nlines - 2)) {
          cat(xnames[new2A[(10 * i + 1):(10 * (i + 1))]])
          cat("\n")
        }
        cat(xnames[new2A[(10 * (nlines - 1) + 1):length(new2A)]])
        cat("\n")
      }
    }
  }
  if (!is.null(colnames(x))) {
    rownames(betahat) <- colnames(x)
  }
  if (q > 1 & !is.null(colnames(y))) {
    colnames(betahat) <- colnames(y)
  }
  object <- list(x = x, y = y, betahat = betahat, A = A, betamat = betamat, 
                 new2As = new2As, mu = mu, meanx = meanx, normx = normx, 
                 normy = normy, eta = eta, K = K, kappa = kappa, select = select, 
                 fit = fit, projection = NA, plsmod=plsfit)
  class(object) <- "spls"
  object
}

ust=function (b, eta) 
{
  b.ust <- matrix(0, length(b), 1)
  if (eta < 1) {
    valb <- abs(b) - eta * max(abs(b))
    b.ust[valb >= 0] <- valb[valb >= 0] * (sign(b))[valb >= 
                                                      0]
  }
  return(b.ust)
}

spls.dv <- function (Z, eta, kappa, eps, maxstep) 
{
  p <- nrow(Z)
  q <- ncol(Z)
  Znorm1 <- median(abs(Z))
  Z <- Z/Znorm1
  if (q == 1) {
    c <- ust(Z, eta)
  }
  if (q > 1) {
    M <- Z %*% t(Z)
    dis <- 10
    i <- 1
    if (kappa == 0.5) {
      c <- matrix(10, p, 1)
      c.old <- c
      while (dis > eps & i <= maxstep) {
        mcsvd <- svd(M %*% c)
        a <- mcsvd$u %*% t(mcsvd$v)
        c <- ust(M %*% a, eta)
        dis <- max(abs(c - c.old))
        c.old <- c
        i <- i + 1
      }
    }
    if (kappa > 0 & kappa < 0.5) {
      kappa2 <- (1 - kappa)/(1 - 2 * kappa)
      c <- matrix(10, p, 1)
      c.old <- c
      h <- function(lambda) {
        alpha <- solve(M + lambda * diag(p)) %*% M %*% 
          c
        obj <- t(alpha) %*% alpha - 1/kappa2^2
        return(obj)
      }
      if (h(eps) * h(1e+30) > 0) {
        while (h(eps) <= 1e+05) {
          M <- 2 * M
          c <- 2 * c
        }
      }
      while (dis > eps & i <= maxstep) {
        if (h(eps) * h(1e+30) > 0) {
          while (h(eps) <= 1e+05) {
            M <- 2 * M
            c <- 2 * c
          }
        }
        lambdas <- uniroot(h, c(eps, 1e+30))$root
        a <- kappa2 * solve(M + lambdas * diag(p)) %*% 
          M %*% c
        c <- ust(M %*% a, eta)
        dis <- max(abs(c - c.old))
        c.old <- c
        i <- i + 1
      }
    }
  }
  return(c)
}

pls.cox=function(X, Y, ncomp = 2, mode = c("regression", "canonical", "invariant", "classic"), max.iter = 500, tol = 1e-06, scale.X=TRUE, scale.Y=TRUE, ...) 
{
  force(ncomp)
  if (length(dim(X)) != 2) 
    stop("'X' must be a numeric matrix.")
  X = as.matrix(X)
  Y = as.matrix(Y)
  if (!is.numeric(X) || !is.numeric(Y)) 
    stop("'X' and/or 'Y' must be a numeric matrix.")
  n = nrow(X)
  q = ncol(Y)
  if ((n != nrow(Y))) 
    stop("unequal number of rows in 'X' and 'Y'.")
  if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0) 
    stop("invalid number of variates, 'ncomp'.")
  nzv = mixOmics::nearZeroVar(X, ...)
  if (length(nzv$Position > 0)) {
    warning("Zero- or near-zero variance predictors. \n  Reset predictors matrix to not near-zero variance predictors.\n  See $nzv for problematic predictors.")
    X = X[, -nzv$Position]
  }
  p = ncol(X)
  ncomp = round(ncomp)
  if (ncomp > p) {
    warning("Reset maximum number of variates 'ncomp' to ncol(X) = ", 
            p, ".")
    ncomp = p
  }
  mode = match.arg(mode)
  X.names = dimnames(X)[[2]]
  if (is.null(X.names)) 
    X.names = paste("X", 1:p, sep = "")
  if (dim(Y)[2] == 1) 
    Y.names = "Y"
  else {
    Y.names = dimnames(Y)[[2]]
    if (is.null(Y.names)) 
      Y.names = paste("Y", 1:q, sep = "")
  }
  ind.names = dimnames(X)[[1]]
  if (is.null(ind.names)) {
    ind.names = dimnames(Y)[[1]]
    rownames(X) = ind.names
  }
  if (is.null(ind.names)) {
    ind.names = 1:n
    rownames(X) = rownames(Y) = ind.names
  }
  if(scale.X){X = scale(X, center = TRUE, scale = TRUE)}
  if(scale.Y){Y = scale(Y, center = TRUE, scale = TRUE)}
  X.temp = X
  Y.temp = Y
  mat.t = matrix(nrow = n, ncol = ncomp)
  mat.u = matrix(nrow = n, ncol = ncomp)
  mat.a = matrix(nrow = p, ncol = ncomp)
  mat.b = matrix(nrow = q, ncol = ncomp)
  mat.c = matrix(nrow = p, ncol = ncomp)
  mat.d = matrix(nrow = q, ncol = ncomp)
  mat.e = matrix(nrow = q, ncol = ncomp)
  n.ones = rep(1, n)
  p.ones = rep(1, p)
  q.ones = rep(1, q)
  na.X = FALSE
  na.Y = FALSE
  is.na.X = is.na(X)
  is.na.Y = is.na(Y)
  if (any(is.na.X)) 
    na.X = TRUE
  if (any(is.na.Y)) 
    na.Y = TRUE
  for (h in 1:ncomp) {
    u = Y.temp[, 1]
    if (any(is.na(u))) 
      u[is.na(u)] = 0
    a.old = 0
    b.old = 0
    iter = 1
    if (na.X) {
      X.aux = X.temp
      X.aux[is.na.X] = 0
    }
    if (na.Y) {
      Y.aux = Y.temp
      Y.aux[is.na.Y] = 0
    }
    repeat {
      if (na.X) {
        a = crossprod(X.aux, u)
        U = drop(u) %o% p.ones
        U[is.na.X] = 0
        u.norm = crossprod(U)
        a = a/diag(u.norm)
        a = a/drop(sqrt(crossprod(a)))
        t = X.aux %*% a
        A = drop(a) %o% n.ones
        A[t(is.na.X)] = 0
        a.norm = crossprod(A)
        t = t/diag(a.norm)
      }
      else {
        a = crossprod(X.temp, u)/drop(crossprod(u))
        a = a/drop(sqrt(crossprod(a)))
        t = X.temp %*% a/drop(crossprod(a))
      }
      if (na.Y) {
        b = crossprod(Y.aux, t)
        T = drop(t) %o% q.ones
        T[is.na.Y] = 0
        t.norm = crossprod(T)
        b = b/diag(t.norm)
        u = Y.aux %*% b
        B = drop(b) %o% n.ones
        B[t(is.na.Y)] = 0
        b.norm = crossprod(B)
        u = u/diag(b.norm)
      }
      else {
        b = crossprod(Y.temp, t)/drop(crossprod(t))
        u = Y.temp %*% b/drop(crossprod(b))
      }
      if (crossprod(a - a.old) < tol) 
        break
      if (iter == max.iter) {
        warning(paste("Maximum number of iterations reached for dimension", 
                      h), call. = FALSE)
        break
      }
      a.old = a
      b.old = b
      iter = iter + 1
    }
    if (na.X) {
      X.aux = X.temp
      X.aux[is.na.X] = 0
      c = crossprod(X.aux, t)
      T = drop(t) %o% p.ones
      T[is.na.X] = 0
      t.norm = crossprod(T)
      c = c/diag(t.norm)
    }
    else {
      c = crossprod(X.temp, t)/drop(crossprod(t))
    }
    X.temp = X.temp - t %*% t(c)
    if (mode == "canonical") {
      if (na.Y) {
        Y.aux = Y.temp
        Y.aux[is.na.Y] = 0
        e = crossprod(Y.aux, u)
        U = drop(u) %o% q.ones
        U[is.na.Y] = 0
        u.norm = crossprod(U)
        e = e/diag(u.norm)
      }
      else {
        e = crossprod(Y.temp, u)/drop(crossprod(u))
      }
      Y.temp = Y.temp - u %*% t(e)
    }
    if (mode == "classic") 
      Y.temp = Y.temp - t %*% t(b)
    if (mode == "regression") {
      if (na.Y) {
        Y.aux = Y.temp
        Y.aux[is.na.Y] = 0
        d = crossprod(Y.aux, t)
        T = drop(t) %o% q.ones
        T[is.na.Y] = 0
        t.norm = crossprod(T)
        d = d/diag(t.norm)
      }
      else {
        d = crossprod(Y.temp, t)/drop(crossprod(t))
      }
      Y.temp = Y.temp - t %*% t(d)
    }
    if (mode == "invariant") 
      Y.temp = Y
    mat.t[, h] = t
    mat.u[, h] = u
    mat.a[, h] = a
    mat.b[, h] = b
    mat.c[, h] = c
    if (mode == "regression") 
      mat.d[, h] = d
    if (mode == "canonical") 
      mat.e[, h] = e
  }
  rownames(mat.a) = rownames(mat.c) = X.names
  rownames(mat.b) = Y.names
  rownames(mat.t) = rownames(mat.u) = ind.names
  comp = paste("comp", 1:ncomp)
  colnames(mat.t) = colnames(mat.u) = comp
  colnames(mat.a) = colnames(mat.b) = colnames(mat.c) = comp
  cl = match.call()
  cl[[1]] = as.name("pls")
  result = list(call = cl, X = X, Y = Y, ncomp = ncomp, mode = mode, 
                mat.c = mat.c, mat.d = mat.d, mat.e = mat.e, variates = list(X = mat.t, 
                                                                             Y = mat.u), loadings = list(X = mat.a, Y = mat.b), 
                names = list(X = X.names, Y = Y.names, indiv = ind.names))
  if (length(nzv$Position > 0)) 
    result$nzv = nzv
  class(result) = "pls"
  return(invisible(result))
}


predict.pls.cox=function(object, newdata, scale.X=TRUE, scale.Y=TRUE,...) 
{
  if (missing(newdata)) 
    stop("No new data available.")
  X = object$X
  Y = object$Y
  q = ncol(Y)
  p = ncol(X)
  if (length(dim(newdata)) == 2) {
    if (ncol(newdata) != p) 
      stop("'newdata' must be a numeric matrix with ncol = ", 
           p, " or a vector of length = ", p, ".")
  }
  if (length(dim(newdata)) == 0) {
    if (length(newdata) != p) 
      stop("'newdata' must be a numeric matrix with ncol = ", 
           p, " or a vector of length = ", p, ".")
    dim(newdata) = c(1, p)
  }
  ncomp = object$ncomp
  a = object$loadings$X
  b = object$loadings$Y
  c = object$mat.c
  if(scale.X){means.X = attr(X, "scaled:center")}
  if(scale.Y){means.Y = attr(Y, "scaled:center")}
  if(scale.X){sigma.X = attr(X, "scaled:scale")}
  if(scale.Y){sigma.Y = attr(Y, "scaled:scale")}
  newdata = as.matrix(newdata)
  ones = matrix(rep(1, nrow(newdata)), ncol = 1)
  B.hat = array(0, dim = c(p, q, ncomp))
  Y.hat = array(0, dim = c(nrow(newdata), q, ncomp))
  t.pred = array(0, dim = c(nrow(newdata), ncomp))
  for (h in 1:ncomp) {
    W = a[, 1:h] %*% solve(t(c[, 1:h]) %*% a[, 1:h])
    B = W %*% drop(t(b[, 1:h]))
    if(scale.Y){B = scale(B, center = FALSE, scale = 1/sigma.Y)}
    if(scale.X){B = as.matrix(scale(t(B), center = FALSE, scale = sigma.X))}
    if(!scale.X){B = as.matrix(t(B))}
    if(scale.X&scale.Y){intercept = -scale(B, center = FALSE, scale = 1/means.X)
                        intercept = matrix(apply(intercept, 1, sum) + means.Y, 
                                           nrow = 1)
                        Y.hat[, , h] = newdata %*% t(B) + ones %*% intercept}
    if(scale.X&!scale.Y){intercept = -scale(B, center = FALSE, scale = 1/means.X)
                         intercept = matrix(apply(intercept, 1, sum), nrow = 1)
                         Y.hat[, , h] = newdata %*% t(B) + ones %*% intercept}
    if(!scale.X&scale.Y){intercept = -B
                         intercept = matrix(apply(intercept, 1, sum) + means.Y, 
                                            nrow = 1)
                         Y.hat[, , h] = newdata %*% t(B) + ones %*% intercept}
    if(!scale.X&!scale.Y){Y.hat[, , h] = newdata %*% t(B)}
    if(!scale.X){t.pred[, h] = newdata %*% W[, h]}
    if(scale.X){t.pred[, h] = scale(newdata, center = means.X, scale = sigma.X) %*% W[, h]}
    B.hat[, , h] = B
  }
  rownames(t.pred) = rownames(newdata)
  colnames(t.pred) = paste("dim", c(1:ncomp), sep = " ")
  rownames(Y.hat) = rownames(newdata)
  colnames(Y.hat) = colnames(Y)
  return(invisible(list(predict = Y.hat, variates = t.pred, 
                        B.hat = B.hat)))
}

