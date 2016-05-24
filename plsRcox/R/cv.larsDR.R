cv.larsDR =
  function (data, method = c("efron", "breslow"), nfold = 5, fraction = seq(0, 1, length = 100), plot.it = TRUE, 
            se = TRUE, givefold, scaleX = TRUE, scaleY = FALSE, folddetails = FALSE, allCVcrit=FALSE, 
            details=FALSE, namedataset="data", save=FALSE, verbose=TRUE,...)
  {
    try(attachNamespace("survival"),silent=TRUE)
    #on.exit(try(unloadNamespace("survival"),silent=TRUE))
    try(attachNamespace("rms"),silent=TRUE)
    on.exit(try(unloadNamespace("rms"),silent=TRUE))
    suppressMessages(try(attachNamespace("lars"),silent=TRUE))
    on.exit(try(unloadNamespace("lars"),silent=TRUE))

    cv.error1<-NULL;cv.error2<-NULL;cv.error3<-NULL;cv.error4<-NULL;cv.error5<-NULL;cv.error6<-NULL;cv.error7<-NULL;cv.error8<-NULL;cv.error9<-NULL;cv.error10<-NULL;cv.error11<-NULL;cv.error12<-NULL;cv.error13<-NULL;cv.error14<-NULL
    cv.se1<-NULL;cv.se2<-NULL;cv.se3<-NULL;cv.se4<-NULL;cv.se5<-NULL;cv.se6<-NULL;cv.se7<-NULL;cv.se8<-NULL;cv.se9<-NULL;cv.se10<-NULL;cv.se11<-NULL;cv.se12<-NULL;cv.se13<-NULL;cv.se14<-NULL
    lamin1<-NULL;lamin2<-NULL;lamin3<-NULL;lamin4<-NULL;lamin5<-NULL;lamin6<-NULL;lamin7<-NULL;lamin8<-NULL;lamin9<-NULL;lamin10<-NULL;lamin11<-NULL;lamin12<-NULL;lamin13<-NULL;lamin14<-NULL
    completed.cv1<-NULL;completed.cv2<-NULL;completed.cv3<-NULL;completed.cv4<-NULL;completed.cv5<-NULL;completed.cv6<-NULL;completed.cv7<-NULL;completed.cv8<-NULL;completed.cv9<-NULL;completed.cv10<-NULL;completed.cv11<-NULL;completed.cv12<-NULL;completed.cv13<-NULL;completed.cv14<-NULL
    
    method <- match.arg(method)
    x <- data$x
    time <- data$time
    status <- data$status
    n <- length(time)
    
    if(missing(givefold)){
      folds <- split(sample(seq(n)), rep(1:nfold, length = n))} else
      {
        folds <- givefold
      }
    
    number_ind = 14
    titlesCV = c("Cross-validated log-partial-likelihood","van Houwelingen Cross-validated log-partial-likelihood","iAUC_CD","iAUC_hc","iAUC_sh","iAUC_Uno","iAUC_hz.train","iAUC_hz.test","iAUC_survivalROC.train","iAUC_survivalROC.test","iBrierScore unw","iSchmidScore (robust BS) unw","iBrierScore w","iSchmidScore (robust BS) w")
    ylabsCV = c(rep("Minus log-partial-likelihood",2),rep("iAUC",8),rep("Prediction Error",4))
    xlabsCV = c(rep("L1 norm fraction",14))
    signCVerror = c(rep(1,2),rep(-1,8),rep(1,4))
    show_nbr_var = TRUE
    for(ind in 1:number_ind) {
      assign(paste("errormat",ind,sep=""),matrix(NA, length(fraction), nfold))
    }
    
    for (i in seq(nfold)) {
      for(ind in 1:number_ind) {
        assign(paste("pred",ind,sep=""),rep(NA, length(fraction)))
      }
      
      omit <- folds[[i]]
      trdata <- list(x = x[-omit, ], time = time[-omit], status = status[-omit])
      tsdata <- list(x = x[omit, ], time = time[omit], status = status[omit])
      
      YCsurv<-Surv(time[-omit],status[-omit])
      coxDR=coxph(YCsurv~1)
      DR_coxph=residuals(coxDR,type="deviance")
      
      if(!file.exists(paste("cv.larsDR_",namedataset,"_folds_",i,".RData",sep=""))){
        assign(paste("cv.larsDR_",namedataset,"_folds_",i,sep=""),lars(x=trdata$x,y=DR_coxph,use.Gram=FALSE,type="lasso",...))
        if(save){save(list=c(paste("cv.larsDR_",namedataset,"_folds_",i,sep="")),file=paste("cv.larsDR_",namedataset,"_folds_",i,".RData",sep=""))}
      } else {
        load(paste("cv.larsDR_",namedataset,"_folds_",i,".RData",sep=""))
      }
    
    coefmat <- coef(get(paste("cv.larsDR_",namedataset,"_folds_",i,sep="")),s=fraction,mode="fraction")
    predict.trainmat <- predict(get(paste("cv.larsDR_",namedataset,"_folds_",i,sep="")),s=fraction,mode="fraction",newx=trdata$x)$fit
    predictmat <- predict(get(paste("cv.larsDR_",namedataset,"_folds_",i,sep="")),s=fraction,mode="fraction",newx=tsdata$x)$fit
    nzb <- apply((!(coefmat==0)),1,sum)
    for(jj in 1:length(fraction)){
      coeffit <- t(coefmat[jj,,drop=FALSE])
      
      pred1[jj] <- logplik(x=tsdata$x, time=tsdata$time, status=tsdata$status, b=coeffit, method = method,return.all = FALSE)   #"efron"  match with loglik of coxph
      
      plfull <- logplik(x=x, time=time, status=status, b=coeffit, method = method,return.all = FALSE)   #"efron" 
      plminusk <- logplik(x=trdata$x, time=trdata$time, status=trdata$status, b=coeffit, method = method,return.all = FALSE)   #"efron" 
      pred2[jj] = plfull - plminusk
      
      
      Xlp <- rep(NA,length(time))
      Xlp[-omit] <- predict.trainmat[,jj]
      Xlp[omit] <- predictmat[,jj]
      
      assign(paste("dataset_",namedataset,"_",jj,sep=""),as.data.frame(cbind(time=time,status=status,Xlp=Xlp)))
      
      TR <- get(paste("dataset_",namedataset,"_",jj,sep=""))[-omit,]
      TE <- get(paste("dataset_",namedataset,"_",jj,sep=""))[omit,]
      survival.time <- get(paste("dataset_",namedataset,"_",jj,sep=""))[,"time"]
      survival.status <- get(paste("dataset_",namedataset,"_",jj,sep=""))[,"status"]
      tr.survival.time <- trdata$time
      tr.survival.status <- trdata$status
      te.survival.time <- tsdata$time
      te.survival.status <- tsdata$status
      
      #require(survival)
      Surv.rsp <- Surv(tr.survival.time, tr.survival.status)
      Surv.rsp.new <- Surv(te.survival.time, te.survival.status)
      
      train.fit <- coxph(Surv(time,status) ~ Xlp, x=TRUE, y=TRUE, method=method, data=TR, iter.max=0, init=1) #offset
      #library(rms)
      train.fit.cph <- cph(Surv(time,status) ~ Xlp, x=TRUE, y=TRUE, method=method, data=TR, iter.max=0, init=1) #offset
      lp <- predict(train.fit)
      lpnew <- predict(train.fit, newdata=TE)
      
      if(allCVcrit){
      AUCs <- getIndicCV(lp,lpnew,Surv.rsp,Surv.rsp.new,times.auc=seq(0,max(time),length.out=1000),times.prederr=seq(0,max(time),length.out=1000)[-(990:1000)],train.fit,plot.it=FALSE)
      pred3[jj] = AUCs$AUC_CD$iauc 
      pred4[jj] = AUCs$AUC_hc$iauc 
      pred5[jj] = AUCs$AUC_sh$iauc 
      pred6[jj] = AUCs$AUC_Uno$iauc
      pred7[jj] = AUCs$AUC_hz.train$iauc
      pred8[jj] = AUCs$AUC_hz.test$iauc
      pred9[jj] = AUCs$AUC_survivalROC.train$iauc
      pred10[jj] = AUCs$AUC_survivalROC.test$iauc
      pred11[jj] = AUCs$prederr$brier.unw$ierror
      pred12[jj] = AUCs$prederr$robust.unw$ierror
      pred13[jj] = AUCs$prederr$brier.w$ierror
      pred14[jj] = AUCs$prederr$robust.w$ierror
    }
    }
    if(allCVcrit){
    if(is.na(pred10[1])){pred10[1]<-.5}
    }
    if (length(omit) == 1){
      for(ind in 1:number_ind) {
        assign(paste("pred",ind,sep=""),matrix(get(paste("pred",ind,sep="")), nrow = 1))
      }
    }
    #if(any(is.na(pred10))){save(list=c("pred10"),file=paste(Predspath,"/failed.fold.cv.",typemodel,"_",namedataset,"_folds_",i,".RData",sep=""))}
    
    if(allCVcrit){
    errormat1[, i] <- ifelse(is.finite(pred1),-pred1/length(omit),NA)
    errormat2[, i] <- ifelse(is.finite(pred2),-pred2/length(omit),NA)
    errormat3[, i] <- ifelse(is.finite(pred3),pred3,NA)
    errormat4[, i] <- ifelse(is.finite(pred4),pred4,NA)
    errormat5[, i] <- ifelse(is.finite(pred5),pred5,NA)
    errormat6[, i] <- ifelse(is.finite(pred6),pred6,NA)
    errormat7[, i] <- ifelse(is.finite(pred7),pred7,NA)
    errormat8[, i] <- ifelse(is.finite(pred8),pred8,NA)
    errormat9[, i] <- ifelse(is.finite(pred9),pred9,NA)
    errormat10[, i] <- ifelse(is.finite(pred10),pred10,NA)
    errormat11[, i] <- ifelse(is.finite(pred11),pred11,NA)
    errormat12[, i] <- ifelse(is.finite(pred12),pred12,NA)
    errormat13[, i] <- ifelse(is.finite(pred13),pred13,NA)
    errormat14[, i] <- ifelse(is.finite(pred14),pred14,NA)
      } else {
        errormat2[, i] <- ifelse(is.finite(pred2),pred2,NA)        
      }
    if(verbose){cat("CV Fold", i, "\n")}
    rm(list=c(paste("cv.larsDR_",namedataset,"_folds_",i,sep="")))
  }
  if(allCVcrit){
  for(ind in 1:number_ind) {
    assign(paste("cv.error",ind,sep=""),apply(get(paste("errormat",ind,sep="")), 1, mean, na.rm=TRUE))
    assign(paste("completed.cv",ind,sep=""),is.finite(get(paste("errormat",ind,sep=""))))
    assign(paste("cv.se",ind,sep=""),sqrt(apply(get(paste("errormat",ind,sep="")), 1, var, na.rm=TRUE))/nfold)
    assign(paste("lamin",ind,sep=""),getmin2(fraction,signCVerror[ind]*get(paste("cv.error",ind,sep="")),get(paste("cv.se",ind,sep=""))))
    }} else {
      ind=2
      assign(paste("cv.error",ind,sep=""),apply(get(paste("errormat",ind,sep="")), 1, mean, na.rm=TRUE))
      assign(paste("completed.cv",ind,sep=""),is.finite(get(paste("errormat",ind,sep=""))))
      assign(paste("cv.se",ind,sep=""),sqrt(apply(get(paste("errormat",ind,sep="")), 1, var, na.rm=TRUE))/nfold)
      assign(paste("lamin",ind,sep=""),getmin2(fraction,signCVerror[ind]*get(paste("cv.error",ind,sep="")),get(paste("cv.se",ind,sep=""))))
    }
  YCsurv<-Surv(time,status)
  coxDR=coxph(YCsurv~1)
  DR_coxph=residuals(coxDR,type="deviance")
  
  
  assign(paste("system.time_larsDR_",namedataset,sep=""),system.time(assign("larsmodfull"   
  ,lars(x=x,y=DR_coxph,use.Gram=FALSE,type="lasso",...)),gcFirst=TRUE))
  
  coefmatfull <- coef(larsmodfull,s=fraction,mode="fraction")
  nzb <- apply((!(coefmatfull==0)),1,sum)

    sign.lambda=1

    if(allCVcrit){
    object <- list(fraction = fraction, cv.error1 = cv.error1, cv.error2 = cv.error2, cv.error3 = cv.error3, cv.error4 = cv.error4, cv.error5 = cv.error5, cv.error6 = cv.error6, cv.error7 = cv.error7, cv.error8 = cv.error8, cv.error9 = cv.error9, cv.error10 = cv.error10, cv.error11 = cv.error11, cv.error12 = cv.error12, cv.error13 = cv.error13, cv.error14 = cv.error14,
                   cv.se1 = cv.se1, cv.se2 = cv.se2, cv.se3 = cv.se3, cv.se4 = cv.se4, cv.se5 = cv.se5, cv.se6 = cv.se6, cv.se7 = cv.se7, cv.se8 = cv.se8, cv.se9 = cv.se9, cv.se10 = cv.se10, cv.se11 = cv.se11, cv.se12 = cv.se12, cv.se13 = cv.se13, cv.se14 = cv.se14, 
                   folds = folds, 
                   lambda.min1 = lamin1[[1]], lambda.1se1 = lamin1[[2]], lambda.min2 = lamin2[[1]], lambda.1se2 = lamin2[[2]], lambda.min3 = lamin3[[1]], lambda.1se3 = lamin3[[2]], lambda.min4 = lamin4[[1]], lambda.1se4 = lamin4[[2]], 
                   lambda.min5 = lamin5[[1]], lambda.1se5 = lamin5[[2]], lambda.min6 = lamin6[[1]], lambda.1se6 = lamin6[[2]], lambda.min7 = lamin7[[1]], lambda.1se7 = lamin7[[2]], lambda.min8 = lamin8[[1]], lambda.1se8 = lamin8[[2]], 
                   lambda.min9 = lamin9[[1]], lambda.1se9 = lamin9[[2]], lambda.min10 = lamin10[[1]], lambda.1se10 = lamin10[[2]], lambda.min11 = lamin11[[1]], lambda.1se11 = lamin11[[2]], lambda.min12 = lamin12[[1]], lambda.1se12 = lamin12[[2]], 
                   lambda.min13 = lamin13[[1]], lambda.1se13 = lamin13[[2]], lambda.min14 = lamin14[[1]], lambda.1se14 = lamin14[[2]], nzb=nzb, larsmodfull=larsmodfull 
                   )#sign.lambda=sign.lambda
    if(folddetails){object <- c(object,list(errormat1 = errormat1, errormat2 = errormat2, errormat3 = errormat3, errormat4 = errormat4, errormat5 = errormat5, errormat6 = errormat6, errormat7 = errormat7, errormat8 = errormat8, errormat9 = errormat9, errormat10 = errormat10, errormat11 = errormat11, errormat12 = errormat12, errormat13 = errormat13, errormat14 = errormat14, 
                   completed.cv1 = completed.cv1, completed.cv2 = completed.cv2, completed.cv3 = completed.cv3, completed.cv4 = completed.cv4, completed.cv5 = completed.cv5, completed.cv6 = completed.cv6, completed.cv7 = completed.cv7,
                   completed.cv8 = completed.cv8, completed.cv9 = completed.cv9, completed.cv10 = completed.cv10, completed.cv11 = completed.cv11, completed.cv12 = completed.cv12, completed.cv13 = completed.cv13, completed.cv14 = completed.cv14))}
    if(details){object <- c(object,list(All_indics=AUCs))}
    } else {
      object <- list(fraction = fraction,cv.error2=cv.error2,cv.se2=cv.se2,folds=folds,lambda.min2=lamin2[[1]],lambda.1se2=lamin2[[2]],nzb=nzb,larsmodfull=larsmodfull)
      if(folddetails){object <- c(object,list(errormat2 = errormat2, completed.cv2 = completed.cv2))}
    }
    
    
    if (plot.it) {
      if(allCVcrit){
        for(ind in 1:number_ind) {
        if((ind%% 4)==1){dev.new();layout(matrix(1:4,nrow=2))}
        plot(sign.lambda*(fraction[!is.nan(get(paste("cv.error",ind,sep="")))]), get(paste("cv.error",ind,sep=""))[!is.nan(get(paste("cv.error",ind,sep="")))], type = "l", xlim =range(fraction), ylim = range(c(get(paste("cv.error",ind,sep=""))[!is.nan(get(paste("cv.error",ind,sep="")))] -
                                                                                                                                                                                                                  get(paste("cv.se",ind,sep=""))[!is.nan(get(paste("cv.error",ind,sep="")))], get(paste("cv.error",ind,sep=""))[!is.nan(get(paste("cv.error",ind,sep="")))] + get(paste("cv.se",ind,sep=""))[!is.nan(get(paste("cv.error",ind,sep="")))])), xlab = xlabsCV[ind], ylab = ylabsCV[ind],
             main = titlesCV[ind]
        )
        abline(v = sign.lambda*getElement(object,paste("lambda.min",ind,sep="")), lty = 3)
        abline(v = sign.lambda*getElement(object,paste("lambda.1se",ind,sep="")), lty = 3, col="red")
        if(show_nbr_var){axis(side = 3, at = sign.lambda*object$fraction, labels = paste(object$nzb), tick = FALSE, line = -1)}
        if (se)
          segments(sign.lambda*(fraction[!is.nan(get(paste("cv.error",ind,sep="")))]), get(paste("cv.error",ind,sep=""))[!is.nan(get(paste("cv.error",ind,sep="")))] - get(paste("cv.se",ind,sep=""))[!is.nan(get(paste("cv.error",ind,sep="")))], sign.lambda*(fraction[!is.nan(get(paste("cv.error",ind,sep="")))]), get(paste("cv.error",ind,sep=""))[!is.nan(get(paste("cv.error",ind,sep="")))] +
                     get(paste("cv.se",ind,sep=""))[!is.nan(get(paste("cv.error",ind,sep="")))])
      }        
      layout(1)
    } else {
      ind=2
      plot(sign.lambda*(fraction[!is.nan(get(paste("cv.error",ind,sep="")))]), get(paste("cv.error",ind,sep=""))[!is.nan(get(paste("cv.error",ind,sep="")))], type = "l", xlim =range(fraction), ylim = range(c(get(paste("cv.error",ind,sep=""))[!is.nan(get(paste("cv.error",ind,sep="")))] -
                                                                                                                                                                                                                  get(paste("cv.se",ind,sep=""))[!is.nan(get(paste("cv.error",ind,sep="")))], get(paste("cv.error",ind,sep=""))[!is.nan(get(paste("cv.error",ind,sep="")))] + get(paste("cv.se",ind,sep=""))[!is.nan(get(paste("cv.error",ind,sep="")))])), xlab = xlabsCV[ind], ylab = ylabsCV[ind],
           main = titlesCV[ind]
      )
      abline(v = sign.lambda*getElement(object,paste("lambda.min",ind,sep="")), lty = 3)
      abline(v = sign.lambda*getElement(object,paste("lambda.1se",ind,sep="")), lty = 3, col="red")
      if(show_nbr_var){axis(side = 3, at = sign.lambda*object$fraction, labels = paste(object$nzb), tick = FALSE, line = -1)}
      if (se)
        segments(sign.lambda*(fraction[!is.nan(get(paste("cv.error",ind,sep="")))]), get(paste("cv.error",ind,sep=""))[!is.nan(get(paste("cv.error",ind,sep="")))] - get(paste("cv.se",ind,sep=""))[!is.nan(get(paste("cv.error",ind,sep="")))], sign.lambda*(fraction[!is.nan(get(paste("cv.error",ind,sep="")))]), get(paste("cv.error",ind,sep=""))[!is.nan(get(paste("cv.error",ind,sep="")))] +
                   get(paste("cv.se",ind,sep=""))[!is.nan(get(paste("cv.error",ind,sep="")))])
    }
    }
    invisible(object)
  }
