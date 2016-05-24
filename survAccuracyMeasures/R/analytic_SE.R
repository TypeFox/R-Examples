##### SEMI-PARAMETRIC

WGT.FUN <- function(newdata, data, w.ptb=NULL, t0)
{
  ## ====================================##
  ## KM Estimator of Censoring Survival  ##
  ## ====================================##
  Ghat.FUN <- function(tt, Ti, Di,type='fl',w.ptb=NULL)
  {
    tmpind <- rank(tt); if(is.null(w.ptb)){w.ptb=rep(1,length(Ti))}
    summary(survfit(Surv(Ti,1-Di)~1, se.fit=F, type=type, weight=w.ptb), sort(tt))$surv[tmpind]
  }
  Ti = data[,1]; Di = data[,2]; tj = newdata[,1]; Wj = dj = newdata[,2]
  Wj[tj<=t0] = dj[tj<=t0]/Ghat.FUN(tj[tj<=t0],Ti,Di)
  Wj[tj >t0] = 1/Ghat.FUN(t0,Ti,Di)
  Wj
}


Est.Wexp.cpp <-
  function(data,N,RT.out,predict.time,uu0Vec,typexVec,typeyVec, resid.sco, fit.var, cutoffs) {
    
    if(missing(data))      { stop("Est.Wexp0: data not specified") }  
    
    #if( !("status" %in% names(data)) )  { stop(sprintf(errFormat,"status")) }
    #if( !("times" %in% names(data)) )  { stop(sprintf(errFormat,"times")) }
    
    numCuts = nrow(data)
    nr = numCuts
    if(!"wi" %in% names(data)) {
      data$weights=1
    }
    
    # First, fit the survival model    
    data = data[order(data$linearY),]   ## it is sorted it before this function; 
    
    Y  <- as.matrix(data[,!is.element(names(data), c("times", "zi", "status", "wi", "vi","Sy","linearY"))])
    
    np = dim(Y)[2]
    # fit  = coxph(Surv(data$times,data$status)~Y, 
    #              method="breslow", weight=data$weights)   
    
    # Doing riskmat, haz0 and time by hand since coxph.detail appears 
    #  to be a newer R feature & some users may have not updated their R.
    #    Note: this hazard is frequently normalized,
    #    by multiplying by exp(mean(data$Y)*fit$coef), but that is 
    #    not necessary here, as our haz0 below doesn't want it.
    
    
    #dataD    =  subset(data[order(data$times),],status==1)  

    if(is.na(cutoffs)[1]){
    Wexp.all <- getWEXP(as.matrix(data), as.matrix(Y), N, as.matrix(RT.out), predict.time, c(resid.sco), fit.var);
    }else{
      
      cutpos = sum.I(cutoffs,">=", data$linearY)
      
      Y.sub <- as.matrix(Y[cutpos,])
      subdata <- data[cutpos,]
     
      ncut = nrow(subdata)
      
      np = dim(Y)[2]
      

    Wexp.all <- getWEXPcutoff(as.matrix(data),
                               as.matrix(subdata),
                               Y = as.matrix(Y),
                               Y.sub,
                               N, as.matrix(RT.out), 
                               predict.time, c(resid.sco), fit.var, 
                               cutoffs);

      
    }
    ## now get iid expansion for other accuracy summaries
    ## global summaries 
    ## AUC = sum(RT.out$TPR*(RT.out$FPR-c(RT.out$FPR[-1],0)))
    mmm = length(RT.out$TPR)
    #ITPR = sum(RT.out$TPR*(RT.out$RiskT-c(0,RT.out$RiskT[-mmm])))
    #IFPR = sum(RT.out$FPR*(RT.out$RiskT-c(0,RT.out$RiskT[-mmm])))
    #IDI = ITPR - IFPR
    
    #Wexp.ITPR = Wexp.all[[4]]%*%(RT.out$RiskT-c(0,RT.out$RiskT[-mmm]))+
    #             (Wexp.all[[1]]-cbind(0,Wexp.all[[1]][,-mmm]))%*%RT.out$TPR
    #Wexp.IFPR = Wexp.all[[3]]%*%(RT.out$RiskT-c(0,RT.out$RiskT[-mmm]))+
    #             (Wexp.all[[1]]-cbind(0,Wexp.all[[1]][,-mmm]))%*%RT.out$FPR
    #Wexp.IDI= Wexp.ITPR - Wexp.IFPR 	
    #Wexp.AUC = Wexp.all[[4]]%*%(RT.out$FPR-c(RT.out$FPR[-1],0))+(Wexp.all[[3]]-cbind(Wexp.all[[3]][,-1],0))%*%RT.out$TPR
    Wexp.AUC = cbind(0,Wexp.all[[4]])%*%(c(1,RT.out$FPR)-c(RT.out$FPR,0))+
      (cbind(0,Wexp.all[[3]])-cbind(Wexp.all[[3]],0))%*%c(1,RT.out$TPR)
    
    
    
    
    if(!is.null(uu0Vec)){
      nvp = length(uu0Vec)
      Wexp.vp  = matrix(0,nr,nvp)
      
      for (pp in 1:nvp) {	
        uu0 = uu0Vec[pp]   
        
        uuk = sort(RT.out[,1]); 
        tmpind = sum.I(uu0,">=",uuk)
        ind0.y = match(typeyVec[pp],c("RiskT","v","FPR","TPR","rho","NPV","PPV"))
        
        Wexp.vp[,pp] = Wexp.all[[ind0.y]][,tmpind]
      }
      
    }else{
      Wexp.vp = NULL
    }
    
    
    list(Wexp.beta = Wexp.all[[8]], Wexp.AUC = Wexp.AUC,Wexp.vp=Wexp.vp)   
  }


##non parametric



