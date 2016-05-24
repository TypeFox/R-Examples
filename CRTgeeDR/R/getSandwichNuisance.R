### Calculate the sandwich estimator accounting for estimation of nuisance parameters in PS and OM.
getSandwichNuisance = function(Y, X,X.t,X.c, B, beta,off, id, R.alpha.inv, phi, InvLinkDeriv, InvLink, VarFun, hessMat, StdErr, dInvLinkdEta, BlockDiag, sqrtW,W,included,typeweights,pi.a,nameTRT,propensity.score,om.t,om.c,data,nameY,nameMISS,print.log){
  
  eta <- as.vector(X%*%beta) + off
  diag(dInvLinkdEta) <- InvLinkDeriv(eta)
  mu <- InvLink(eta)			
  diag(StdErr) <- sqrt(1/VarFun(mu))
  
  nn<-length(Y)
  StdErr.c <- Diagonal(nn)
  dInvLinkdEta.c <- Diagonal(nn)
  eta.c <- as.vector(X.c%*%beta) + off    
  diag(dInvLinkdEta.c) <- InvLinkDeriv(eta.c)
  mu.c <- InvLink(eta.c)	    
  diag(StdErr.c) <- sqrt(1/VarFun(mu.c))
  
  nn<-length(Y)
  StdErr.t <- Diagonal(nn)
  dInvLinkdEta.t <- Diagonal(nn)
  eta.t <- as.vector(X.t%*%beta) + off
  diag(dInvLinkdEta.t) <- InvLinkDeriv(eta.t)
  mu.t <- InvLink(eta.t)      
  diag(StdErr.t) <- sqrt(1/VarFun(mu.t))
  
  if((!is.null(om.t))&(!is.null(om.c))){
    
    scoreDiag <- Diagonal(x= Y - B[,"Bi"])
    scoreDiag.t <-Diagonal(x= B[,"B.t"]-mu.t)
    scoreDiag.c <-Diagonal(x= B[,"B.c"]-mu.c)
  }
  scoreDiag.ps <- Diagonal(x= Y - mu)
  
  design.weights<-NULL
  if(!is.null(propensity.score)){   
    ######## Nuisance parameters for weights
    design.weights<-model.matrix(propensity.score)
    var.w <- getfam(propensity.score$family)$VarFun
    HW<--crossprod(design.weights,Diagonal(x=var.w(fitted(propensity.score)))%*%design.weights)
    SW<-t(design.weights)%*%Diagonal(x=((1-data[,nameMISS])-fitted(propensity.score))) #%*%Diagonal(x=((1-data[,nameMISS])-fitted(propensity.score)))%*%design.weights
  }
  
  design.om.trt<-NULL
  design.om.ctrl<-NULL
  design.om.trt.all<-NULL
  design.om.trt.all<-NULL
  if((!is.null(om.t))&(!is.null(om.c))){
    ######## Nuisance parameters for outcome
    design.om.trt<-model.matrix(om.t)
    design.om.ctrl<-model.matrix(om.c)  
    design.om.trt.all<-model.matrix(as.formula(paste("~",as.character(om.t$formula[3]))),data=data)
    design.om.ctrl.all<-model.matrix(as.formula(paste("~",as.character(om.c$formula[3]))),data=data)
    
    var.om <- getfam(om.t$family)$VarFun
    HB.trt<--crossprod(design.om.trt,Diagonal(x=var.om(fitted(om.t)))%*%design.om.trt)
    HB.ctrl<--crossprod(design.om.ctrl,Diagonal(x=var.om(fitted(om.c)))%*%design.om.ctrl)
    
    YsansNA<-ifelse(is.na(data[,nameY]),0,data[,nameY])
    predtrt<-as.data.frame(InvLink(as.matrix(design.om.trt.all)%*%om.t$coefficients))
    diag.trt<-unlist((YsansNA-predtrt)*as.numeric(I(!is.na(data[,nameY]))))
    predctrl<-as.data.frame(InvLink(as.matrix(design.om.ctrl.all)%*%om.c$coefficients))
    diag.ctrl<-unlist((YsansNA-predctrl)*as.numeric(I(!is.na(data[,nameY]))))
    SB.trt.all<-t(design.om.trt.all)%*%Diagonal(x=as.numeric(diag.trt))#%*%Diagonal(x=(data.trt[which(!is.na(data.trt[,nameY])),nameY]-fitted(om.t)))%*%design.om.trt
    SB.ctrl.all<-t(design.om.ctrl.all)%*%Diagonal(x=as.vector(diag.ctrl))#%*%Diagonal(x=(data.ctrl[which(!is.na(data.ctrl[,nameY])),nameY]-fitted(om.c)))%*%design.om.ctrl
  }
  
  if((!is.null(om.t))&(!is.null(om.c))&(!is.null(propensity.score))){
    
    nuisance<-c(beta,propensity.score$coefficients,om.t$coefficients,om.c$coefficients)
    nuisance<-replace(nuisance, is.na(nuisance), 0)
    
    jacobian.nuisance<-jacobian(funcoptvarnuisanceOMPS, nuisance,method="Richardson",sqrtW=sqrtW,
                                design.weights=design.weights,design.om.trt.all=design.om.trt.all,design.om.ctrl.all=design.om.ctrl.all,
                                ,design.om.trt=design.om.trt,design.om.ctrl=design.om.ctrl,
                                nameTRT=nameTRT,nameY=nameY,nameMISS=nameMISS,
                                data=data,X=X,X.c=X.c,X.t=X.t,InvLinkDeriv=InvLinkDeriv, InvLink=InvLink,
                                VarFun=VarFun,dInvLinkdEta=dInvLinkdEta,StdErr=StdErr,R.alpha.inv=R.alpha.inv,Y=Y,off=off,pi.a=pi.a,B=B,typeweights=typeweights)    
    
    if(typeweights=="GENMOD"){
      SEE<-crossprod(sqrtW %*% StdErr %*% dInvLinkdEta %*% X, R.alpha.inv %*% sqrtW %*% StdErr %*% scoreDiag)+
        (pi.a)*crossprod( StdErr.t %*% dInvLinkdEta.t %*% X.t, R.alpha.inv %*% StdErr.t %*% scoreDiag.t)+
        ((1-pi.a))*crossprod( StdErr.c %*% dInvLinkdEta.c %*% X.c, R.alpha.inv  %*% StdErr.c %*% scoreDiag.c)
    }else{
      SEE<-crossprod( StdErr %*% dInvLinkdEta %*% X, R.alpha.inv %*% StdErr %*% W %*% scoreDiag)+
        (pi.a)*crossprod( StdErr.t %*% dInvLinkdEta.t %*% X.t, R.alpha.inv %*% StdErr.t %*% scoreDiag.t)+
        ((1-pi.a))*crossprod( StdErr.c %*% dInvLinkdEta.c %*% X.c, R.alpha.inv  %*% StdErr.c %*% scoreDiag.c)      
    }
    
    
    stackEE<-rbind(as.matrix(SEE),as.matrix(SW),as.matrix(SB.trt.all),as.matrix(SB.ctrl.all))    
    stackEE<-replace(stackEE, is.na(stackEE), 0)
    
    stackEEprod<-stackEE%*%BlockDiag%*%t(stackEE)
    
    sandadjvar<-t(ginv(jacobian.nuisance))%*%stackEEprod%*%ginv(jacobian.nuisance)
    
  }
  
  if((!is.null(om.t))&(!is.null(om.c))&(is.null(propensity.score))){
    
    nuisance<-c(beta,om.t$coefficients,om.c$coefficients)
    nuisance<-replace(nuisance, is.na(nuisance), 0)
    jacobian.nuisance<-jacobian(funcoptvarnuisanceOM, nuisance,method="Richardson",sqrtW=sqrtW,
                                design.weights=design.weights,design.om.trt.all=design.om.trt.all,design.om.ctrl.all=design.om.ctrl.all,
                                ,design.om.trt=design.om.trt,design.om.ctrl=design.om.ctrl,
                                nameTRT=nameTRT,nameY=nameY,nameMISS=nameMISS,
                                data=data,X=X,X.c=X.c,X.t=X.t,InvLinkDeriv=InvLinkDeriv, InvLink=InvLink,
                                VarFun=VarFun,dInvLinkdEta=dInvLinkdEta,StdErr=StdErr,R.alpha.inv=R.alpha.inv,Y=Y,off=off,pi.a=pi.a,B=B,typeweights=typeweights)    
    if(!is.null(typeweights)){
      if(typeweights=="GENMOD"){ 
        SEE<-crossprod(sqrtW %*% StdErr %*% dInvLinkdEta %*% X, R.alpha.inv %*% sqrtW %*% StdErr %*% scoreDiag)+
          (pi.a)*crossprod( StdErr.t %*% dInvLinkdEta.t %*% X.t, R.alpha.inv %*% StdErr.t %*% scoreDiag.t)+
          ((1-pi.a))*crossprod( StdErr.c %*% dInvLinkdEta.c %*% X.c, R.alpha.inv  %*% StdErr.c %*% scoreDiag.c)
      }else{
        SEE<-crossprod(StdErr %*% dInvLinkdEta %*% X, R.alpha.inv %*% sqrtW %*% sqrtW %*% StdErr %*% scoreDiag)+
          (pi.a)*crossprod( StdErr.t %*% dInvLinkdEta.t %*% X.t, R.alpha.inv %*% StdErr.t %*% scoreDiag.t)+
          ((1-pi.a))*crossprod( StdErr.c %*% dInvLinkdEta.c %*% X.c, R.alpha.inv  %*% StdErr.c %*% scoreDiag.c)
      }
    }else{
      SEE<-crossprod(StdErr %*% dInvLinkdEta %*% X, R.alpha.inv %*% StdErr %*% scoreDiag)+
        (pi.a)*crossprod( StdErr.t %*% dInvLinkdEta.t %*% X.t, R.alpha.inv %*% StdErr.t %*% scoreDiag.t)+
        ((1-pi.a))*crossprod( StdErr.c %*% dInvLinkdEta.c %*% X.c, R.alpha.inv  %*% StdErr.c %*% scoreDiag.c)
      
    }
    
    stackEE<-rbind(as.matrix(SEE),as.matrix(SB.trt.all),as.matrix(SB.ctrl.all))    
    stackEE<-replace(stackEE, is.na(stackEE), 0)
    
    stackEEprod<-stackEE%*%BlockDiag%*%t(stackEE)
    
    sandadjvar<-t(ginv(jacobian.nuisance))%*%stackEEprod%*%ginv(jacobian.nuisance)
    
  }    
  
  if((is.null(om.t))&(is.null(om.c))&(!is.null(propensity.score))){
    
    nuisance<-c(beta,propensity.score$coefficients)
    nuisance<-replace(nuisance, is.na(nuisance), 0)
    
    jacobian.nuisance<-jacobian(funcoptvarnuisancePS, nuisance,method="Richardson",sqrtW=sqrtW,
                                design.weights=design.weights,design.om.trt.all=design.om.trt.all,design.om.ctrl.all=design.om.ctrl.all,
                                design.om.trt=design.om.trt,design.om.ctrl=design.om.ctrl,
                                nameTRT=nameTRT,nameY=nameY,nameMISS=nameMISS,
                                data=data,X=X,X.c=X.c,X.t=X.t,InvLinkDeriv=InvLinkDeriv, InvLink=InvLink,
                                VarFun=VarFun,dInvLinkdEta=dInvLinkdEta,StdErr=StdErr,R.alpha.inv=R.alpha.inv,Y=Y,off=off,pi.a=pi.a,B=B,typeweights=typeweights)    
    
    if(is.null(B)){
      if(typeweights=="GENMOD"){
        SEE<-crossprod(sqrtW %*% StdErr %*% dInvLinkdEta %*% X, R.alpha.inv %*% sqrtW %*% StdErr %*% scoreDiag.ps)
      }else{
        SEE<-crossprod( StdErr %*% dInvLinkdEta %*% X, R.alpha.inv %*% StdErr %*% W %*% scoreDiag.ps)
      }
    }else{
      if(typeweights=="GENMOD"){
        SEE<-crossprod(sqrtW %*% StdErr %*% dInvLinkdEta %*% X, R.alpha.inv %*% sqrtW %*% StdErr %*% scoreDiag)+
          (pi.a)*crossprod( StdErr.t %*% dInvLinkdEta.t %*% X.t, R.alpha.inv %*% StdErr.t %*% scoreDiag.t)+
          ((1-pi.a))*crossprod( StdErr.c %*% dInvLinkdEta.c %*% X.c, R.alpha.inv  %*% StdErr.c %*% scoreDiag.c)
      }else{
        SEE<-crossprod( StdErr %*% dInvLinkdEta %*% X, R.alpha.inv  %*% StdErr %*%W %*% scoreDiag)+
          (pi.a)*crossprod( StdErr.t %*% dInvLinkdEta.t %*% X.t, R.alpha.inv %*% StdErr.t %*% scoreDiag.t)+
          ((1-pi.a))*crossprod( StdErr.c %*% dInvLinkdEta.c %*% X.c, R.alpha.inv  %*% StdErr.c %*% scoreDiag.c)       
      }
    }
    
    
    stackEE<-rbind(as.matrix(SEE),as.matrix(SW))    
    stackEE<-replace(stackEE, is.na(stackEE), 0)
    
    stackEEprod<-stackEE%*%BlockDiag%*%t(stackEE)
    
    sandadjvar<-t(ginv(jacobian.nuisance))%*%stackEEprod%*%ginv(jacobian.nuisance)
  }
  
  if(print.log){
    print("Nuisance-adjusted sandwich")
    print(sandadjvar)
  }
  
  sandadjvar<-sandadjvar[1: dim(X)[2],1: dim(X)[2]]  
  
  
  return(list(jacobian.nuisance=jacobian.nuisance,sandadjvar=sandadjvar))
}


### Function computing the join score equation (EE,PS,OM) depending on main and nuisance parameters values
funcoptvarnuisanceOMPS <- function(nuisance,sqrtW,
                                   design.weights,design.om.trt.all,design.om.ctrl.all,design.om.trt,design.om.ctrl,
                                   nameTRT,nameY,nameMISS,
                                   data,X,X.c,X.t,InvLinkDeriv, InvLink, VarFun,dInvLinkdEta,StdErr,R.alpha.inv,Y,off,pi.a,B,typeweights){
  
  etaEE<-nuisance[1:ncol(X)]+1
  etaW<-nuisance[(ncol(X)+1):(ncol(X)+ncol(design.weights))]
  etaBtrt<-nuisance[(ncol(X)+ncol(design.weights)+1):(ncol(X)+ncol(design.weights)+ncol(design.om.trt.all))]
  etaBctrl<-nuisance[(ncol(X)+ncol(design.weights)+ncol(design.om.trt.all)+1):(ncol(X)+ncol(design.weights)+ncol(design.om.trt.all)+ncol(design.om.ctrl.all))]
  
  
  eta <- as.vector(X%*%etaEE) + off
  diag(dInvLinkdEta) <- InvLinkDeriv(eta)
  mu <- InvLink(eta)  		
  diag(StdErr) <- sqrt(1/VarFun(mu))
  
  nn<-length(Y)
  StdErr.c <- Diagonal(nn)
  dInvLinkdEta.c <- Diagonal(nn)
  eta.c <- as.vector(X.c%*%etaEE) + off    
  diag(dInvLinkdEta.c) <- InvLinkDeriv(eta.c)
  mu.c <- InvLink(eta.c)	    
  diag(StdErr.c) <- sqrt(1/VarFun(mu.c))
  
  nn<-length(Y)
  StdErr.t <- Diagonal(nn)
  dInvLinkdEta.t <- Diagonal(nn)
  eta.t <- as.vector(X.t%*%etaEE) + off
  diag(dInvLinkdEta.t) <- InvLinkDeriv(eta.t)
  mu.t <- InvLink(eta.t)      
  diag(StdErr.t) <- sqrt(1/VarFun(mu.t))
  
  temp<-sqrt(as.numeric(diag(sqrtW)>0)/(exp((design.weights%*%etaW))/(1+exp((design.weights%*%etaW)))) )
  temp<-ifelse(is.na(temp),0,temp)
  sqrtW.temp<-Diagonal(x= temp)
  temp<-(((1-data[,nameMISS])-(exp((design.weights%*%etaW))/(1+exp((design.weights%*%etaW))))))
  SW<-as.vector(t(design.weights)%*%ifelse(is.na(temp),0,temp))
  
  B.temp<-as.data.frame(InvLink(as.matrix(design.om.trt.all)%*%etaBtrt))
  names(B.temp)<-c("B.t")
  B.temp[,"B.c"]<-as.data.frame(InvLink(as.matrix(design.om.ctrl.all)%*%etaBctrl))
  temp<-as.data.frame(cbind(as.data.frame(X)[,colnames(as.data.frame(X))==nameTRT],B.temp[,"B.c"],B.temp[,"B.t"])) 
  names(temp)<-c(nameTRT,"B.c","B.t")
  Bi.temp<-ifelse(temp[,nameTRT]==1,temp[,"B.t"],temp[,"B.c"])
  B.temp<-as.data.frame(cbind(temp[,"B.c"],temp[,"B.t"],Bi.temp))
  names(B.temp)<-c("B.c","B.t","Bi")
  
  predtrt<-as.data.frame(InvLink(as.matrix(design.om.trt)%*%etaBtrt))
  predctrl<-as.data.frame(InvLink(as.matrix(design.om.ctrl)%*%etaBctrl))  
  SB.trt<-as.vector(t(design.om.trt)%*%as.matrix((data[which((!is.na(data[,nameY]))&(data[,nameTRT]==1)),nameY]-predtrt)))
  SB.ctrl<-as.vector(t(design.om.ctrl)%*%as.matrix((data[which((!is.na(data[,nameY]))&(data[,nameTRT]==0)),nameY]-predctrl)))
  
  if(typeweights=="GENMOD"){
    scoreEE<-as.vector(crossprod(sqrtW.temp %*% StdErr %*%dInvLinkdEta %*%X , R.alpha.inv %*% sqrtW.temp %*% StdErr %*% (Y - B.temp[,"Bi"])) +
                         (1-pi.a)*crossprod(  StdErr.c %*%dInvLinkdEta.c%*%X.c , R.alpha.inv %*%  StdErr.c %*% (B.temp[,"B.c"]-mu.c))+
                         (pi.a)*crossprod(   StdErr.t %*%dInvLinkdEta.t %*%X.t , R.alpha.inv %*%  StdErr.t %*% (B.temp[,"B.t"]-mu.t)))
  }else{
    scoreEE<-as.vector(crossprod( StdErr %*%dInvLinkdEta %*%X , R.alpha.inv %*% StdErr %*% sqrtW.temp %*% sqrtW.temp %*% (Y - B.temp[,"Bi"])) +
                         (1-pi.a)*crossprod(  StdErr.c %*%dInvLinkdEta.c%*%X.c , R.alpha.inv %*%  StdErr.c %*% (B.temp[,"B.c"]-mu.c))+
                         (pi.a)*crossprod(   StdErr.t %*%dInvLinkdEta.t %*%X.t , R.alpha.inv %*%  StdErr.t %*% (B.temp[,"B.t"]-mu.t)))    
  }
  funcoptvarnuisance<-c(scoreEE,SW,SB.trt,SB.ctrl)
} 


### Function computing the join score equation (EE,OM) depending on main and nuisance parameters values
funcoptvarnuisanceOM <- function(nuisance,sqrtW,
                                 design.weights,design.om.trt.all,design.om.ctrl.all,design.om.trt,design.om.ctrl,
                                 nameTRT,nameY,nameMISS,
                                 data,X,X.c,X.t,InvLinkDeriv, InvLink, VarFun,dInvLinkdEta,StdErr,R.alpha.inv,Y,off,pi.a,B,typeweights){
  
  etaEE<-nuisance[1:ncol(X)]+1
  etaBtrt<-nuisance[(ncol(X)+1):(ncol(X)+ncol(design.om.trt.all))]
  etaBctrl<-nuisance[(ncol(X)+ncol(design.om.trt.all)+1):(ncol(X)+ncol(design.om.trt.all)+ncol(design.om.ctrl.all))]
  
  
  eta <- as.vector(X%*%etaEE) + off
  diag(dInvLinkdEta) <- InvLinkDeriv(eta)
  mu <- InvLink(eta)  		
  diag(StdErr) <- sqrt(1/VarFun(mu))
  
  nn<-length(Y)
  StdErr.c <- Diagonal(nn)
  dInvLinkdEta.c <- Diagonal(nn)
  eta.c <- as.vector(X.c%*%etaEE) + off    
  diag(dInvLinkdEta.c) <- InvLinkDeriv(eta.c)
  mu.c <- InvLink(eta.c)	    
  diag(StdErr.c) <- sqrt(1/VarFun(mu.c))
  
  nn<-length(Y)
  StdErr.t <- Diagonal(nn)
  dInvLinkdEta.t <- Diagonal(nn)
  eta.t <- as.vector(X.t%*%etaEE) + off
  diag(dInvLinkdEta.t) <- InvLinkDeriv(eta.t)
  mu.t <- InvLink(eta.t)      
  diag(StdErr.t) <- sqrt(1/VarFun(mu.t))
  
  
  B.temp<-as.data.frame(InvLink(as.matrix(design.om.trt.all)%*%etaBtrt))
  names(B.temp)<-c("B.t")
  B.temp[,"B.c"]<-as.data.frame(InvLink(as.matrix(design.om.ctrl.all)%*%etaBctrl))
  temp<-as.data.frame(cbind(as.data.frame(X)[,colnames(as.data.frame(X))==nameTRT],B.temp[,"B.c"],B.temp[,"B.t"])) 
  names(temp)<-c(nameTRT,"B.c","B.t")
  Bi.temp<-ifelse(temp[,nameTRT]==1,temp[,"B.t"],temp[,"B.c"])
  B.temp<-as.data.frame(cbind(temp[,"B.c"],temp[,"B.t"],Bi.temp))
  names(B.temp)<-c("B.c","B.t","Bi")
  
  predtrt<-as.data.frame(InvLink(as.matrix(design.om.trt)%*%etaBtrt))
  predctrl<-as.data.frame(InvLink(as.matrix(design.om.ctrl)%*%etaBctrl))  
  SB.trt<-as.vector(t(design.om.trt)%*%as.matrix((data[which((!is.na(data[,nameY]))&(data[,nameTRT]==1)),nameY]-predtrt)))
  SB.ctrl<-as.vector(t(design.om.ctrl)%*%as.matrix((data[which((!is.na(data[,nameY]))&(data[,nameTRT]==0)),nameY]-predctrl)))
  
  if(!is.null(typeweights)){
    if(typeweights=="GENMOD"){
      scoreEE<-as.vector(crossprod(sqrtW %*% StdErr %*%dInvLinkdEta %*%X , R.alpha.inv %*% sqrtW %*% StdErr %*% (Y - B.temp[,"Bi"])) +
                           (1-pi.a)*crossprod(  StdErr.c %*%dInvLinkdEta.c%*%X.c , R.alpha.inv %*%  StdErr.c %*% (B.temp[,"B.c"]-mu.c))+
                           (pi.a)*crossprod(   StdErr.t %*%dInvLinkdEta.t %*%X.t , R.alpha.inv %*%  StdErr.t %*% (B.temp[,"B.t"]-mu.t)))
    }else{
      scoreEE<-as.vector(crossprod(StdErr %*%dInvLinkdEta %*%X , R.alpha.inv %*% sqrtW %*% sqrtW %*% StdErr %*% (Y - B.temp[,"Bi"])) +
                           (1-pi.a)*crossprod(  StdErr.c %*%dInvLinkdEta.c%*%X.c , R.alpha.inv %*%  StdErr.c %*% (B.temp[,"B.c"]-mu.c))+
                           (pi.a)*crossprod(   StdErr.t %*%dInvLinkdEta.t %*%X.t , R.alpha.inv %*%  StdErr.t %*% (B.temp[,"B.t"]-mu.t)))  
    }
  }else{
    scoreEE<-as.vector(crossprod(StdErr %*%dInvLinkdEta %*%X , R.alpha.inv %*% StdErr %*% (Y - B.temp[,"Bi"])) +
                         (1-pi.a)*crossprod(  StdErr.c %*%dInvLinkdEta.c%*%X.c , R.alpha.inv %*%  StdErr.c %*% (B.temp[,"B.c"]-mu.c))+
                         (pi.a)*crossprod(   StdErr.t %*%dInvLinkdEta.t %*%X.t , R.alpha.inv %*%  StdErr.t %*% (B.temp[,"B.t"]-mu.t)))     
  }
  funcoptvarnuisance<-c(scoreEE,SB.trt,SB.ctrl)
} 

### Function computing the join score equation (EE,PS) depending on main and nuisance parameters values
funcoptvarnuisancePS <- function(nuisance,sqrtW,
                                 design.weights,design.om.trt.all,design.om.ctrl.all,design.om.trt,design.om.ctrl,
                                 nameTRT,nameY,nameMISS,
                                 data,X,X.c,X.t,InvLinkDeriv, InvLink, VarFun,dInvLinkdEta,StdErr,R.alpha.inv,Y,off,pi.a,B,typeweights){
  
  etaEE<-nuisance[1:ncol(X)]+1
  etaW<-nuisance[(ncol(X)+1):(ncol(X)+ncol(design.weights))]
  
  eta <- as.vector(X%*%etaEE) + off
  diag(dInvLinkdEta) <- InvLinkDeriv(eta)
  mu <- InvLink(eta)  		
  diag(StdErr) <- sqrt(1/VarFun(mu))
  
  nn<-length(Y)
  StdErr.c <- Diagonal(nn)
  dInvLinkdEta.c <- Diagonal(nn)
  eta.c <- as.vector(X.c%*%etaEE) + off    
  diag(dInvLinkdEta.c) <- InvLinkDeriv(eta.c)
  mu.c <- InvLink(eta.c)	    
  diag(StdErr.c) <- sqrt(1/VarFun(mu.c))
  
  nn<-length(Y)
  StdErr.t <- Diagonal(nn)
  dInvLinkdEta.t <- Diagonal(nn)
  eta.t <- as.vector(X.t%*%etaEE) + off
  diag(dInvLinkdEta.t) <- InvLinkDeriv(eta.t)
  mu.t <- InvLink(eta.t)      
  diag(StdErr.t) <- sqrt(1/VarFun(mu.t))
  
  
  temp<-sqrt(as.numeric(diag(sqrtW)>0)/(exp((design.weights%*%etaW))/(1+exp((design.weights%*%etaW)))) )
  temp<-ifelse(is.na(temp),0,temp)
  sqrtW.temp<-Diagonal(x= temp)
  temp<-(((1-data[,nameMISS])-(exp((design.weights%*%etaW))/(1+exp((design.weights%*%etaW))))))
  SW<-as.vector(t(design.weights)%*%ifelse(is.na(temp),0,temp))
  
  if(is.null(B)){
    if(typeweights=="GENMOD"){
      scoreEE<-as.vector(crossprod(sqrtW %*% StdErr %*%dInvLinkdEta %*%X , R.alpha.inv %*% sqrtW %*% StdErr %*% (Y - mu)))
    }else{
      scoreEE<-as.vector(crossprod(StdErr %*%dInvLinkdEta %*%X , R.alpha.inv %*% sqrtW %*% sqrtW %*% StdErr %*% (Y - mu)))
    }
  }else{
    if(typeweights=="GENMOD"){
      scoreEE<-as.vector(crossprod(sqrtW.temp %*% StdErr %*%dInvLinkdEta %*%X , R.alpha.inv %*% sqrtW.temp %*% StdErr %*% (Y - B[,"Bi"])) +
                           (1-pi.a)*crossprod(  StdErr.c %*%dInvLinkdEta.c%*%X.c , R.alpha.inv %*%  StdErr.c %*% (B[,"B.c"]-mu.c))+
                           (pi.a)*crossprod(   StdErr.t %*%dInvLinkdEta.t %*%X.t , R.alpha.inv %*%  StdErr.t %*% (B[,"B.t"]-mu.t)))
    }else{
      scoreEE<-as.vector(crossprod(StdErr %*%dInvLinkdEta %*%X , R.alpha.inv %*% sqrtW.temp%*% sqrtW.temp %*% StdErr %*% (Y - B[,"Bi"])) +
                           (1-pi.a)*crossprod(  StdErr.c %*%dInvLinkdEta.c%*%X.c , R.alpha.inv %*%  StdErr.c %*% (B[,"B.c"]-mu.c))+
                           (pi.a)*crossprod(   StdErr.t %*%dInvLinkdEta.t %*%X.t , R.alpha.inv %*%  StdErr.t %*% (B[,"B.t"]-mu.t)))    
    }
  }
  
  funcoptvarnuisance<-c(scoreEE,SW)
} 
