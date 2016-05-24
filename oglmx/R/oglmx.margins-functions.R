margins.oglmx<-function(object,Vars=NULL,outcomes="All",atmeans=TRUE,AME=FALSE,ascontinuous=FALSE,dummyzero=FALSE,location=NULL){
  # take object of type oglmx, return marginal effects.
  # outcomes specifies which outcomes to calculate marginal effects for, input "All" or a integer vector listing each outcome to consider.
  # atmeans indicates that marginal effects are to be calculated at the means for all variables.
  # dummyzero indicates that marginal effects are to be taken with binary variables set to zero, overrides atmean for relevant variables.
  # in the case that atmean is FALSE location specifies the vector where the marginal effect is to be taken.
  # ascontinuous specifies whether binary variables should be treated as continuous 
  # Vars is a character vector with the name of the variables for which a marginal effect is sought. If NULL all are collected.
  
  if (outcomes[1]=="All"){
    listoutcomes<-c(1:object$NOutcomes)
    outcomes<-object$Outcomes
  } else {
    if (sum(is.na(match(outcomes,object$Outcomes)))>0){
      stop("Invalid Outcome Specified")
    }
    listoutcomes<-c(1:object$NOutcomes)[match(outcomes,object$Outcomes)]
  }
  nameoutcomes<-object$Outcomes[match(outcomes,object$Outcomes)]
  
  link<-object$link
  if (link=="logit"){
    ProbFunc<-function(p){plogis(p)}
    ProbFuncD<-function(p){dlogis(p)}
    ProbFuncDD<-function(p){dlogis(p)*(1-2*plogis(p))}
  }
  if (link=="probit"){
    ProbFunc<-function(p){pnorm(p)}
    ProbFuncD<-function(p){dnorm(p)}
    ProbFuncDD<-function(p){-p*dnorm(p)}
  }
  if (link=="cauchit"){
    ProbFunc<-function(p){pcauchy(p)}
    ProbFuncD<-function(p){dcauchy(p)}
    ProbFuncDD<-function(p){-2*p*dcauchy(p)/(1+p^2)}
  }
  if (link=="loglog"){
    ProbFunc<-function(p){exp(-exp(-p))}
    ProbFuncD<-function(p){exp(-(exp(-p)+p))}
    ProbFuncDD<-function(p){exp(-(exp(-p)+p))*(1+exp(-p))}
  }
  if (link=="cloglog"){
    ProbFunc<-function(p){1-exp(-exp(p))}
    ProbFuncD<-function(p){exp(p-exp(p))}
    ProbFuncDD<-function(p){exp(p-exp(p))*(1-exp(p))}
  }
  
  sdmodfirstderiv<-D(object$sdmodel,"z")
  sdmodsecondderiv<-D(sdmodfirstderiv,"z")
  
  # find the location of the relevant parameters related to the marginal effects.
  if (is.null(Vars)){
    betalocs<-c(1:length(object$allparams$beta))
    if (object$Hetero){
      deltalocs<-c(1:length(object$allparams$delta))   
    } else {
      deltalocs<-NULL
    }
  } else {
    
  }
  
  # first get elements that are in both equations and specify locations
  # then elements that are only in beta
  # then elements that are only in delta
  
  betaboth<-match(object$BothEq$meanandvarLOC,betalocs)[!is.na(match(object$BothEq$meanandvarLOC,betalocs))] # gives the location of beta parameters
  deltaboth<-object$BothEq$meanandvarLOCZ[match(betaboth,object$BothEq$meanandvarLOC)]
  betaonly<-betalocs[is.na(match(betalocs,betaboth)) & names(object$varMeans[[1]])!="(Intercept)" ] #
  deltaonly<-deltalocs[is.na(match(deltalocs,deltaboth)) & names(object$varMeans[[2]])!="(Intercept)"] # 
  no.effects<-length(betaboth)+length(betaonly)+length(deltaonly)
  EffectNames<-object$BothEq$meanandvarNAME
  if (length(betaonly)>0){
    EffectNames<-c(EffectNames,names(object$varMeans[[1]])[betaonly])
  }
  if (length(deltaonly)>0){
    EffectNames<-c(EffectNames,names(object$varMeans[[2]])[deltaonly])
  }
  
  if (atmeans){
    XVar<-object$varMeans[[1]]
    if (dummyzero){XVar[object$varBinary[[1]]]<-0}
    if (object$Hetero){
      ZVar<-object$varMeans[[2]]
      if (dummyzero){ZVar[object$varBinary[[2]]]<-0}
    }
  }
  if (!object$Hetero){ZVar<-as.matrix(1)}
  # need to setup for the case of user specified location.
  
  BinaryVar<-vector("logical",no.effects)
  if (length(betaboth)>0){
    BinaryVar[1:length(betaboth)]<-as.logical(object$varBinary[[1]][betaboth])
  }
  if (length(betaonly)>0){
    BinaryVar[(length(betaboth)+1):(length(betaboth)+length(betaonly))]<-object$varBinary[[1]][betaonly]
  }
  if (length(deltaonly)>0){
    BinaryVar[(length(betaboth)+length(betaonly)+1):(length(betaboth)+length(betaonly)+length(deltaonly))]<-object$varBinary[[2]][deltaonly]
  }
  
  inputdata<-list()
  for (j in 1:length(listoutcomes)){
    inputdata[[j]]<-cbind(c(betaboth,betaonly,rep(NA,length(deltaonly))),c(deltaboth,rep(NA,length(betaonly)),deltaonly),BinaryVar)
  }
  
  # collect terms that are used commonly in expressions for both the marginal effect and the standard error for the 
  # marginal effect.
  bparams<-object$allparams$beta
  dparams<-object$allparams$delta
  tparams<-object$allparams$threshparam
  
  if (!AME){
    Xb<-XVar%*%bparams
    Zd<-ZVar%*%dparams
    Zdinv<-1/eval({z<-Zd;object$sdmodel})
    gDZd<-eval({z<-Zd;sdmodfirstderiv})
  }
  est.beta<-c(1:length(object$Est.Parameters$beta))[object$Est.Parameters$beta]
  est.delta<-c(1:length(object$Est.Parameters$delta))[object$Est.Parameters$delta]
  est.alpha<-c(1:length(object$Est.Parameters$alpha))[object$Est.Parameters$alpha]
              
  atmeans.func<-function(x){
  
    # x is vector in which the
    # first element specifies the location of the relevant variable in the beta vector
    # second element specifies the location of the relevant variable in the delta vector
    # third element element specifies whether the variable is binary
    
    # collect vector of derivatives of function that gives marginal effects: for application of delta method
    # to calculate standard errors.
    # length of vector should match the number of columns in the hessian returned via the call to oglmx
        
    Gderiv<-vector("numeric",ncol(object$hessian))
    # calculate marginal effects for non-binary explanatories
    if (x[3]==0 | ascontinuous){ # | !dummyzero 
      if (!is.na(x[1])){
        value1<- pdfdiff*bparams[x[1]]*-Zdinv
        start<-1
        if (length(est.beta)>0){
          Gderiv[1:length(est.beta)]<-bparams[x[1]]*XVar[est.beta]*pdfDdiff*(Zdinv^2)
          for (b in est.beta){
            if (b==x[1]){
              Gderiv[start]<-Gderiv[start]-Zdinv*pdfdiff
            }
            start<-start+1
          }
        }
        if (length(est.delta)>0){
          Gderiv[(1+length(est.beta)):(length(est.beta)+length(est.delta))]<-bparams[x[1]]*ZVar[est.delta]*(Zdinv^2)*gDZd*(pdfdiff+xpdfDdiff)
          start<-start+length(est.delta)
        }
        if (length(est.alpha)>0){
          for (a in est.alpha){
            if (a==j){
              Gderiv[start]<- Gderiv[start]-bparams[x[1]]*(Zdinv^2)*ProbFuncDD((tparams[a]-Xb)*Zdinv)
            } else if (a==j-1){
              Gderiv[start]<- Gderiv[start]+bparams[x[1]]*(Zdinv^2)*ProbFuncDD((tparams[a]-Xb)*Zdinv)
            }
            start<-start+1
          }
        }
      } else {
        value1<-0
      }
      if (!is.na(x[2])){
        value2<-xpdfdiff*dparams[x[2]]*gDZd*-Zdinv
        start<-1
        if (length(est.beta)>0){
          Gderiv[1:length(est.beta)]<-Gderiv[1:length(est.beta)]+dparams[x[2]]*XVar[est.beta]*(Zdinv^2)*gDZd*(pdfdiff+xpdfDdiff)
          start<-start+length(est.beta)
        }
        if (length(est.delta)>0){
          Gderiv[start:(length(est.beta)+length(est.delta))]<-Gderiv[start:(length(est.beta)+length(est.delta))]+dparams[x[2]]*ZVar[est.delta]*Zdinv*(2*Zdinv*(gDZd^2)*xpdfdiff-eval({z<-Zd;sdmodsecondderiv})*xpdfdiff+Zdinv*(gDZd^2)*xpdfD2diff)
          for (d in est.delta){
            if (d==x[2]){
              Gderiv[start]<-Gderiv[start]-Zdinv*gDZd*xpdfdiff
            }
            start<-start+1
          }
        }
        if (length(est.alpha)>0){
          start<-length(est.beta)+length(est.delta)+1
          for (a in est.alpha){
            if (a==j){
              Gderiv[start]<-Gderiv[start]-dparams[x[2]]*(Zdinv^2)*gDZd*(ProbFuncD((tparams[a]-Xb)*Zdinv)+(tparams[a]-Xb)*Zdinv*ProbFuncDD((tparams[a]-Xb)*Zdinv))
            } else if (a==j-1){
              Gderiv[start]<-Gderiv[start]+dparams[x[2]]*(Zdinv^2)*gDZd*(ProbFuncD((tparams[a]-Xb)*Zdinv)+(tparams[a]-Xb)*Zdinv*ProbFuncDD((tparams[a]-Xb)*Zdinv))
            }
            start<-start+1
          }
        }
      } else {
        value2<-0
      }
      value<-value1+value2
    
    } else if (x[3]==1){  #& dummyzero
      # calculate marginal effects for variables that are binary
      if (!is.na(x[1])){
        NoTreatX<-TreatX<-XVar
        TreatX[x[1]]<-1
        NoTreatX[x[1]]<-0
        XbTreat<-TreatX%*%bparams
        XbNoTreat<-NoTreatX%*%bparams
      } else {
        NoTreatX<-TreatX<-XVar
        XbNoTreat<-XbTreat<-Xb
      }
      if (!is.na(x[2])){
        NoTreatZ<-TreatZ<-ZVar
        TreatZ[x[2]]<-1
        NoTreatZ[x[2]]<-0
        ZdTreat<-TreatZ%*%dparams
        ZdNoTreat<-NoTreatZ%*%dparams
        ZdTreatInv<-1/eval({z<-ZdTreat;object$sdmodel})
        ZdNoTreatInv<-1/eval({z<-ZdNoTreat;object$sdmodel})
      } else {
        ZdTreatInv<-ZdNoTreatInv<-Zdinv
        ZdTreat<-ZdNoTreat<-Zd
        NoTreatZ<-TreatZ<-ZVar
      }
      gDZdTreat<-eval({z<-ZdTreat;sdmodfirstderiv})
      gDZdNoTreat<-eval({z<-ZdNoTreat;sdmodfirstderiv})
      
      if (j==1){
        value<-ProbFunc((tparams[j]-XbTreat)*ZdTreatInv)-ProbFunc((tparams[j]-XbNoTreat)*ZdNoTreatInv)
        Treatpdfdiff<-ProbFuncD((tparams[j]-XbTreat)*ZdTreatInv)
        NoTreatpdfdiff<-ProbFuncD((tparams[j]-XbNoTreat)*ZdNoTreatInv)
        Treatxpdfdiff<-ProbFuncD((tparams[j]-XbTreat)*ZdTreatInv)*(tparams[j]-XbTreat)*ZdTreatInv
        NoTreatxpdfdiff<-ProbFuncD((tparams[j]-XbNoTreat)*ZdNoTreatInv)*(tparams[j]-XbNoTreat)*ZdNoTreatInv
      } else if (j==object$NOutcomes){
        value<- 1-ProbFunc((tparams[j-1]-XbTreat)*ZdTreatInv)-(1-ProbFunc((tparams[j-1]-XbNoTreat)*ZdNoTreatInv))
        Treatpdfdiff <- -ProbFuncD((tparams[j-1]-XbTreat)*ZdTreatInv)
        NoTreatpdfdiff <- -ProbFuncD((tparams[j-1]-XbNoTreat)*ZdNoTreatInv)
        Treatxpdfdiff<- -ProbFuncD((tparams[j-1]-XbTreat)*ZdTreatInv)*(tparams[j-1]-XbTreat)*ZdTreatInv
        NoTreatxpdfdiff<- -ProbFuncD((tparams[j-1]-XbNoTreat)*ZdNoTreatInv)*(tparams[j-1]-XbNoTreat)*ZdNoTreatInv
      } else {
        value<- ProbFunc((tparams[j]-XbTreat)*ZdTreatInv)-ProbFunc((tparams[j-1]-XbTreat)*ZdTreatInv)-(ProbFunc((tparams[j]-XbNoTreat)*ZdNoTreatInv)-ProbFunc((tparams[j-1]-XbNoTreat)*ZdNoTreatInv))
        Treatpdfdiff<-ProbFuncD((tparams[j]-XbTreat)*ZdTreatInv)-ProbFuncD((tparams[j-1]-XbTreat)*ZdTreatInv)
        NoTreatpdfdiff<-ProbFuncD((tparams[j]-XbNoTreat)*ZdNoTreatInv)-ProbFuncD((tparams[j-1]-XbNoTreat)*ZdNoTreatInv)
        Treatxpdfdiff<-ProbFuncD((tparams[j]-XbTreat)*ZdTreatInv)*(tparams[j]-XbTreat)*ZdTreatInv-ProbFuncD((tparams[j-1]-XbTreat)*ZdTreatInv)*(tparams[j-1]-XbTreat)*ZdTreatInv
        NoTreatxpdfdiff<-ProbFuncD((tparams[j]-XbNoTreat)*ZdNoTreatInv)*(tparams[j]-XbNoTreat)*ZdNoTreatInv-ProbFuncD((tparams[j-1]-XbNoTreat)*ZdNoTreatInv)*(tparams[j-1]-XbNoTreat)*ZdNoTreatInv
      }
      
      if (length(est.beta)>0){
        Gderiv[1:length(est.beta)]<- -TreatX[est.beta]*ZdTreatInv*Treatpdfdiff+NoTreatX[est.beta]*ZdNoTreatInv*NoTreatpdfdiff
      }
      
      if (length(est.delta)>0){
        Gderiv[(1+length(est.beta)):(length(est.beta)+length(est.delta))]<- -TreatZ[est.delta]*ZdTreatInv*gDZdTreat*Treatxpdfdiff+NoTreatZ[est.delta]*ZdNoTreatInv*gDZdNoTreat*NoTreatxpdfdiff
      }
      
      if (length(est.alpha)>0){
        start<-length(est.beta)+length(est.delta)+1
        for (a in est.alpha){
          if (a==j){
            Gderiv[start]<- ZdTreatInv*ProbFuncD((tparams[a]-XbTreat)*ZdTreatInv)-ZdNoTreatInv*ProbFuncD((tparams[a]-XbNoTreat)*ZdNoTreatInv)
          } else if (a==j-1){
            Gderiv[start]<- -ZdTreatInv*ProbFuncD((tparams[a]-XbTreat)*ZdTreatInv)+ZdNoTreatInv*ProbFuncD((tparams[a]-XbNoTreat)*ZdNoTreatInv)
          }
          start<-start+1
        }
      }
      
    }
    attr(value,"GDerivatives")<-Gderiv
    return(value)
  }
  
if (AME & is.null(object$modelframes)){
  stop("Model frame required to obtain average marginal effects, reestimate with savemodelframes set to TRUE.")
} else if (AME){
  XVar<-object$modelframes[[1]]
  ZVar<-object$modelframes[[2]]
  Xb<-XVar%*%bparams
  Zd<-ZVar%*%dparams
  Zdinv<-1/eval({z<-object$modelframes[[2]]%*%dparams;object$sdmodel})
  gDZd<-eval({z<-Zd;sdmodfirstderiv})
}

  AME.func<-function(x){
     #x<-inputdata[[1]][1,]
    Zdinv<-as.vector(Zdinv)
    pdfdiff<-as.vector(pdfdiff)
    xpdfdiff<-as.vector(xpdfdiff)
    pdfDdiff<-as.vector(pdfDdiff)
    gDZd<-as.vector(gDZd)
    xpdfDdiff<-as.vector(xpdfDdiff)
    xpdfD2diff<-as.vector(xpdfD2diff)
    
    Gderiv<-matrix(vector("numeric",length(Xb)*ncol(object$hessian)),nrow=length(Xb),ncol=ncol(object$hessian))
    if (x[3]==0 | ascontinuous){
      if (!is.na(x[1])){
        value1<- -bparams[x[1]]*Zdinv*pdfdiff
        start<-1
        if (length(est.beta)>0){
          Gderiv[,1:length(est.beta)]<-bparams[x[1]]*(Zdinv^2)*XVar[,est.beta]*pdfDdiff
          for (b in est.beta){
            if (b==x[1]){
              Gderiv[,start]<-Gderiv[,start]-Zdinv*pdfdiff
            }
            start<-start+1
          }
        }
        if (length(est.delta)>0){
          Gderiv[,(1+length(est.beta)):(length(est.beta)+length(est.delta))]<-bparams[x[1]]*ZVar[,est.delta]*(Zdinv^2)*gDZd*(pdfdiff+xpdfDdiff)
          start<-start+length(est.delta)
        }
        if (length(est.alpha)>0){
          for (a in est.alpha){
            if (a==j){
              Gderiv[,start]<-Gderiv[,start]-bparams[x[1]]*(Zdinv^2)*ProbFuncDD((tparams[a]-Xb)*Zdinv)
            } else if (a==j-1){
              Gderiv[,start]<-Gderiv[,start]+bparams[x[1]]*(Zdinv^2)*ProbFuncDD((tparams[a]-Xb)*Zdinv)
            }
            start<-start+1
          }
        }
      } else {
        value1<-vector("numeric",length(Xb))
      }  
      if (!is.na(x[2])){
        value2<- -dparams[x[2]]*Zdinv*gDZd*xpdfdiff
        start<-1
        if (length(est.beta)>0){
          Gderiv[,1:length(est.beta)]<-Gderiv[,1:length(est.beta)]+dparams[x[2]]*XVar[,est.beta]*(Zdinv^2)*gDZd*(pdfdiff+xpdfDdiff)
          start<-start+length(est.beta)
        }
        if (length(est.delta)>0){
          Gderiv[,start:(length(est.beta)+length(est.delta))]<-Gderiv[,start:(length(est.beta)+length(est.delta))]+dparams[x[2]]*ZVar[,est.delta]*Zdinv*(2*Zdinv*(gDZd^2)*xpdfdiff-as.vector(eval({z<-Zd;sdmodsecondderiv}))*xpdfdiff+Zdinv*(gDZd^2)*xpdfD2diff)
          for (d in est.delta){
            if (d==x[2]){
              Gderiv[,start]<-Gderiv[,start]-Zdinv*gDZd*xpdfdiff
            }
            start<-start+1
          }
        }
        if (length(est.alpha)>0){
          start<-length(est.beta)+length(est.delta)+1
          for (a in est.alpha){
            if (a==j){
              Gderiv[,start]<-Gderiv[,start]-dparams[x[2]]*(Zdinv^2)*gDZd*(ProbFuncD((tparams[a]-Xb)*Zdinv)+(tparams[a]-Xb)*Zdinv*ProbFuncDD((tparams[a]-Xb)*Zdinv))
            } else if (a==j-1){
              Gderiv[,start]<-Gderiv[,start]+dparams[x[2]]*(Zdinv^2)*gDZd*(ProbFuncD((tparams[a]-Xb)*Zdinv)+(tparams[a]-Xb)*Zdinv*ProbFuncDD((tparams[a]-Xb)*Zdinv))
            }
            start<-start+1
          }
        }
      } else {
        value2<-vector("numeric",length(Xb))
      }
      value<-mean(value1+value2)
      Gderiv<-apply(Gderiv,2,mean)
    } else if (x[3]==1){
      stop("Code not written to deal with binary variables when calculating Average Marginal Effect.")
    }
    
    attr(value,"GDerivatives")<-Gderiv
    return(value)
  }
  
  VCOV<-vcov.oglmx(object,tol=1e-40)
  calcSE<-function(x){
    gprime<-attr(x,"GDerivatives")
    outputSE<-(gprime%*%VCOV%*%gprime)^0.5
  }
  collectmfx<-function(x){outputmfx<-x[[1]]}
  
  results<-list()
  count<-1
  for (j in listoutcomes){
    input<-list()
    for (k in 1:nrow(inputdata[[count]])){
      input[[k]]<-inputdata[[count]][k,]
    }
    
    if (j==1){
      frac1<-(tparams[j]-Xb)*Zdinv
      pdfdiff<-ProbFuncD(frac1)
      pdfDdiff<-ProbFuncDD(frac1)
      xpdfdiff<-ProbFuncD(frac1)*frac1
      xpdfDdiff<-ProbFuncDD(frac1)*frac1
      xpdfD2diff<-ProbFuncDD(frac1)*(frac1)^2
    } else if (j==object$NOutcomes){
      frac0<-(tparams[j-1]-Xb)*Zdinv
      pdfdiff<- -ProbFuncD(frac0)
      pdfDdiff<- -ProbFuncDD(frac0)
      xpdfdiff<- -ProbFuncD(frac0)*frac0
      xpdfDdiff<- -ProbFuncDD(frac0)*frac0
      xpdfD2diff<- -ProbFuncDD(frac0)*(frac0)^2
    } else {
      frac1<-(tparams[j]-Xb)*Zdinv
      frac0<-(tparams[j-1]-Xb)*Zdinv
      pdfdiff<-ProbFuncD(frac1)-ProbFuncD(frac0)
      pdfDdiff<-ProbFuncDD(frac1)-ProbFuncDD(frac0)
      xpdfdiff<-ProbFuncD(frac1)*frac1-ProbFuncD(frac0)*frac0
      xpdfDdiff<-ProbFuncDD(frac1)*frac1-ProbFuncDD(frac0)*frac0
      xpdfD2diff<-ProbFuncDD(frac1)*(frac1)^2-ProbFuncDD(frac0)*(frac0)^2
    }
    
    if (!AME){
      output<-lapply(input,atmeans.func)
    }
    if (AME){
      output<-lapply(input,AME.func)
    }
    mfx<-vapply(output,collectmfx,vector("numeric",1))
    mfxSE<-vapply(output,calcSE,vector("numeric",1))
    mfxt<-mfx/mfxSE
    mfxProb<- 2*pnorm( -abs( mfxt))
    results[[count]]<-cbind("Marg. Eff"=mfx,"Std. error"=mfxSE,"t value"=mfxt,"Pr(>|t|)"=mfxProb)
    row.names(results[[count]])<-EffectNames
    count<-count+1
  }
  output<-results
  names(output)<-as.character(nameoutcomes)
  class(output)<-c("oglmx.margins")
  if (!AME){
    marginat<-list(XVar,ZVar)
  } else {
    marginat<-NULL
  }
  attr(output,"Margin Type")<-list(AME=AME,at=marginat)
  class(output)<-c("margins.oglmx")
  output
}

print.margins.oglmx<-function(x, ... ){
  for (m in 1:length(x)){
    cat("Marginal Effects on Pr(Outcome==",names(x)[m],")","\n",sep="")
    if (m==length(x)){
      printCoefmat(x[[m]])
    } else {
      printCoefmat(x[[m]],signif.legend=FALSE)
      cat("------------------------------------","\n")
    }
  }
}