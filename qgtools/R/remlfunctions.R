


##A function to generate a simulated data set
lmm.simudata=function(formula, data = list(),v,b,SimuNum=NULL){
    if (is.null(data))gdata = mixed.data(formula)
    else gdata = mixed.data(formula, data)
    if(is.null(SimuNum))SimuNum=200
    gdata=SimuData(gdata,v=v,b=b,SimuNum=SimuNum)
    YS=gdata$Y
    return(YS)
}

lmm=function(formula,data = list(),method=NULL,ALPHA=NULL){
    if(is.null(method))method=c("minque")
    if(is.null(ALPHA))ALPHA=0.05
    if (is.null(data))gdata = mixed.data(formula)
    else gdata = mixed.data(formula, data)
    #gdata = mixed.data(formula)
    mn=length(method)
    RES=list()
    for(i in 1:mn){
       if(method[i]=="minque")RES[[i]]=mq0(gdata)
       if(method[i]=="reml")RES[[i]]=reml0(gdata)
    }
    mk=length(gdata$U)
    Uk=numeric()
    BigU=NULL
    for(u in 1:mk){
       Uk[u]=ncol(gdata$U[[u]])
       BigU=cbind(BigU,gdata$U[[u]])
    }
    mt=sum(Uk)
    xk=length(gdata$X)
    Xk=numeric()
    X=NULL
    for(i in 1:xk){
       Xk[i]=ncol(gdata$X[[i]])
       X=cbind(X,gdata$X[[i]])
    }
    TraitNum=ncol(gdata$Y)
    for(i in 1:mn){
       Residual=NULL
       Fitted=NULL
       for(j in 1:TraitNum){
           a=0
           re=RES[[i]]$RandomEffect[[j]][,1]
           fe=RES[[i]]$FixedEffect[[j]][,1]
           if(mk>0)a=a+BigU%*%re
           a=a+X%*%fe
           Fitted=cbind(Fitted,a)
           a=gdata$Y[,j]-a
           Residual=cbind(Residual,a)
       }
       colnames(Residual)=colnames(gdata$Y)
       colnames(Fitted)=colnames(gdata$Y)
       RES[[i]]$Residual=Residual
       RES[[i]]$Fitted=Fitted
    }

    names(RES)=method
    RES$ALPHA=ALPHA
    return(RES)
}

lmm.jack=function(formula,data = list(),method=NULL,JacNum=NULL,JacRep=NULL,ALPHA=NULL){
    if(is.null(method))method=c("minque")
    if(is.null(JacNum))JacNum=10
    if(is.null(JacRep))JacRep=1
    if(is.null(ALPHA))ALPHA=0.05
    if (is.null(data))gdata = mixed.data(formula)
    else gdata = mixed.data(formula, data)
    mn=length(method)
    RES=list()
    for(i in 1:mn){
       if(method[i]=="minque")RES[[i]]=genmod.jack(gdata,JacNum=JacNum,JacRep=JacRep,ALPHA=ALPHA)
       if(method[i]=="reml")RES[[i]]=genmod.reml.jack(gdata,JacNum=JacNum,JacRep=JacRep,ALPHA=ALPHA)
    }
    mk=length(gdata$U)
    Uk=numeric()
    BigU=NULL
    for(u in 1:mk){
       Uk[u]=ncol(gdata$U[[u]])
       BigU=cbind(BigU,gdata$U[[u]])
    }
    mt=sum(Uk)
    xk=length(gdata$X)
    Xk=numeric()
    X=NULL
    for(i in 1:xk){
       Xk[i]=ncol(gdata$X[[i]])
       X=cbind(X,gdata$X[[i]])
    }
    TraitNum=ncol(gdata$Y)
    for(i in 1:mn){
       Residual=NULL
       Fitted=NULL
       for(j in 1:TraitNum){
           a=0
           re=RES[[i]]$RandomEffect[[j]][,1]
           fe=RES[[i]]$FixedEffect[[j]][,1]
           if(mk>0)a=a+BigU%*%re
           a=a+X%*%fe
           Fitted=cbind(Fitted,a)
           a=gdata$Y[,j]-a
           Residual=cbind(Residual,a)
       }
       colnames(Residual)=colnames(gdata$Y)
       colnames(Fitted)=colnames(gdata$Y)
       RES[[i]]$Residual=Residual
       RES[[i]]$Fitted=Fitted
    }

    names(RES)=method
    RES$ALPHA=ALPHA
    return(RES)
}

lmm.jack1=function(formula,lmlist=NULL){
    # data = list(),method=NULL,JacNum=NULL,JacRep=NULL,ALPHA=NULL){
    if(is.null(lmlist)){
        method=c("minque")
        JacNum=10
        JacRep=1
        ALPHA=0.05
        gdata = mixed.data(formula)
    }
    else{
      
       if(is.null(lmlist$method))method=c("minque")
       if(is.null(lmlist$JacNum))JacNum=10
       if(is.null(lmlist$JacRep))JacRep=1
       if(is.null(lmlist$ALPHA))ALPHA=0.05
       if (is.null(lmlist$data))gdata = mixed.data(formula)
       else gdata = mixed.data(formula, data=lmlist$data)
    }
    mn=length(method)
    RES=list()
    for(i in 1:mn){
       if(method[i]=="minque")RES[[i]]=genmod.jack(gdata,JacNum=JacNum,JacRep=JacRep,ALPHA=ALPHA)
       if(method[i]=="reml")RES[[i]]=genmod.reml.jack(gdata,JacNum=JacNum,JacRep=JacRep,ALPHA=ALPHA)
    }
    mk=length(gdata$U)
    Uk=numeric()
    BigU=NULL
    for(u in 1:mk){
       Uk[u]=ncol(gdata$U[[u]])
       BigU=cbind(BigU,gdata$U[[u]])
    }
    mt=sum(Uk)
    xk=length(gdata$X)
    Xk=numeric()
    X=NULL
    for(i in 1:xk){
       Xk[i]=ncol(gdata$X[[i]])
       X=cbind(X,gdata$X[[i]])
    }
    TraitNum=ncol(gdata$Y)
    for(i in 1:mn){
       Residual=NULL
       Fitted=NULL
       for(j in 1:TraitNum){
           a=0
           re=RES[[i]]$RandomEffect[[j]][,1]
           fe=RES[[i]]$FixedEffect[[j]][,1]
           if(mk>0)a=a+BigU%*%re
           a=a+X%*%fe
           Fitted=cbind(Fitted,a)
           a=gdata$Y[,j]-a
           Residual=cbind(Residual,a)
       }
       colnames(Residual)=colnames(gdata$Y)
       colnames(Fitted)=colnames(gdata$Y)
       RES[[i]]$Residual=Residual
       RES[[i]]$Fitted=Fitted
    }

    names(RES)=method
    RES$ALPHA=ALPHA
    return(RES)
}




lmm.simu=function(formula,method=NULL,ALPHA=NULL){
    if(is.null(method))method=c("minque")
    if(is.null(ALPHA))ALPHA=0.05
    gdata = mixed.data(formula)
    mn=length(method)
    RES=list()
    for(i in 1:mn){
       if(method[i]=="minque")RES[[i]]=mq0(gdata)
       if(method[i]=="reml")RES[[i]]=reml0(gdata)
    }
    for(i in 1:mn){
       res=RES[[i]]
       V=NULL
       Var=res$Var
       SimuNum=length(Var)
       ml=length(Var[[1]][,1])
       P=numeric(ml)
       vm=numeric(ml)
       for(j in 1:SimuNum){
           vm=vm+Var[[j]][,1]
           V=cbind(V,Var[[j]][,1])
           id=which(Var[[j]][,4]<=ALPHA)
           P[id]=P[id]+1
       }
       SE=numeric(ml)
       for(j in 1:ml)SE[j]=sqrt(var(V[j,])/SimuNum)
       vm=vm/SimuNum
       P=P/SimuNum
       a=data.frame(vm,SE,P)
       colnames(a)=c("Estimate","SE","Power")
       rownames(a)=rownames(Var[[1]])
       RES[[i]]=a
    }
    names(RES)=method
    RES$ALPHA=ALPHA
    return(RES)
}

lmm.simu.jack=function(formula,method=NULL,JacNum=NULL,JacRep=NULL,ALPHA=NULL){
    if(is.null(method))method=c("minque")
    if(is.null(JacNum))JacNum=10
    if(is.null(JacRep))JacRep=1
    if(is.null(ALPHA))ALPHA=0.05
    gdata = mixed.data(formula)
    mn=length(method)
    RES=list()
    for(i in 1:mn){
       if(method[i]=="minque")RES[[i]]=genmod.jack(gdata,JacNum=JacNum,JacRep=JacRep,ALPHA=ALPHA)
       if(method[i]=="reml")RES[[i]]=genmod.reml.jack(gdata,JacNum=JacNum,JacRep=JacRep,ALPHA=ALPHA)
    }
    for(i in 1:mn){
       res=RES[[i]]
       V=NULL
       Var=res$Var
       SimuNum=length(Var)
       ml=length(Var[[1]][,1])
       P=numeric(ml)
       vm=numeric(ml)
       for(j in 1:SimuNum){
           vm=vm+Var[[j]][,1]
           V=cbind(V,Var[[j]][,1])
           id=which(Var[[j]][,3]<=ALPHA)
           P[id]=P[id]+1
       }
       SE=numeric(ml)
       for(j in 1:ml)SE[j]=sqrt(var(V[j,])/SimuNum)
       vm=vm/SimuNum
       P=P/SimuNum
       a=data.frame(vm,SE,P)
       colnames(a)=c("Estimate","SE","Power")
       rownames(a)=rownames(Var[[1]])
       RES[[i]]=a
    }
    names(RES)=method
    RES$ALPHA=ALPHA
    return(RES)
}

  


      

lmm.reml.simu=function(formula, data = list(),v0,b0,SimuNum=NULL,JacNum=NULL,JacRep=NULL,ALPHA=NULL){
    if (is.null(data))gdata = mixed.data(formula)
    else gdata = mixed.data(formula, data)
    if(is.null(SimuNum))SimuNum=200
    if(is.null(JacNum))JacNum=10
    if(is.null(JacRep))JacRep=1
    if(is.null(ALPHA))ALPHA=0.05
    result=genmod.reml.simu(gdata,v=v0,b=b0,SimuNum=SimuNum,JacNum=JacNum,JacRep=JacRep,ALPHA=ALPHA)
    return(result)
}    

genmod.reml.simu=function(gdata,v,b,SimuNum=NULL,JacNum=NULL,JacRep=NULL,ALPHA=NULL,...){
   if(is.null(SimuNum))SimuNum=200
   if(is.null(ALPHA))ALPHA=0.05
   gdata=SimuData(gdata,v=v,b=b,SimuNum=SimuNum)
   if(is.null(JacNum))JacNum=10
   if(is.null(JacRep))JacRep=1
   
   jac=genmod.reml.jack(gdata,JacNum=JacNum,JacRep=JacRep,LUP=1)
   ml=length(v)
   Var=jac$Var
   P=numeric(ml)
   vm=numeric(ml)
   for(i in 1:SimuNum){
      vm=vm+Var[[i]][,1]
      id=which(Var[[i]][,3]<=ALPHA)
      P[id]=P[id]+1
   }
   vm=vm/SimuNum
   P=P/SimuNum
   bias=vm-v
   a=data.frame(v,vm,bias,P)
   colnames(a)=c("True","Estimate","Bias","Power")
   rownames(a)=rownames(Var[[1]])
   res=list(P=a,ALPHA=ALPHA)
   return(res)
}




lmm.reml=function(formula,data=list(),criterion=NULL){
   if(is.null(criterion))criterion=1e-3
   if(is.null(data))gdata=mixed.data(formula)
   else gdata=mixed.data(formula,data)
   res=reml0(gdata,criterion)
   #res$FixedEffect=res$FixEffect
   return(res)
}

lmm.reml.jack=function(formula,data=list(),criterion=NULL,JacNum=NULL,JacRep=NULL,ALPHA=NULL){
   if(is.null(criterion))criterion=1e-3
   if(is.null(JacNum))JacNum=10
   if(is.null(JacRep))JacRep=1
   if(is.null(ALPHA))ALPHA=0.05
   if(is.null(data))gdata=mixed.data(formula)
   else gdata=mixed.data(formula,data)
   res=genmod.reml.jack(gdata,criterion=criterion,JacNum=JacNum,JacRep=JacRep,ALPHA=ALPHA)
   #res$FixedEffect=res$FixEffect
   return(res)
}

lmm.mq=function(formula,data=list()){
   if(is.null(data))gdata=mixed.data(formula)
   else gdata=mixed.data(formula,data)
   res=mq0(gdata)
   #res$FixedEffect=res$FixEffect
   return(res)
}

#gdata$Y

get_HR=function(gdata,Qa,v){
   mk=length(gdata$U)
   h=NULL
   dim(t(gdata$U[[1]]))
   n=ncol(Qa)  
   #V=Diagonal(n)*v[mk+1]
   ##i=1
   #for(i in 1:mk){
   #    vv=gdata$U[[i]]%*%t(gdata$U[[i]])*v[i]
   #    V=V+vv
   #}
   #Q=Qa%*%V%*%Qa
   Q=Qa
 
   for(i in 1:mk){
      m=t(gdata$U[[i]])%*%Q%*%gdata$U[[i]]
      m=m*(v[i])^2
      #m=svdinv(m)
      if(i==1)h=diag(m)
      else h=c(h,diag(m))
   }
   summary(h)
   h=sqrt(h)
   return(h)
}
genmod.mq=function(gdata){
   res0=genmod(gdata)
   res=genmod(gdata,au=res0$Var[,1])
   ##Statistical tests for variance components
   h=diag(svdinv(res$HV))*2
   Chi_sq=(res0$Var[,1])^2/h
   P_value=(1-pchisq(Chi_sq,1))/2
   Var=data.frame(res0$Var,sqrt(h),Chi_sq,P_value)
   colnames(Var)=c("Est","SE","Chi_sq","P_value")
   res$Var=Var
   
   ##Statistical tests for fixed effects
   #res0$FixEffect[,1]
   h=diag(svdinv(res$HF[[1]]))
   h=sqrt(h)
   z=res0$FixEffect[,1]/h
   P_value=1-pchisq(z^2,1)
   attributes(res)
   FixEffect=data.frame(res0$FixEffect,h,z,P_value)
   colnames(FixEffect)=c("Est","SE","z_value","P_value")
   res$FixEffect=FixEffect
   
   ##Statistical tests for random effect  
   mk=length(gdata$U)
   if(mk>0){
   
     h=get_HR(gdata,res$Qa,res0$Var[,1])
     z=res0$RandomEffect[,1]/h
     res0$RandomEffect[,1][1]/h[1]
   
     P_value=1-pchisq(z^2,1)
     RandomEffect=data.frame(res0$RandomEffect,h,z,P_value)
     colnames(RandomEffect)=c("Pre","SE","z_value","P_value")
     res$RandomEffect=RandomEffect
   }
  
   return(res)
}

genmod.reml=function(gdata,criterion=NULL){

   if(is.null(criterion))DIFF=1e-3
   else DIFF=criterion
   d=100
   res=genmod(gdata)
   mk=length(gdata$U)
   au0=res$Var[,1]
   id=which(au0<=0)
   au0[id]=0
   its=1
   while(d>DIFF&&its<10){
      
      res=genmod(gdata,au=au0)
      (v1=res$Var[,1])
      id=which(v1<=0)
      v1[id]=0
      #print(v1)
      #cat("\n")
      
      d=sum(au0-v1)^2
      d=sqrt(d)
      d=d/sum(au0)
      au0=v1
      #print(d)
      #cat("\n")
      its=its+1
   }
   ##Statistical tests for variance components
   h=diag(svdinv(res$HV))*2
   Chi_sq=(res$Var[,1])^2/h
   P_value=(1-pchisq(Chi_sq,1))/2
   Var=data.frame(res$Var,sqrt(h),Chi_sq,P_value)
   colnames(Var)=c("Est","SE","Chi_sq","P_value")
   res$Var=Var
   
   ##Statistical tests for fixed effects
   res$FixEffect[,1]
   h=diag(svdinv(res$HF[[1]]))
   h=sqrt(h)
   z=res$FixEffect[,1]/h
   P_value=1-pchisq(z^2,1)
   attributes(res)
   FixEffect=data.frame(res$FixEffect,h,z,P_value)
   colnames(FixEffect)=c("Est","SE","z_value","P_value")
   res$FixEffect=FixEffect
   
   ##Statistical tests for random effect  
   mk=length(gdata$U)
   if(mk>0){
     h=get_HR(gdata,res$Qa,res$Var[,1])
     z=res$RandomEffect[,1]/h
     res$RandomEffect[,1][1]/h[1]
     P_value=1-pchisq(z^2,1)
     RandomEffect=data.frame(res$RandomEffect,h,z,P_value)
     colnames(RandomEffect)=c("Pre","SE","z_value","P_value")
     res$RandomEffect=RandomEffect
   }
  
   return(res)
}
##This function is used without jackknife process
reml0=function(gdata,criterion=NULL){
   if(is.null(criterion))criterion=1e-3
   Y0=gdata$Y
   TraitNum=ncol(Y0)
   TraitNames=colnames(Y0)
   Var=list()
   RE=NULL
   FE=list()
   mk=length(gdata$U)
   for(i in 1:TraitNum){
      gdata$Y=as.matrix(Y0[,i])
      res=genmod.reml(gdata,criterion)
      Var[[i]]=res$Var
      if(mk>0)RE[[i]]=res$RandomEffect
      #else RE[[i]]=NULL
      FE[[i]]=res$FixEffect
   }
   
   names(Var)=TraitNames
   if(!is.null(RE))names(RE)=TraitNames
   names(FE)=TraitNames
   res=list(Var=Var,FixedEffect=FE,RandomEffect=RE)
   gdata$Y=Y0
   return(res)
}


##This function is used for jackknife process
reml=function(gdata,criterion=NULL){
   if(is.null(criterion))criterion=1e-3
   Y0=as.matrix(gdata$Y)
   TraitNum=ncol(Y0)
   TraitNames=colnames(Y0)
   Var=NULL
   RE=NULL
   FE=NULL
   mk=length(gdata$U)
   for(i in 1:TraitNum){
      gdata$Y=as.matrix(Y0[,i])
      res=genmod.reml(gdata,criterion)
      Var=cbind(Var,res$Var[,1])
      if(mk>0)RE=cbind(RE,res$RandomEffect[,1])
      FE=cbind(FE,res$FixEffect[,1])
   }
   res$FixEffect
   Var=as.matrix(Var)
   colnames(Var)=TraitNames
   rownames(Var)=rownames(res$Var)
   if(mk>0){
     RE=as.matrix(RE)
     colnames(RE)=TraitNames
     rownames(RE)=rownames(res$RandomEffect)
   }
   FE=as.matrix(FE)
   colnames(FE)=TraitNames
   rownames(FE)=rownames(res$FixEffect)
   res=list(Var=Var,FixedEffect=FE,RandomEffect=RE)
   gdata$Y=Y0
   return(res)
}

#A=NULL
#A[[1]]=c(1:2)
##This MINQUE function is used without jackknife process
mq0=function(gdata){
   Y0=as.matrix(gdata$Y)
   #Y0=as.matrix(Y0)
   TraitNum=ncol(Y0)
   TraitNames=colnames(Y0)
   Var=list()
   RE=NULL
   FE=list()
   #i=1
   mk=length(gdata$U)
   for(i in 1:TraitNum){
      gdata$Y=as.matrix(Y0[,i])
      res=genmod.mq(gdata)
      Var[[i]]=res$Var
      if(mk>0)RE[[i]]=res$RandomEffect
      #else RE[[i]]=NULL
      FE[[i]]=res$FixEffect
   }
   #res$FixEffect
   names(Var)=TraitNames
   if(!is.null(RE))names(RE)=TraitNames
   names(FE)=TraitNames
   res=list(Var=Var,FixedEffect=FE,RandomEffect=RE)
   gdata$Y=Y0
   return(res)
}


##This MINQUE function is used with jackknife
mq=function(gdata){
   Y0=gdata$Y
   TraitNum=ncol(Y0)
   TraitNames=colnames(Y0)
   Var=NULL
   RE=NULL
   FE=NULL
   for(i in 1:TraitNum){
      gdata$Y=as.matrix(Y0[,i])
      res=genmod(gdata)
      Var=cbind(Var,res$Var[,1])
      RE=cbind(RE,res$RandomEffect[,1])
      FE=cbind(FE,res$FixEffect[,1])
   }
   res$FixEffect
   Var=as.matrix(Var)
   colnames(Var)=TraitNames
   rownames(Var)=rownames(res$Var)
   RE=as.matrix(RE)
   colnames(RE)=TraitNames
   rownames(RE)=rownames(res$RandomEffect)
   FE=as.matrix(FE)
   colnames(FE)=TraitNames
   rownames(FE)=rownames(res$FixEffect)
   res=list(Var=Var,FixedEffect=FE,RandomEffect=RE)
   gdata$Y=Y0
   return(res)
}
   

genmod.reml.jack=function(gdata,criterion=NULL,JacNum=NULL,JacRep=NULL,ALPHA=NULL,au=NULL,LUP=NULL){
    if(is.null(criterion))criterion=1e-3
    if(is.null(JacNum))JacNum=10
    if(is.null(JacRep))JacRep=1
    if(is.null(ALPHA))ALPHA=0.05
    if(is.null(au))au=NULL
    if(is.null(LUP))LUP=NULL
    a=1-(1-ALPHA)^(1/(JacNum-2))
    t=qt(1-a/2,JacNum-1)
    Y=as.matrix(gdata$Y)
    TraitNum=ncol(Y)
    TraitNames=colnames(Y)
    n=nrow(Y)
    #colnames(Y)
    mk=length(gdata$U)
    xk=length(gdata$X)
    JB=rep(1:JacNum,length=n)
    gdata0=gdata
    ml=mk+1
    nx=0
    for(i in 1:xk)nx=nx+ncol(gdata0$X[[i]])
    res0=reml(gdata0,criterion)
    V0=res0$Var
    B0=res0$FixedEffect
    P0=res0$RandomEffect
    if(mk>0){
      mt=nrow(P0)
      Str=rownames(P0)
   }
   VNames=rownames(V0)
   XNames=rownames(B0)

   ##############################################################################  
   RES0=NULL
   VPower=numeric(ml)
   V0=matrix(0,TraitNum,ml)
   if(mk>0)P0=matrix(0,TraitNum,mt)
   #colnames(P0)
   Vp=matrix(0,TraitNum,ml)
   Cv=matrix(0,TraitNum,ml)
   Cv1=matrix(0,TraitNum,ml)
   Cv2=matrix(0,TraitNum,ml)
   Vm=matrix(0,TraitNum,ml)

   VCp=matrix(0,TraitNum,ml)
   VCv=matrix(0,TraitNum,ml)
   VCv1=matrix(0,TraitNum,ml)
   VCv2=matrix(0,TraitNum,ml)
   VCm=matrix(0,TraitNum,ml)
   
   if(mk>0){
      Pp=matrix(0,TraitNum,mt)
      Cp=matrix(0,TraitNum,mt)
      Cp1=matrix(0,TraitNum,mt)
      Cp2=matrix(0,TraitNum,mt)
      Pm=matrix(0,TraitNum,mt)
   }
   
   Cb1=matrix(0,TraitNum,nx)
   Cb2=matrix(0,TraitNum,nx)
   Cb=matrix(0,TraitNum,nx)
   Bm=matrix(0,TraitNum,nx)
   Bp=matrix(0,TraitNum,nx)
 
   
   VJR=matrix(0,nrow=ml,ncol=JacRep)
   JacVar=matrix(0,nrow=JacNum,ncol=ml)
   if(mk>0)JacPre=matrix(0,nrow=JacNum,ncol=mt)
   JacB=matrix(0,nrow=JacNum,ncol=nx)

   
   JV=matrix(0,nrow=ml,ncol=JacNum)
   if(mk>0)JP=matrix(0,nrow=mt,ncol=JacNum)
   #JB=matrix(0,nrow=nx,ncol=JacNum)

   df=JacNum-1

   JACVAR=matrix(0,TraitNum,ml)
   if(mk>0)JACPRE=matrix(0,TraitNum,mt)
   JACB=matrix(0,TraitNum,nx)

   JACKPRE=list()
   for(k in 1:JacRep){
      #k=1
      JB=sample(JB)
      U1=NULL
      X1=list()
      JacRes=NULL
      JAC1=NULL
      JAC2=NULL
      JAC3=NULL
      jackpre=list()

      for(i in 1:JacNum){
         #i=1
         index=which(JB==i)
         if(mk>0){
            for(j in 1:mk){
              m=gdata0$U[[j]]
              m=m[-index,]
              U1[[j]]=m
            }
         }
         gdata$U=U1
         names(gdata$U)=names(gdata0$U)
         for(j in 1:xk){
           m=gdata0$X[[j]]
           m=as.matrix(m[-index,])
           X1[[j]]=m
           colnames(X1[[j]])=colnames(gdata0$X[[j]])

         }
         gdata$X=X1
         #colnames(gdata$X[[1]])
         YD=as.matrix(Y[-index,])
         colnames(YD)=colnames(Y)
         gdata$Y=YD
         res=reml(gdata,criterion)
         JACVAR=t(res$Var)
         JACB=t(res$FixedEffect)
         if(mk>0){
            JACPRE=t(res$RandomEffect)
            jackpre[[i]]=JACPRE
         }
         JAC1[[i]]=as.matrix(JACVAR)
         if(mk>0)JAC2[[i]]=as.matrix(JACPRE)
         #if(nx==1)JACB
         JAC3[[i]]=as.matrix(JACB)
         

      }
      if(mk>0)JACKPRE[[k]]=jackpre
      for(s in 1:TraitNum){
          #s=1
          for(j in 1:JacNum){
             if(mk==0&&TraitNum==1)JacVar[j,]=JAC1[[j]]
             else JacVar[j,]=JAC1[[j]][s,]
             if(mk>0)JacPre[j,]=JAC2[[j]][s,]
             if(nx==1&&TraitNum==1)JacB[j,]=JAC3[[j]]
             else JacB[j,]=JAC3[[j]][s,]

          }
          for(j in 1:JacNum){
             index=which(JacVar[j,]<0)
             JacVar[j,index]=0
          }
          ####Calculate the variance components and their proportions from the jackknife estimates
          vt=as.vector(apply(JacVar,1,sum))
          for(i in 1:ml){
             m=mean(JacVar[,i])
             se=sqrt(var(JacVar[,i]))
             if(m<=0){
               m=0
               t1=0
               p=1.0
             }
             else if(m>0){
               t1=m/se
               p=pt(-abs(t1),df)
               p=1-(1-p)^(JacNum-2)
             }
             Vp[s,i]=Vp[s,i]+p
             Vm[s,i]=Vm[s,i]+m
             Cv[s,i]=Cv[s,i]+se
             m=mean(JacVar[,i]/vt)
             se=sqrt(var(JacVar[,i]/vt))
             if(m<=0){
               m=0
               t1=0
               p=1.0
             }
             else if(m>0){
               t1=m/se
               p=pt(-abs(t1),df)
               p=1-(1-p)^(JacNum-2)
             }
             VCp[s,i]=VCp[s,i]+p
             VCm[s,i]=VCm[s,i]+m
             VCv[s,i]=VCv[s,i]+se
             
          }
          #if(is.null(SIMUYES)){ 
          if(mk>0){
             for(i in 1:mt){
               m=mean(JacPre[,i])
               se=sqrt(var(JacPre[,i]))
               t1=m/se
             
               p=pt(-abs(t1),df)
               if(se<=0.0000001)p=1.0
               #p=a$p.value
               p=1-(1-p)^(JacNum-2)
               Pp[s,i]=Pp[s,i]+p
               Pm[s,i]=Pm[s,i]+m
               Cp[s,i]=Cp[s,i]+se
             }
          }
          #}
          for(i in 1:nx){
             m=mean(JacB[,i])
             se=sqrt(var(JacB[,i]))
             t1=m/se
             p=pt(-abs(t1),df)
             if(se<=0.0000001)p=1.0
             #p=a$p.value
             p=1-(1-p)^(JacNum-2)
             Bp[s,i]=Bp[s,i]+p
             Bm[s,i]=Bm[s,i]+m
             Cb[s,i]=Cb[s,i]+se
          }

      }
   }
   ### Variance components
   Vm=Vm/JacRep
   Vp=Vp/JacRep
   Cv=Cv/JacRep
   Cv1=Vm-t*Cv
   Cv2=Vm+t*Cv

   ### Proportional variance components
   VCm=VCm/JacRep
   VCp=VCp/JacRep
   VCv=VCv/JacRep
   VCv1=VCm-t*VCv
   VCv2=VCm+t*VCv
   #if(is.null(SIMUYES)){
   if(mk>0){
     ##Random effects   
     Pm=Pm/JacRep
     Pp=Pp/JacRep
     Cp=Cp/JacRep
     Cp1=Pm-t*Cp
     Cp2=Pm+t*Cp
   }
   #}
   ##Fixed effwaects
   Bm=Bm/JacRep
   Bp=Bp/JacRep
   Cb=Cb/JacRep
   Cb1=Bm-t*Cb
   Cb2=Bm+t*Cb
   
   ##

   VAR=NULL
   PRE=NULL
   VC=NULL
   B=NULL
   CL=c(ALPHA/2,1-ALPHA/2)*100
   CL2=c("%LL","%UL")
   CL=paste(CL,CL2,sep="")
   CNames=c("Estimate","SE","PValue",CL)

   #CNames=c("Estimate","SE","PValue","2.5%CL","97.5%CU")
   for(i in 1: TraitNum){
      VAR[[i]]=data.frame(Vm[i,],Cv[i,],Vp[i,],Cv1[i,],Cv2[i,])
      VC[[i]]=data.frame(VCm[i,],VCv[i,],VCp[i,],VCv1[i,],VCv2[i,])
      if(mk>0)PRE[[i]]=data.frame(Pm[i,],Cp[i,],Pp[i,],Cp1[i,],Cp2[i,])
      B[[i]]=data.frame(Bm[i,],Cb[i,],Bp[i,],Cb1[i,],Cb2[i,])

      rownames(VAR[[i]])=VNames
      colnames(VAR[[i]])=CNames
      rownames(VC[[i]])=paste(VNames,"/VP",sep="")
      colnames(VC[[i]])=CNames
      #if(is.null(SIMUYES)){
        if(mk>0){
           rownames(PRE[[i]])=Str
           Rnames=CNames
           Rnames[1]="Pre"
           colnames(PRE[[i]])=Rnames
        }
      #}
      colnames(B[[i]])=CNames
      rownames(B[[i]])=XNames

   }
   
   names(VAR)=TraitNames
   names(VC)=TraitNames
   if(mk>0)names(PRE)=TraitNames
   names(B)=TraitNames
   ##a=list(Var=VAR,PVar=VC,RandomEffect=PRE,FixedEffect=B,Method=Method,JacNum=JacNum,JacRep=JacRep,VNames=VNames,Y=Y,X=X,U=U)
   #if(is.null(SIMUYES))a=list(Var=VAR,PVar=VC,RandomEffect=PRE,FixedEffect=B,JackPre=JACKPRE)##,Method=Method,JacNum=JacNum,JacRep=JacRep,VNames=VNames,Y=Y,X=X,U=U)
   #else a=list(Var=VAR,PVar=VC,FixedEffect=B)
   if(is.null(LUP))Prediction="AUP"
   else Prediction="LUP"
   #a=list(Var=VAR,PVar=VC,RandomEffect=PRE,FixedEffect=B,JackPre=JACKPRE,Prediction=Prediction)
   #a=list(Var=VAR,PVar=VC,RandomEffect=PRE,FixedEffect=B,Prediction=Prediction)
   a=list(Var=VAR,PVar=VC,RandomEffect=PRE,FixedEffect=B)
   #class(a)="jackknife.result"
   #a$call=match.call()
   #a$Var
   return(a)
}

lmm.mq.perm <-
function(formula,data=list(),PermNum=NULL){
   if(is.null(PermNum))PermNum=500
   if(is.null(data))gdata=mixed.data(formula)
   else gdata=mixed.data(formula,data)
   res=genmod.perm(gdata=gdata,au=NULL,PermNum=PermNum)
   return(res)
}

lmm.reml.perm <-
function(formula,data=list(),PermNum=NULL){
   if(is.null(PermNum))PermNum=500
   ##if(is.null(JacRep))JacRep=1
   if(is.null(data))gdata=mixed.data(formula)
   else gdata=mixed.data(formula,data)
   res=reml.perm(gdata=gdata,au=NULL,PermNum=PermNum,LUP=NULL)
   return(res)
}



lmm.perm <-function(formula,data=list(),method=NULL,PermNum=NULL){
   if(is.null(PermNum))PermNum=100
   if(is.null(method))method=c("minque")
   if(is.null(data))gdata=mixed.data(formula)
   else gdata=mixed.data(formula,data)
   mn=length(method)
   RES=list()
   for(i in 1:mn){
       if(method[i]=="minque")RES[[i]]=genmod.perm(gdata,au=NULL,PermNum=PermNum)
       if(method[i]=="reml")RES[[i]]=reml.perm(gdata,au=NULL,PermNum=PermNum)
   }
   names(RES)=method
   return(RES)
}


lmm.ga=function(formula,data=list()){
     if(is.null(data))gdata=mixed.data(formula)
     else gdata=mixed.data(formula=formula,data=data)
     res=genmod.mq(gdata)
}

lmm.gls=function(formula,data=list(),type=NULL,subject=NULL,repeated=NULL,group=NULL){
    if(is.null(au))au=NULL
    if(is.null(data))data=NULL
    gdata=mixed.data(formula,data)
    gdata$au=au
    gdata$data=data
    if(!is.null(subject)){
       if(is.null(type))type="cs"
       #subx=model.frame(subject,data=data)
       #sub.id=CFactor(subx,v=c(1:ncol(subx)))
       aa=getglsid(subject,repeated,data)
       rep.id=aa$rep.id
       sub.id=aa$sub.id
       gdata$type=type
       gdata$rep.id=rep.id
       gdata$sub.id=sub.id
   }
   lmm.gls.eval(gdata)
}
##solve
lmm.gls.eval=function(gdata){
    type=gdata$type
    repid=gdata$rep.id
    subid=gdata$sub.id
    Ur=IndMatrix(subid)
    #dim(Ur)
    res0=gls.mq0(gdata,au=NULL,R=NULL)
    resid0=res0$Residual[,1]
    v=res0$Var[[1]][,1]
    ve=v[length(v)]
    rr=0.50
    e=1
    while(rr<=0.9999&e<=20){
        e=e+1
        mat=getrmmat(resid0,subid,repid)
        rc=ncol(mat)
        rmat=matrix(0,rc,rc)
        if(type=="cs"){
           r=cs_r(mat)
           for(i in 1:rc){
              for(j in 1:rc){
                 if(i==j)rmat[i,j]=ve
                 else rmat[i,j]=ve*r
              }
           }
        }
   
        if(type=="ar1"){
           r=ar1_r(mat)
           for(i in 1:rc){
              for(j in 1:rc){
                 rmat[i,j]=ve*r^abs(i-j)
              }
           }
        }
        if(type=="un"){ok=complete.cases(mat);rmat=cov(mat[ok,])}
        if(type=="het")rmat=un_v(mat)
        V=getR(Ur,rmat,subid,repid)
        res1=gls.mq0(gdata,au=v,R=V)
        resid1=res1$Residual[,1]
        rr=cor(resid0,resid1)
        resid0=resid1
        v=res1$Var[[1]][,1]
        ve=v[length(v)]
    }
    res1$type=type
    res1$rmat=rmat
    res1$rr=rr
    return(res1)
}

lmm.gls.jack=function(formula,data=list(),type=NULL,subject=NULL,repeated=NULL,JacNum=NULL,JacRep=NULL,ALPHA=NULL){
   if(is.null(JacNum))JacNum=10
   if(is.null(JacRep))JacRep=1
   if(is.null(ALPHA))ALPHA=0.05
   if(is.null(data))data=NULL
   if(is.null(data))gdata=mixed.data(formula)
   else gdata=mixed.data(formula,data)
   if(!is.null(subject)){
      aa=getglsid(subject,repeated,data)
      gdata$rep.id=aa$rep.id
      gdata$sub.id=aa$sub.id
      if(is.null(type))type="cs"
      gdata$type=type
   }
   res=gls.jack(gdata,JacNum=JacNum,JacRep=JacRep,ALPHA=ALPHA)
   return(res)
}


gls.jack=function(gdata,JacNum=NULL,JacRep=NULL,ALPHA=NULL,LUP=NULL){
    if(is.null(JacNum))JacNum=10
    if(is.null(JacRep))JacRep=1
    if(is.null(ALPHA))ALPHA=0.05
    #if(is.null(au))
    au=NULL
    if(is.null(LUP))LUP=NULL
    a=1-(1-ALPHA)^(1/(JacNum-2))
    t=qt(1-a/2,JacNum-1)
    Y=gdata$Y
    TraitNum=ncol(Y)
    TraitNames=colnames(Y)
    n=nrow(Y)
    #colnames(Y)
    mk=length(gdata$U)
    xk=length(gdata$X)
    gdata0=gdata
    rep.id=gdata0$rep.id
    sub.id=gdata0$sub.id
    #if(!is.null(gdata$sub.id))JB=rep(1:JacNum,length=length(unique(sub.id)))
    #else 
    JB=rep(1:JacNum,length=n)

    ml=mk+1
    nx=0
    for(i in 1:xk)nx=nx+ncol(gdata0$X[[i]])
    res0=lmm.gls.eval(gdata0)
    V0=res0$Var[[1]]
    B0=res0$FixedEffect[[1]]
    P0=res0$RandomEffect[[1]]
    if(mk>0){
       mt=nrow(P0)
       Str=rownames(P0)
    }
    VNames=rownames(V0)
    XNames=rownames(B0)

    ##############################################################################  
    RES0=NULL
    VPower=numeric(ml)
    V0=matrix(0,TraitNum,ml)
    if(mk>0)P0=matrix(0,TraitNum,mt)
    #colnames(P0)
    Vp=matrix(0,TraitNum,ml)
    Cv=matrix(0,TraitNum,ml)
    Cv1=matrix(0,TraitNum,ml)
    Cv2=matrix(0,TraitNum,ml)
    Vm=matrix(0,TraitNum,ml)

    VCp=matrix(0,TraitNum,ml)
    VCv=matrix(0,TraitNum,ml)
    VCv1=matrix(0,TraitNum,ml)
    VCv2=matrix(0,TraitNum,ml)
    VCm=matrix(0,TraitNum,ml)
   
    if(mk>0){
       Pp=matrix(0,TraitNum,mt)
       Cp=matrix(0,TraitNum,mt)
       Cp1=matrix(0,TraitNum,mt)
       Cp2=matrix(0,TraitNum,mt)
       Pm=matrix(0,TraitNum,mt)
    }
   
    Cb1=matrix(0,TraitNum,nx)
    Cb2=matrix(0,TraitNum,nx)
    Cb=matrix(0,TraitNum,nx)
    Bm=matrix(0,TraitNum,nx)
    Bp=matrix(0,TraitNum,nx)
 
   
    VJR=matrix(0,nrow=ml,ncol=JacRep)
    JacVar=matrix(0,nrow=JacNum,ncol=ml)
    if(mk>0)JacPre=matrix(0,nrow=JacNum,ncol=mt)
    JacB=matrix(0,nrow=JacNum,ncol=nx)

   
    JV=matrix(0,nrow=ml,ncol=JacNum)
    if(mk>0)JP=matrix(0,nrow=mt,ncol=JacNum)
    #JB=matrix(0,nrow=nx,ncol=JacNum)

    df=JacNum-1

    JACVAR=matrix(0,TraitNum,ml)
    if(mk>0)JACPRE=matrix(0,TraitNum,mt)
    JACB=matrix(0,TraitNum,nx)

    JACKPRE=list()
    for(k in 1:JacRep){
      #k=1
      JB=sample(JB)
      U1=NULL
      X1=list()
      JacRes=NULL
      JAC1=NULL
      JAC2=NULL
      JAC3=NULL
      jackpre=list()

      for(i in 1:JacNum){
         #i=2
         index=which(JB==i)
         #indx=getglsindex(sub.id,rep.id,index)
         #index=indx
         #rep.id1=rep.id[-index]
         #sub.id1=sub.id[-index]
         #ok=order(sub.id1,rep.id1)
         gdata$rep.id=gdata0$rep.id[-index]
         gdata$sub.id=gdata0$sub.id[-index]
         #data.frame(gdata$sub.id,gdata$rep.id)
         if(mk>0){
            for(j in 1:mk){
              m=gdata0$U[[j]]
              m=m[-index,]
              #m=m[ok,]
              U1[[j]]=m
            }
         }
         gdata$U=U1
         names(gdata$U)=names(gdata0$U)
         for(j in 1:xk){
           m=gdata0$X[[j]]
           m=as.matrix(m[-index,])
           #m=m[ok,]
           X1[[j]]=as.matrix(m)
           colnames(X1[[j]])=colnames(gdata0$X[[j]])

         }
         gdata$X=X1
         #colnames(gdata$X[[1]])
         YD=as.matrix(Y[-index,])
         #YD=as.matrix(YD[ok,])
         colnames(YD)=colnames(Y)
         gdata$Y=YD
         res=lmm.gls.eval(gdata)
         res$Var
         for(hh in 1:TraitNum){
             if(hh==1){
                 jvar=cbind(res$Var[[hh]][,1])
                 jfe=cbind(res$FixedEffect[[hh]][,1])
             }
             else{
                 jvar=cbind(jvar,res$Var[[hh]][,1])
                 jfe=cbind(jfe,res$FixedEffect[[hh]][,1])
             }

         }
         jvar=as.matrix(jvar)
         jfe=as.matrix(jfe)
         JACVAR=t(jvar)
         JACB=t(jfe)
         if(mk>0){
             for(hh in 1:TraitNum){
               if(hh==1)jre=cbind(res$RandomEffect[[hh]][,1])
               else jre=cbind(jre,res$RandomEffect[[hh]][,1])

             }
             jre=as.matrix(jre)

            JACPRE=t(jre)
            jackpre[[i]]=JACPRE
         }
         JAC1[[i]]=as.matrix(JACVAR)
         if(mk>0)JAC2[[i]]=as.matrix(JACPRE)
         #if(nx==1)JACB
         JAC3[[i]]=as.matrix(JACB)
         

      }
      if(mk>0)JACKPRE[[k]]=jackpre
      for(s in 1:TraitNum){
          #s=1
          for(j in 1:JacNum){
             if(mk==0&&TraitNum==1)JacVar[j,]=JAC1[[j]]
             else JacVar[j,]=JAC1[[j]][s,]
             if(mk>0)JacPre[j,]=JAC2[[j]][s,]
             if(nx==1&&TraitNum==1)JacB[j,]=JAC3[[j]]
             else JacB[j,]=JAC3[[j]][s,]

          }
          for(j in 1:JacNum){
             index=which(JacVar[j,]<0)
             JacVar[j,index]=0
          }
          ####Calculate the variance components and their proportions from the jackknife estimates
          vt=as.vector(apply(JacVar,1,sum))
          for(i in 1:ml){
             m=mean(JacVar[,i])
             se=sqrt(var(JacVar[,i]))
             if(m<=0){
               m=0
               t1=0
               p=1.0
             }
             else if(m>0){
               t1=m/se
               p=pt(-abs(t1),df)
               p=1-(1-p)^(JacNum-2)
             }
             Vp[s,i]=Vp[s,i]+p
             Vm[s,i]=Vm[s,i]+m
             Cv[s,i]=Cv[s,i]+se
             m=mean(JacVar[,i]/vt)
             se=sqrt(var(JacVar[,i]/vt))
             if(m<=0){
               m=0
               t1=0
               p=1.0
             }
             else if(m>0){
               t1=m/se
               p=pt(-abs(t1),df)
               p=1-(1-p)^(JacNum-2)
             }
             VCp[s,i]=VCp[s,i]+p
             VCm[s,i]=VCm[s,i]+m
             VCv[s,i]=VCv[s,i]+se
             
          }
          #if(is.null(SIMUYES)){ 
          if(mk>0){
             for(i in 1:mt){
               m=mean(JacPre[,i])
               se=sqrt(var(JacPre[,i]))
               t1=m/se
             
               p=pt(-abs(t1),df)
               if(se<=0.0000001)p=1.0
               #p=a$p.value
               p=1-(1-p)^(JacNum-2)
               Pp[s,i]=Pp[s,i]+p
               Pm[s,i]=Pm[s,i]+m
               Cp[s,i]=Cp[s,i]+se
             }
          }
          #}
          for(i in 1:nx){
             m=mean(JacB[,i])
             se=sqrt(var(JacB[,i]))
             t1=m/se
             p=pt(-abs(t1),df)
             if(se<=0.0000001)p=1.0
             #p=a$p.value
             p=1-(1-p)^(JacNum-2)
             Bp[s,i]=Bp[s,i]+p
             Bm[s,i]=Bm[s,i]+m
             Cb[s,i]=Cb[s,i]+se
          }

      }
   }


   ### Variance components
   Vm=Vm/JacRep
   Vp=Vp/JacRep
   Cv=Cv/JacRep
   Cv1=Vm-t*Cv
   Cv2=Vm+t*Cv

   ### Proportional variance components
   VCm=VCm/JacRep
   VCp=VCp/JacRep
   VCv=VCv/JacRep
   VCv1=VCm-t*VCv
   VCv2=VCm+t*VCv
   #if(is.null(SIMUYES)){
   if(mk>0){
     ##Random effects   
     Pm=Pm/JacRep
     Pp=Pp/JacRep
     Cp=Cp/JacRep
     Cp1=Pm-t*Cp
     Cp2=Pm+t*Cp
   }
   #}
   ##Fixed effwaects
   Bm=Bm/JacRep
   Bp=Bp/JacRep
   Cb=Cb/JacRep
   Cb1=Bm-t*Cb
   Cb2=Bm+t*Cb
   
   ##

   VAR=NULL
   PRE=NULL
   VC=NULL
   B=NULL
   CL=c(ALPHA/2,1-ALPHA/2)*100
   CL2=c("%LL","%UL")
   CL=paste(CL,CL2,sep="")
   CNames=c("Estimate","SE","PValue",CL)

   #CNames=c("Estimate","SE","PValue","2.5%CL","97.5%CU")
   for(i in 1: TraitNum){
      VAR[[i]]=data.frame(Vm[i,],Cv[i,],Vp[i,],Cv1[i,],Cv2[i,])
      VC[[i]]=data.frame(VCm[i,],VCv[i,],VCp[i,],VCv1[i,],VCv2[i,])
      if(mk>0)PRE[[i]]=data.frame(Pm[i,],Cp[i,],Pp[i,],Cp1[i,],Cp2[i,])
      B[[i]]=data.frame(Bm[i,],Cb[i,],Bp[i,],Cb1[i,],Cb2[i,])

      rownames(VAR[[i]])=VNames
      colnames(VAR[[i]])=CNames
      rownames(VC[[i]])=paste(VNames,"/VP",sep="")
      colnames(VC[[i]])=CNames
      #if(is.null(SIMUYES)){
        if(mk>0){
           rownames(PRE[[i]])=Str
           Rnames=CNames
           Rnames[1]="Pre"
           colnames(PRE[[i]])=Rnames
        }
      #}
      colnames(B[[i]])=CNames
      rownames(B[[i]])=XNames

   }
   
   names(VAR)=TraitNames
   names(VC)=TraitNames
   if(mk>0)names(PRE)=TraitNames
   names(B)=TraitNames
   ##a=list(Var=VAR,PVar=VC,RandomEffect=PRE,FixedEffect=B,Method=Method,JacNum=JacNum,JacRep=JacRep,VNames=VNames,Y=Y,X=X,U=U)
   #if(is.null(SIMUYES))a=list(Var=VAR,PVar=VC,RandomEffect=PRE,FixedEffect=B,JackPre=JACKPRE)##,Method=Method,JacNum=JacNum,JacRep=JacRep,VNames=VNames,Y=Y,X=X,U=U)
   #else a=list(Var=VAR,PVar=VC,FixedEffect=B)
   if(is.null(LUP))Prediction="AUP"
   else Prediction="LUP"
   #a=list(Var=VAR,PVar=VC,RandomEffect=PRE,FixedEffect=B,JackPre=JACKPRE,Prediction=Prediction)
   #a=list(Var=VAR,PVar=VC,RandomEffect=PRE,FixedEffect=B,Prediction=Prediction)
   a=list(Var=VAR,PVar=VC,RandomEffect=PRE,FixedEffect=B)
   #class(a)="jackknife.result"
   #a$call=match.call()
   #a$Var
   return(a)
}

gls.mq0=function(gdata,au=NULL,R=NULL){
   if(is.null(R))R=NULL
   if(is.null(au))au=NULL
   Y0=as.matrix(gdata$Y)
   #Y0=as.matrix(Y0)
   TraitNum=ncol(Y0)
   TraitNames=colnames(Y0)
   Var=list()
   RE=NULL
   FE=list()
   #i=1
   mk=length(gdata$U)
   for(i in 1:TraitNum){
      #i=1
      gdata$Y=as.matrix(Y0[,i])
      res=gls.mq(gdata,au=au,R=R)
      Var[[i]]=res$Var
      if(mk>0)RE[[i]]=res$RandomEffect
      #else RE[[i]]=NULL
      FE[[i]]=res$FixEffect
   }
   #res$FixEffect
   names(Var)=TraitNames
   if(!is.null(RE))names(RE)=TraitNames
   names(FE)=TraitNames
   res=list(Var=Var,FixedEffect=FE,RandomEffect=RE)
   mk=length(gdata$U)
   Uk=numeric()
   BigU=NULL
   for(u in 1:mk){
       Uk[u]=ncol(gdata$U[[u]])
       BigU=cbind(BigU,gdata$U[[u]])
   }
   mt=sum(Uk)
   xk=length(gdata$X)
   Xk=numeric()
   X=NULL
   for(i in 1:xk){
      Xk[i]=ncol(gdata$X[[i]])
      X=cbind(X,gdata$X[[i]])
   }
   TraitNum=ncol(gdata$Y)
   Residual=NULL
   for(j in 1:TraitNum){
       a=0
       re=res$RandomEffect[[j]][,1]
       fe=res$FixedEffect[[j]][,1]
       if(mk>0)a=a+BigU%*%re
       a=a+X%*%fe
       a=gdata$Y[,j]-a
       Residual=cbind(Residual,a)
   }
   colnames(Residual)=colnames(gdata$Y)
   res$Residual=Residual
   gdata$Y=Y0
   return(res)
}



gls.mq=function(gdata,au=NULL,R=NULL){
   if(is.null(R))R=NULL
   if(is.null(au))au=NULL
   res0=gls.est(gdata,R=R,au=au)
   res=gls.est(gdata,au=res0$Var[,1],R=R)
   ##Statistical tests for variance components
   h=diag(svdinv(res$HV))*2
   Chi_sq=(res0$Var[,1])^2/h
   P_value=(1-pchisq(Chi_sq,1))/2
   Var=data.frame(res0$Var,sqrt(h),Chi_sq,P_value)
   colnames(Var)=c("Est","SE","Chi_sq","P_value")
   res$Var=Var
   
   ##Statistical tests for fixed effects
   res0$FixEffect[,1]
   h=diag(svdinv(res$HF[[1]]))
   h=sqrt(h)
   z=res0$FixEffect[,1]/h
   P_value=1-pchisq(z^2,1)
   attributes(res)
   FixEffect=data.frame(res0$FixEffect,h,z,P_value)
   colnames(FixEffect)=c("Est","SE","z_value","P_value")
   res$FixEffect=FixEffect
   
   ##Statistical tests for random effect  
   mk=length(gdata$U)
   if(mk>0){
     h=get_HR(gdata,res$Qa,res0$Var[,1])
     z=res0$RandomEffect[,1]/h
     res0$RandomEffect[,1][1]/h[1]
     P_value=1-pchisq(z^2,1)
     RandomEffect=data.frame(res0$RandomEffect,h,z,P_value)
     colnames(RandomEffect)=c("Pre","SE","z_value","P_value")
     res$RandomEffect=RandomEffect
   }
   return(res)
}

gls.est=function(gdata,au=NULL,R=NULL){
   mk=length(gdata$U)
   n=nrow(gdata$Y)
   Y=gdata$Y
   U=gdata$U
   C=gdata$C
   xk=length(gdata$X)
   if(xk==1)X=gdata$X[[1]]
   if(xk>1)X=BigX(gdata$X)
   if(is.null(au))au=rep(1,(mk+1))
   if(is.null(R))R=NULL
   if(is.null(gdata$X)){X=matrix(1,nrow=n,ncol=1);colnames(X)=c("mu")}
   if(is.null(C))C=NULL
   res=gls.its(X=X,U=U,Y=Y,au=au,C=C,R=R)
   return(res)
}

gls.its=function(X=NULL,U,Y,au=NULL,C=NULL,R=NULL){
    mk=length(U)
    n=nrow(Y)
    if(is.null(au))au=rep(1,(mk+1))
    if(is.null(R))R=NULL
    if(is.null(X)){X=matrix(1,nrow=n,ncol=1);colnames(X)=c("mu")}
    if(is.null(C))C=NULL
    xnames=colnames(X)
    if(mk==0)vnames="V(e)"
    else if(mk>0){
       rnames=NULL
       vnames=c(paste("V(",names(U),")",sep=""),"V(e)")
       for(i in 1:mk)rnames=c(rnames,colnames(U[[i]]))
    }
    aa=GetQa.gls(au=au,U=U,X=X,C=C,R=R)
    Qa=aa$Qa
    VI=aa$VI
    ML=MINQUE_L(Qa,U)
    if(class(Y)=="numeric")Y=data.frame(Y)
    TraitNum=ncol(Y)
    TraitNames=colnames(Y)
    Var=matrix(0,(mk+1),TraitNum)
    B=matrix(0,ncol(X),TraitNum)
    hf=list()
    if(mk>0)Pre=matrix(0,length(rnames),TraitNum)
    #if(is.null(LUP))VI=NULL
    for(i in 1:TraitNum){
       #i=1
       MR=MINQUE_R(Qa,U,Y[,i])
       a=MINQUEVarPre(Qa,X,U,ML,MR,Y[,i],C,VI)
       Var[,i]=a[[1]]
       B[,i]=a[[2]]
       if(mk>0)Pre[,i]=a[[3]]
       hf[[i]]=a$HF
    }
    
    colnames(Var)=TraitNames
    rownames(Var)=vnames
    colnames(B)=TraitNames
    rownames(B)=xnames
    names(hf)=TraitNames
    if(mk>0){
      colnames(Pre)=TraitNames
      rownames(Pre)=rnames
    }
    if(mk>0)result=list(Var=data.frame(Var),
                 FixEffect=data.frame(B),
                 RandomEffect=data.frame(Pre),
                 HV=ML,
                 Qa=Qa,
                 HF=hf)
    if(mk==0)result=list(Var=data.frame(Var),
                 FixEffect=data.frame(B),
                 RandomEffect=NULL,
                 HV=ML,
                 Qa=Qa,
                 HF=hf)
    result$call=match.call()
    return(result)
}

GetQa.gls=function(au,U,X,C=NULL,R=NULL){
  svdinv=LargeInv
  mk=length(U)
  n=nrow(X)
  if(is.null(R))R=NULL
  if(mk>0)VI=Woodbury(au=au,U=U,R=R)
  else if(mk==0){
      if(is.null(R))VI=diag(1/au,n)
      else VI=svdinv(R)
  }
  if(ncol(X)==1)VX=VI%*%X[,1]
  else VX=VI%*%X
  if(is.null(C))Qa=VI-VX%*%solve(t(X)%*%VX)%*%t(VX)
  else Qa=VI-VX%*%svdinv(t(X)%*%VX+C%*%t(C))%*%t(VX)
  res=list(Qa=Qa,VI=VI)
  return(res)

}

form=function(id,repeated,data=list()){
   name.id=paste(deparse(substitute(id)))
   name.repeated=paste(deparse(substitute(repeated)))
   form=paste("~ ",name.id,"+",name.repeated,sep="")
   form=as.formula(form)
}

un_v=function(mat){
   ok=complete.cases(mat)
   mat=mat[ok,]
   cov=cov(mat)
   d=diag(cov)
   vmat=diag(d)
}
cs_r=function(mat){
   nc=ncol(mat)
   mt=NULL
   for(i in 1:(nc-1)){
      for(j in (i+1):nc){
         if(i==1&&j==2)mt=cbind(mat[,i],mat[,j])
         else mt=rbind(mt,cbind(mat[,i],mat[,j]))
      }
    }
    ok=complete.cases(mt)
    mt=mt[ok,]
    r=cor(mt[,1],mt[,2])
    return(r)
}
ar1_r=function(mat){
   nc=ncol(mat)
   mt=NULL
   for(i in 1:(nc-1)){
      j=i+1
      if(i==1)mt=cbind(mat[,i],mat[,j])
      else mt=rbind(mt,cbind(mat[,i],mat[,j]))
    }
    ok=complete.cases(mt)
    mt=mt[ok,]
    r=cor(mt[,1],mt[,2])
    return(r)
}

getRR=function(Uid,rmat,subid,repid){  ##will be discarded
   R=Uid%*%t(Uid)
   rs=nrow(rmat)
   nrow=nrow(R)
   gs=unique(subid)
   gr=unique(repid)
   imat=list()
   for(i in 1:length(gs)){
       id=which(subid==gs[i])
       if(length(id)==length(gr))v1=c(1:length(gr))
       else{
          v1=numeric(length(id))
          s1=repid[id]
          e=1
          for(j in 1:length(gr)){
             if(length(which(s1==gr[j]))>0){v1[e]=j;e=e+1}
          }
       }
       imat[[i]]=v1
   }
   e1=1
   while(i<=nrow){
      e=0
      for(j in i:(i+rs)){
         if(j<=nrow){
            if(R[i,j]==1){
               if(i!=j)R[i,j]=R[j,i]=rmat[1,j-i+1]
               R[j,j]=rmat[j-i+1,j-i+1]
               e=e+1
            }
         }
      }
      i=i+e
   }
   return(R)
}

getR=function(Ur,rmat,subid,repid){
   R=Ur%*%t(Ur)
   rs=nrow(rmat)
   nrow=nrow(R)
   gs=unique(subid)
   gr=unique(repid)
   imat=list()
   for(i in 1:length(gs)){
       id=which(subid==gs[i])
       if(length(id)==length(gr))v1=c(1:length(gr))
       else{
          v1=numeric(length(id))
          s1=repid[id]
          e=1
          for(j in 1:length(gr)){
             if(length(which(s1==gr[j]))>0){v1[e]=j;e=e+1}
          }
       }
       imat[[i]]=v1
   }
   e=0
   for(i in 1:length(gs)){
      v1=imat[[i]]
      for(j1 in 1:length(v1)){
         p1=v1[j1]
         for(j2 in 1:length(v1)){
             p2=v1[j2]
             if(j1==j2)R[e+j1,e+j1]=rmat[p1,p1]
             else R[e+j1,e+j2]=R[e+j2,e+j1]=rmat[p1,p2]
         }
      }
      e=e+length(v1)
   }
   return(R)
}


getrmmat=function(v,id,group){
   dat=data.frame(v,id,group)
   dat=dat[order(id,group),]
   gn=length(unique(group))
   ni=length(unique(id))
   mat=matrix(0,nrow=ni,ncol=gn)
   idnames=unique(id)
   gpnames=unique(group)
   for(i in 1:ni){
       for(j in 1:gn){
          id=which(dat$id==idnames[i]&dat$group==gpnames[j])
          if(length(id)>0)mat[i,j]=dat$v[id]
          else mat[i,j]=NA
       }
   }
   return(mat)
}

getglsid=function(subject,repeated,data=list()){
   if(is.null(data))data=NULL
   #name.rep=paste(deparse(substitute(repeated)))
   #form=paste("~ ",name.rep,sep="")
   #form=as.formula(form)
   rep.id=model.frame(repeated,data=data)
   rep.id=rep.id[,1]
   #if(is.null(type))type="cs"
   subx=model.frame(subject,data=data)
   sub.id=CFactor(subx,v=c(1:ncol(subx)))
   list(sub.id=sub.id,rep.id=rep.id)
}

glsdata.sort=function(subject,repeated,data=list()){
    aa=getglsid(subject,repeated,data)
    subid=aa$sub.id
    repid=aa$rep.id
    #dat=data.frame(subid,repid)
    ok=order(subid,repid)
    data=data[ok,]
}

getglsindex=function(subid,repid,index){
   #subid=gdata0$sub.id
   #repid=gdata0$rep.id
   id=unique(subid)[-index]
   for(i in 1:length(id)){
      if(i==1)nid=which(subid==id[i])
      else nid=c(nid,which(subid==id[i]))
   }
   return(nid)
}


