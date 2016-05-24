#qgmod.mq=mixed.genmod.minque
#qgmod.mq.jack=mixed.genmod.minque.jack
#qgmod.mq.perm=mixed.genmod.minque.perm

#stab.fw=genmod.fw.reg
#stab.mean=genmod.rank
#stab.var=genmod.var
#stab.ammi=genmod.ammi



ad4.mq=function(Y,Ped){
  res=qgmod.mq(Y,Ped,Model="AD",Cross="Four")
  return(res)
}
adc4.mq=function(Y,Ped){
  res=qgmod.mq(Y,Ped,Model="ADC",Cross="Four")
  return(res)

}

ad4.mq.jack=function(Y,Ped,JacNum=NULL,JacRep=NULL){
  if(is.null(JacNum))JacNum=10
  if(is.null(JacRep))JacRep=1
  res=qgmod.mq.jack(Y,Ped,Model="AD",Cross="Four",JacNum=JacNum,JacRep=JacRep)
  return(res)
}

adc4.mq.jack=function(Y,Ped,JacNum=NULL,JacRep=NULL){
  if(is.null(JacNum))JacNum=10
  if(is.null(JacRep))JacRep=1
  res=qgmod.mq.jack(Y,Ped,Model="ADC",Cross="Four",JacNum=JacNum,JacRep=JacRep)
  return(res)
}

ad4.reml=function(Y,Ped){
  res=qgmod.reml(Y,Ped,Model="AD",Cross="Four")
  return(res)
}
adc4.reml=function(Y,Ped){
  res=qgmod.reml(Y,Ped,Model="ADC",Cross="Four")
  return(res)

}

ad4.reml.jack=function(Y,Ped,JacNum=NULL,JacRep=NULL){
  if(is.null(JacNum))JacNum=10
  if(is.null(JacRep))JacRep=1
  res=qgmod.reml.jack(Y,Ped,Model="AD",Cross="Four",JacNum=JacNum,JacRep=JacRep)
  return(res)
}

adc4.reml.jack=function(Y,Ped,JacNum=NULL,JacRep=NULL){
  if(is.null(JacNum))JacNum=10
  if(is.null(JacRep))JacRep=1
  res=qgmod.reml.jack(Y,Ped,Model="ADC",Cross="Four",JacNum=JacNum,JacRep=JacRep)
  return(res)
}



adrc.mq=function(Y,Ped,Row=NULL,Col=NULL){
  gdata=genmod.data(Y,Ped,Model="AD",Cross="Two")
  en=unique(Ped[,1])
  mat=cbind(Row,Col)
  if(length(en)>1){
     if(!is.null(Row))Row=paste(Ped[,1],Row,sep=":")
     if(!is.null(Col))Col=paste(Ped[,1],Col,sep=":")
  }
  if(!is.null(mat))URC=GenerateU(mat)
  gdata$U=c(gdata$U,URC)
  gdata$VC=c(gdata$VC,1,1)
  res=mq0(gdata)
  #res$FixedEffect=res$FixEffect
  #res$RandomFixedEffect=res$FixEffect
  return(res)
}

adrc.mq.jack=function(Y,Ped,Row=NULL,Col=NULL,JacNum=NULL,JacRep=NULL){
  if(is.null(JacNum))JacNum=10
  if(is.null(JacRep))JacRep=1
  gdata=genmod.data(Y,Ped,Model="AD",Cross="Two")
  en=unique(Ped[,1])
  mat=cbind(Row,Col)
  if(length(en)>1){
     if(!is.null(Row))Row=paste(Ped[,1],Row,sep=":")
     if(!is.null(Col))Col=paste(Ped[,1],Col,sep=":")
  }
  if(!is.null(mat))URC=GenerateU(mat)
  gdata$U=c(gdata$U,URC)
  gdata$VC=c(gdata$VC,1,1)
  res=genmod.jack(gdata,JacNum=JacNum,JacRep=JacRep)
  #res$FixedEffect=res$FixEffect
  #res$RandomFixedEffect=res$FixEffect
  return(res)
}


adrc.reml=function(Y,Ped,Row=NULL,Col=NULL){
  gdata=genmod.data(Y,Ped,Model="AD",Cross="Two")
  en=unique(Ped[,1])
  mat=cbind(Row,Col)
  if(length(en)>1){
     if(!is.null(Row))Row=paste(Ped[,1],Row,sep=":")
     if(!is.null(Col))Col=paste(Ped[,1],Col,sep=":")
  }
  if(!is.null(mat))URC=GenerateU(mat)
  gdata$U=c(gdata$U,URC)
  gdata$VC=c(gdata$VC,1,1)
  res=reml0(gdata)
  #res$FixedEffect=res$FixEffect
  #res$RandomFixedEffect=res$FixEffect
  return(res)
}

adrc.reml.jack=function(Y,Ped,Row=NULL,Col=NULL,JacNum=NULL,JacRep=NULL){
  if(is.null(JacNum))JacNum=10
  if(is.null(JacRep))JacRep=1
  gdata=genmod.data(Y,Ped,Model="AD",Cross="Two")
  en=unique(Ped[,1])
  mat=cbind(Row,Col)
  if(length(en)>1){
     if(!is.null(Row))Row=paste(Ped[,1],Row,sep=":")
     if(!is.null(Col))Col=paste(Ped[,1],Col,sep=":")
  }
  if(!is.null(mat))URC=GenerateU(mat)
  gdata$U=c(gdata$U,URC)
  gdata$VC=c(gdata$VC,1,1)
  res=genmod.reml.jack(gdata,JacNum=JacNum,JacRep=JacRep)
  #res$FixedEffect=res$FixEffect
  #res$RandomFixedEffect=res$FixEffect
  return(res)
}

ad.mq=function(Y,Ped){
  res=qgmod.mq(Y,Ped,Model="AD",Cross="Two")
  return(res)
}

adc.mq=function(Y,Ped){
  res=qgmod.mq(Y,Ped,Model="ADC",Cross="Two")
  return(res)

}
adm.mq=function(Y,Ped){
  res=qgmod.mq(Y,Ped,Model="ADM",Cross="Two")
  return(res)

}

adaa.mq=function(Y,Ped){
  res=qgmod.mq(Y,Ped,Model="ADAA",Cross="Two")
  return(res)

}


ad.mq.jack=function(Y,Ped,JacNum=NULL,JacRep=NULL){
  if(is.null(JacNum))JacNum=10
  if(is.null(JacRep))JacRep=1
  res=qgmod.mq.jack(Y,Ped,Model="AD",Cross="Two",JacNum=JacNum,JacRep=JacRep)
  return(res)
}

adc.mq.jack=function(Y,Ped,JacNum=NULL,JacRep=NULL){
  if(is.null(JacNum))JacNum=10
  if(is.null(JacRep))JacRep=1
  res=qgmod.mq.jack(Y,Ped,Model="ADC",Cross="Two",JacNum=JacNum,JacRep=JacRep)
  return(res)
}

adm.mq.jack=function(Y,Ped,JacNum=NULL,JacRep=NULL){
  if(is.null(JacNum))JacNum=10
  if(is.null(JacRep))JacRep=1
  res=qgmod.mq.jack(Y,Ped,Model="ADM",Cross="Two",JacNum=JacNum,JacRep=JacRep)
  return(res)
}

adaa.mq.jack=function(Y,Ped,JacNum=NULL,JacRep=NULL){
  if(is.null(JacNum))JacNum=10
  if(is.null(JacRep))JacRep=1
  res=qgmod.mq.jack(Y,Ped,Model="ADAA",Cross="Two",JacNum=JacNum,JacRep=JacRep)
  return(res)
}

qgmod.mq=function(Y,Ped,Model,Cross){
   gdata=genmod.data(Y,Ped,Model,Cross)
   res=mq0(gdata)
   #res$FixedEffect=res$FixEffect
   return(res)
}

qgmod.mq.jack=function(Y,Ped,Model,Cross,JacNum=NULL,JacRep=NULL){
   if(is.null(JacNum))JacNum=10
   if(is.null(JacRep))JacRep=1
   gdata=genmod.data(Y,Ped,Model,Cross)
   res=genmod.jack(gdata,JacNum=JacNum,JacRep=JacRep)
   return(res)
}

qgmod.reml=function(Y,Ped,Model,Cross){
   gdata=genmod.data(Y,Ped,Model,Cross)
   res=reml0(gdata)
   #res$FixedEffect=res$FixEffect
   return(res)
}

qgmod.reml.jack=function(Y,Ped,Model,Cross,JacNum=NULL,JacRep=NULL){
   if(is.null(JacNum))JacNum=10
   if(is.null(JacRep))JacRep=1
   gdata=genmod.data(Y,Ped,Model,Cross)
   res=genmod.reml.jack(gdata,JacNum=JacNum,JacRep=JacRep)
   return(res)
}


mixed.genmod.minque.perm=function(Y,Ped,Model,Cross,PermNum=NULL){
   if(is.null(PermNum))PermNum=500
   gdata=genmod.data(Y,Ped,Model,Cross)
   res=genmod.perm(gdata=gdata,au=NULL,PermNum=PermNum)
   return(res)
}

ad.reml=function(Y,Ped){
  res=qgmod.reml(Y,Ped,Model="AD",Cross="Two")
  return(res)
}

adc.reml=function(Y,Ped){
  res=qgmod.reml(Y,Ped,Model="ADC",Cross="Two")
  return(res)

}
adm.reml=function(Y,Ped){
  res=qgmod.reml(Y,Ped,Model="ADM",Cross="Two")
  return(res)

}

adaa.reml=function(Y,Ped){
  res=qgmod.reml(Y,Ped,Model="ADAA",Cross="Two")
  return(res)

}


ad.reml.jack=function(Y,Ped,JacNum=NULL,JacRep=NULL){
  if(is.null(JacNum))JacNum=10
  if(is.null(JacRep))JacRep=1
  res=qgmod.reml.jack(Y,Ped,Model="AD",Cross="Two",JacNum=JacNum,JacRep=JacRep)
  return(res)
}

adc.reml.jack=function(Y,Ped,JacNum=NULL,JacRep=NULL){
  if(is.null(JacNum))JacNum=10
  if(is.null(JacRep))JacRep=1
  res=qgmod.reml.jack(Y,Ped,Model="ADC",Cross="Two",JacNum=JacNum,JacRep=JacRep)
  return(res)
}

adm.reml.jack=function(Y,Ped,JacNum=NULL,JacRep=NULL){
  if(is.null(JacNum))JacNum=10
  if(is.null(JacRep))JacRep=1
  res=qgmod.reml.jack(Y,Ped,Model="ADM",Cross="Two",JacNum=JacNum,JacRep=JacRep)
  return(res)
}

adaa.reml.jack=function(Y,Ped,JacNum=NULL,JacRep=NULL){
  if(is.null(JacNum))JacNum=10
  if(is.null(JacRep))JacRep=1
  res=qgmod.reml.jack(Y,Ped,Model="ADAA",Cross="Two",JacNum=JacNum,JacRep=JacRep)
  return(res)
}

ad.simudata=function(Y,Ped,v,b,SimuNum=NULL){
   if(is.null(SimuNum))SimuNum=200
   gdata=genmod.data(Y,Ped,Model="AD",Cross="Two")
   gdata=SimuData(gdata,v=v,b=b,SimuNum=SimuNum)
   YS=gdata$Y
   return(YS)
}

adc.simudata=function(Y,Ped,v,b,SimuNum=NULL){
   if(is.null(SimuNum))SimuNum=200
   gdata=genmod.data(Y,Ped,Model="ADC",Cross="Two")
   gdata=SimuData(gdata,v=v,b=b,SimuNum=SimuNum)
   YS=gdata$Y
   return(YS)
}

adaa.simudata=function(Y,Ped,v,b,SimuNum=NULL){
   if(is.null(SimuNum))SimuNum=200
   gdata=genmod.data(Y,Ped,Model="ADAA",Cross="Two")
   gdata=SimuData(gdata,v=v,b=b,SimuNum=SimuNum)
   YS=gdata$Y
   return(YS)
}

adm.simudata=function(Y,Ped,v,b,SimuNum=NULL){
   if(is.null(SimuNum))SimuNum=200
   gdata=genmod.data(Y,Ped,Model="ADM",Cross="Two")
   gdata=SimuData(gdata,v=v,b=b,SimuNum=SimuNum)
   YS=gdata$Y
   return(YS)
}

ad.simu=function(Y,Ped,method=NULL,ALPHA=NULL){
   if(is.null(ALPHA))ALPHA=0.05
   if(is.null(method))method=c("minque")
   gdata=genmod.data(Y,Ped,Model="AD",Cross="Two")
   ml=length(gdata$U)+1
   gdata$VC=rep(1,ml)
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


adc.simu=function(Y,Ped,method=NULL,ALPHA=NULL){
   if(is.null(ALPHA))ALPHA=0.05
   if(is.null(method))method=c("minque")
   gdata=genmod.data(Y,Ped,Model="ADC",Cross="Two")
   ml=length(gdata$U)+1
   gdata$VC=rep(1,ml)
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

adaa.simu=function(Y,Ped,method=NULL,ALPHA=NULL){
   if(is.null(ALPHA))ALPHA=0.05
   if(is.null(method))method=c("minque")
   gdata=genmod.data(Y,Ped,Model="ADAA",Cross="Two")
   ml=length(gdata$U)+1
   gdata$VC=rep(1,ml)
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

adm.simu=function(Y,Ped,method=NULL,ALPHA=NULL){
   if(is.null(ALPHA))ALPHA=0.05
   if(is.null(method))method=c("minque")
   gdata=genmod.data(Y,Ped,Model="ADM",Cross="Two")
   ml=length(gdata$U)+1
   gdata$VC=rep(1,ml)
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




ad.simu1 <-
function(Y,Ped,v0,b0,SimuNum=NULL,JacNum=NULL,JacRep=NULL,ALPHA=NULL){
    if(is.null(SimuNum))SimuNum=200
    if(is.null(JacNum))JacNum=10
    if(is.null(JacRep))JacRep=1
    if(is.null(ALPHA))ALPHA=0.05
    gdata=genmod.data(Y,Ped,Model="AD",Cross="Two")
    res=genmod.simu(gdata,v0,b0,SimuNum=SimuNum, JacNum=JacNum, JacRep=JacRep)
    return(res)
}

ad.simu.jack=function(Y,Ped,method=NULL,JacNum=NULL,JacRep=NULL,ALPHA=NULL){
    if(is.null(method))method=c("minque")
    if(is.null(JacNum))JacNum=10
    if(is.null(JacRep))JacRep=1
    if(is.null(ALPHA))ALPHA=0.05
    gdata=genmod.data(Y,Ped,Model="AD",Cross="Two")
    ml=length(gdata$U)+1
    gdata$VC=rep(1,ml)
    #gdata = mixed.data(formula)
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

adc.simu.jack=function(Y,Ped,method=NULL,JacNum=NULL,JacRep=NULL,ALPHA=NULL){
    if(is.null(method))method=c("minque")
    if(is.null(JacNum))JacNum=10
    if(is.null(JacRep))JacRep=1
    if(is.null(ALPHA))ALPHA=0.05
    gdata=genmod.data(Y,Ped,Model="ADC",Cross="Two")
    ml=length(gdata$U)+1
    gdata$VC=rep(1,ml)
    #gdata = mixed.data(formula)
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

adaa.simu.jack=function(Y,Ped,method=NULL,JacNum=NULL,JacRep=NULL,ALPHA=NULL){
    if(is.null(method))method=c("minque")
    if(is.null(JacNum))JacNum=10
    if(is.null(JacRep))JacRep=1
    if(is.null(ALPHA))ALPHA=0.05
    gdata=genmod.data(Y,Ped,Model="ADAA",Cross="Two")
    ml=length(gdata$U)+1
    gdata$VC=rep(1,ml)
    #gdata = mixed.data(formula)
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

adm.simu.jack=function(Y,Ped,method=NULL,JacNum=NULL,JacRep=NULL,ALPHA=NULL){
    if(is.null(method))method=c("minque")
    if(is.null(JacNum))JacNum=10
    if(is.null(JacRep))JacRep=1
    if(is.null(ALPHA))ALPHA=0.05
    gdata=genmod.data(Y,Ped,Model="ADM",Cross="Two")
    ml=length(gdata$U)+1
    gdata$VC=rep(1,ml)
    #gdata = mixed.data(formula)
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



  




adaa.simu1 <-
function(Y,Ped,v0,b0,SimuNum=NULL,JacNum=NULL,JacRep=NULL,ALPHA=NULL){
    if(is.null(SimuNum))SimuNum=200
    if(is.null(JacNum))JacNum=10
    if(is.null(JacRep))JacRep=1
    if(is.null(ALPHA))ALPHA=0.05
    gdata=genmod.data(Y,Ped,Model="ADAA",Cross="Two")
    res=genmod.simu(gdata,v0,b0,SimuNum=SimuNum, JacNum=JacNum, JacRep=JacRep)
    return(res)
}

adc.simu1 <-
function(Y,Ped,v0,b0,SimuNum=NULL,JacNum=NULL,JacRep=NULL,ALPHA=NULL){
    if(is.null(SimuNum))SimuNum=200
    if(is.null(JacNum))JacNum=10
    if(is.null(JacRep))JacRep=1
    if(is.null(ALPHA))ALPHA=0.05
    gdata=genmod.data(Y,Ped,Model="ADC",Cross="Two")
    res=genmod.simu(gdata,v0,b0,SimuNum=SimuNum, JacNum=JacNum, JacRep=JacRep,ALPHA=NULL)
    return(res)
}

adm.simu1 <-
function(Y,Ped,v0,b0,SimuNum=NULL,JacNum=NULL,JacRep=NULL,ALPHA=NULL){
    if(is.null(SimuNum))SimuNum=200
    if(is.null(JacNum))JacNum=10
    if(is.null(JacRep))JacRep=1
    if(is.null(ALPHA))ALPHA=0.05
    gdata=genmod.data(Y,Ped,Model="ADM",Cross="Two")
    res=genmod.simu(gdata,v0,b0,SimuNum=SimuNum, JacNum=JacNum, JacRep=JacRep)
    return(res)
}
