
wide_to_long=function(data,names){
   nc=ncol(data)-1
   nr=nrow(data)
   rnames=as.vector(data[,1])
   cnames=colnames(data)[-1]
   loc=rep(cnames,each=nr)
   gen=rep(rnames,nc)
   y=as.matrix(data[,-1])
   y=as.vector(y)
   dat1=data.frame(gen,loc,y)
   colnames(dat1)=names
   return(dat1)
}
#?rep
genmod.fw.cond=function(y,Gen,Env,times=NULL,Rep=NULL,X,alpha=NULL,...){
   if(is.null(alpha))alpha=0.05
   if(is.null(times))times=1000
   dat=data.frame(Gen,Env,X,y)
   ok=complete.cases(dat)
   dat=dat[ok,]
   Y1=dat[,-c(1,2)]
   CY=GetConVariable(Y1)$ConY
   nc=ncol(CY)-1
   mean(CY[,1])
   res=list()
   if(Rep==TRUE){
      for(i in 1:nc){
         #i=1
         dat1=GetGEMean(CY[,i],dat$Gen,dat$Env,CY[,(nc+1)])
         colnames(dat1)=c("Env","Gen","cy","y")
         res[[i]]=fw.cond(dat1,times,alpha)
      }
   }
   else {
       for(i in 1:nc){
         dat1=data.frame(dat$Env,dat$Gen,CY[,i],CY[,(nc+1)])
         colnames(dat1)=c("Env","Gen","cy","y")
         res[[i]]=fw.cond(dat1,times,alpha)
       }
   }
   names(res)=colnames(CY)[1:nc]
   return(res)
}

fw.cond=function(dat1,times,alpha=NULL,...){
   if(is.null(alpha))alpha=0.05
   p1=alpha/2
   p2=1-p1
   names(dat1)
   EI=tapply(dat1$y,dat1$Env,mean)
   name0=names(EI)
   #EI=as.matrix(EI)
   #id=complete.cases(EI)
   #EI=EI[ok]
   #EI=as.matrix(EI[id,])
   #dim(EI)
   #EI.data<-as.data.frame.table(EI)
   #colnames(EI.data)=c("Env","EI")
   #EI.data<-as.data.frame.table(EI)
   #colnames(EI.data)=c("Env","EI")
   #ok=complete.cases(EI.mdat1)
   #mdat1=mdat1[ok,]
   gnames=sort(unique(dat1$Gen))
   gn=length(unique(dat1$Gen))
   r.square=numeric(gn)
   b=numeric(gn)
   #name0=names(EI)
   DAT=list()
   for(i in 1:gn){
      #i=1
      cat("Interation=",i,"\n")
      id=which(dat1$Gen==gnames[i])
      len=length(EI)
      #if(is.null(X)==F){
         if(length(id)==len)dat2=data.frame(dat1[id,],EI)
         else if(length(id)<len(EI)){
             ei=dat1$y[id]
             name1=dat1$Env[id]
             #name1=names(ei)
             intername=intersect(name1,name0)
             id1=numeric(length(intername))
             id2=numeric(length(intername))
             for(j in 1:length(id1)){
                id0=which(name1==intername[j])
                id1[j]=id0
                id0=which(name0==intername[j])
                id2[j]=id0
             }
             ei1=ei[id1]
             ei2=EI[id2]
             dat2=data.frame(ei1,ei2)
         } 
      #}
      colnames(dat2)=c(colnames(dat1),"EI")
      DAT[[i]]=data.frame(dat2)
   }
   RES=list()
   for(i in 1:gn){
      dat2=DAT[[i]]
      y1=dat2[,3]
      X1=dat2[,5]
      RES[[i]]=reg.boot(y1,X1,times=times)
   }
   names(RES)=gnames
   return(RES)
}
genmod.fw.check=function(y,Gen,Env,times,check,Rep,X=NULL,alpha=NULL,...){
   if(is.null(alpha))alpha=0.05
   if(is.null(X))X=NULL
   if(Rep==TRUE){
       dat=GetGEMean(y,Gen,Env,X)
       y1=dat$y
       if(is.null(X))X1=NULL
       else X1=dat[,-c(1:3)]
       Gen1=dat$Gen
       Env1=dat$Env
       result=fw.check(y1,Gen1,Env1,check,times,X1,alpha)
   }
   else result=fw.check(y,Gen,Env,times,check,X,alpha)


   return(result)
}
#?subset
fw.check=function(y,Gen,Env,times,check,X=NULL,alpha=NULL,...){
   if(is.null(alpha))alpha=0.05
   p1=alpha/2
   p2=1-p1
   y=as.vector(y)
   Env=as.vector(Env)
   Gen=as.vector(Gen)
   if(is.null(X)){
        dat1=data.frame(Gen,Env,y)
        colnames(dat1)=c("Gen","Env","y")
   }
   else {
      dat1=data.frame(Gen,Env,y,X)
      colnames(dat1)=c("Gen","Env","y",colnames(X))
   }
   dat1=dat1[order(dat1$Gen,dat1$Env),]
  
   id=get.check.data(dat1,dat1$Gen,check)
   EI=tapply(dat1$y[id],dat1$Env[id],mean) ##get check based EI
   name0=names(EI)
   Gen1=as.vector(dat1$Gen[-id])
   
   gnames=unique(Gen1)
   gn=length(unique(Gen1))
   r.square=numeric(gn)
   b=numeric(gn)
   #name0=names(EI)
   DAT=list()
   for(i in 1:gn){
      #i=1
      cat("Interation=",i,"\n")
      id=which(dat1$Gen==gnames[i])
      len=length(EI)
      #if(is.null(X)==F){
         if(length(id)==len)dat2=data.frame(dat1[id,],EI)
         else if(length(id)<length(EI)){
             ei=dat1$y[id]
             name1=dat1$Env[id]
             #name1=names(ei)
             intername=intersect(name1,name0)
             id1=numeric(length(intername))
             id2=numeric(length(intername))
             for(j in 1:length(id1)){
                id0=which(name1==intername[j])
                id1[j]=id0
                id0=which(name0==intername[j])
                id2[j]=id0
             }
             ei1=ei[id1]
             ei2=EI[id2]
             dat2=data.frame(ei1,ei2)
         } 
      #}
      colnames(dat2)=c(colnames(dat1),"EI")
      DAT[[i]]=data.frame(dat2)
   }
   
   
   #DAT[[1]]
   #ww03[2,]   
   RES=list()
   for(i in 1:gn){
      #i=1
      dat2=DAT[[i]]
      y1=dat2[,3]
      X1=as.matrix(dat2[,-c(1:3)])
      colnames(X1)=colnames(dat2)[-c(1:3)]
      RES[[i]]=reg.boot(y1,X1,times,boot="residual")
   }
   names(RES)=gnames
   return(RES)
}


genmod.fw.reg=function(y,Gen,Env,times,Rep,X=NULL,alpha=NULL,...){
   if(is.null(alpha))alpha=0.05
   if(is.null(X))X=NULL
   if(Rep==TRUE){
       dat=GetGEMean(y,Gen,Env,X)
       y1=dat$y
       if(is.null(X))X1=NULL
       else X1=dat[,-c(1:3)]
       Gen1=dat$Gen
       Env1=dat$Env
       result=fw_reg(y1,Gen1,Env1,times,X1,alpha)
   }
   else result=fw_reg(y,Gen,Env,times,X,alpha)


   return(result)
}

fw=function(y,Gen,Env,times,X=NULL,alpha=NULL,...){
   if(is.null(alpha))alpha=0.05
   p1=alpha/2
   p2=1-p1
   if(is.null(X)){
        dat1=data.frame(Gen,Env,y)
        colnames(dat1)=c("Gen","Env","y")
   }
   else {
      dat1=data.frame(Gen,Env,y,X)
      colnames(dat1)=c("Gen","Env","y",colnames(X))
   }
   dat1=dat1[order(dat1$Gen,dat1$Env),]
   #return(dat1)
   #
   EI=tapply(dat1$y,dat1$Env,mean)
   name0=names(EI)
   gnames=unique(Gen)
   gn=length(unique(Gen))
   r.square=numeric(gn)
   b=numeric(gn)
   
   DAT=list()
   for(i in 1:gn){
      #i=1
      cat("Interation=",i,"\n")
      id=which(dat1$Gen==gnames[i])
      len=length(EI)
      #if(is.null(X)==F){
         if(length(id)==len)dat2=data.frame(dat1[id,],EI)
         else if(length(id)<length(EI)){
             ei=dat1$y[id]
             name1=dat1$Env[id]
             #name1=names(ei)
             intername=intersect(name1,name0)
             id1=numeric(length(intername))
             id2=numeric(length(intername))
             for(j in 1:length(id1)){
                id0=which(name1==intername[j])
                id1[j]=id0
                id0=which(name0==intername[j])
                id2[j]=id0
             }
             ei1=ei[id1]
             ei2=EI[id2]
             dat2=data.frame(ei1,ei2)
         } 
      #}
      colnames(dat2)=c(colnames(dat1),"EI")
      DAT[[i]]=data.frame(dat2)
   }
   #names(DAT)=gnames
   #return(DAT)
#}  
   RES=list()
   for(i in 1:gn){
      dat2=DAT[[i]]
      y1=dat2[,3]
      X1=dat2[,-c(1:3)]
      RES[[i]]=reg.boot(y1,X1,times,boot="residual")
   }
   names(RES)=gnames
   return(RES)
}

fw_reg=function(y,Gen,Env,times,X=NULL,alpha=NULL,...){
   if(is.null(alpha))alpha=0.05
   p1=alpha/2
   p2=1-p1
   if(is.null(X)){
        dat1=data.frame(Gen,Env,y)
        colnames(dat1)=c("Gen","Env","y")
   }
   else {
      dat1=data.frame(Gen,Env,y,X)
      colnames(dat1)=c("Gen","Env","y",colnames(X))
   }
   dat1=dat1[order(dat1$Gen,dat1$Env),]
   #return(dat1)
   #
   EI=tapply(dat1$y,dat1$Env,mean)
   name0=names(EI)
   gnames=unique(Gen)
   gn=length(unique(Gen))
   r.square=numeric(gn)
   b=numeric(gn)
   
   DAT=list()
   for(i in 1:gn){
      #i=1
      cat("Interation=",i,"\n")
      id=which(dat1$Gen==gnames[i])
      len=length(EI)
      #if(is.null(X)==F){
         if(length(id)==len)dat2=data.frame(dat1[id,],EI)
         else if(length(id)<length(EI)){
             ei=dat1$y[id]
             name1=dat1$Env[id]
             #name1=names(ei)
             intername=intersect(name1,name0)
             id1=numeric(length(intername))
             id2=numeric(length(intername))
             for(j in 1:length(id1)){
                id0=which(name1==intername[j])
                id1[j]=id0
                id0=which(name0==intername[j])
                id2[j]=id0
             }
             ei1=ei[id1]
             ei2=EI[id2]
             dat2=data.frame(ei1,ei2)
         } 
      #}
      colnames(dat2)=c(colnames(dat1),"EI")
      DAT[[i]]=data.frame(dat2)
   }
     
   RES=list()
   if(is.null(X)){
      for(i in 1:gn){
        dat2=DAT[[i]]
        y1=dat2[,3]
        X1=dat2[,-c(1:3)]
        RES[[i]]=reg.boot(y1,X1,times,boot="residual")
      }
   }
   else{
      for(i in 1:gn){
       
        dat2=DAT[[i]]
        nc=ncol(dat2)
        #dat2=dat2[,-nc]
        #y1=dat2[,3]
        #X1=dat2[,-c(1:3)]
        y1=dat2[,nc]
        X1=dat2[,c(4:(nc-1))]

        RES[[i]]=reg.boot(y1,X1,times,boot="residual")
      }
   }
      
   names(RES)=gnames
   return(RES)
}



GetGEMean=function(y,Gen,Env,X=NULL,...){
   GE=tapply(y,list(Env,Gen),mean)
   GE=as.data.frame.table(GE)
   colnames(GE)=c("Env","Gen","y")
   if(is.null(X)==F){
       X=as.matrix(X)
       nx=ncol(X)
       for(i in 1:nx){
           GE1=tapply(X[,i],list(Env,Gen),mean)
           GE1=as.data.frame.table(GE1)
           GE=cbind(GE,GE1[,3])
       }
       colnames(GE)=c("Env","Gen","y",colnames(X))
   }
   dat=data.frame(GE)
   ok=complete.cases(dat)
   dat=dat[ok,]
   return(dat)
}
#class=Geno
genmod.rank=function(Y,class,cls2=NULL,resample,times=NULL,alpha=NULL,...){
   Y=as.matrix(Y)
   tn=ncol(Y)
   if(is.null(alpha))alpha=0.05
   if(is.null(times))times=1000
   if(is.null(cls2))cls2=NULL
   if(tn==1)return(Rank2(Y[,1],class,cls2,resample,times,alpha))
   else{
      RES=list()
      for(i in 1:tn)RES[[i]]=Rank2(Y[,i],class,cls2,resample,times,alpha)
      names(RES)=colnames(Y)
      return(RES)
   }
}

genmod.var=function(Y,class,cls2=NULL,resample,times=NULL,alpha=NULL,...){
   Y=as.matrix(Y)
   tn=ncol(Y)
   if(is.null(alpha))alpha=0.05
   if(is.null(times))times=1000
   if(is.null(cls2))cls2=NULL
   if(tn==1)return(Var2(Y[,1],class,cls2,resample,times,alpha))
   else{
      RES=list()
      for(i in 1:tn)RES[[i]]=Var2(Y[,i],class,cls2,resample,times,alpha)
      names(RES)=colnames(Y)
      return(RES)
   }
}

Rank2=function(y,class,cls2=NULL,resample,times=NULL,alpha=NULL,...){
   if(is.null(alpha))alpha=0.05
   if(is.null(times))times=1000
   p1=alpha/2
   p2=1-p1
   m0=tapply(y,class,mean,na.rm=T)
   nc=length(m0)
   Mean=matrix(0,times,nc)
   n=length(y)
   if(is.null(cls2))cls2=rep(1,n)
   cnames=unique(cls2)
   cn=length(cnames)
   for(i in 1:times){
      #i=1
      if(resample=="Boot"){
         y1=NULL
         class1=NULL
         for(j in 1:cn){
            id=which(cls2==cnames[j])
            y0=y[id]
            c0=class[id]
            n0=length(y0)
            id0=sample(n0,replace=T)
            y1=c(y1,y0[id0])
            class1=c(class1,c0[id0])
         }
         # 
         # 
         #id=sample(n,replace=T)
         #y1=y[id]
         #class1=class[id]
      }
      else{
        id=sample(n,replace=F)
        y1=y[id]
        class1=class
      }
      m1=tapply(y1,class1,mean,na.rm=T)
      if(length(m1)<nc)m1=vectcomp(m0,m1)
      Mean[i,]=m1  ##tapply(y1,class1,mean,na.rm=T)
   }
   classnames=names(m0)
   #cn=length(m)
   CI=matrix(0,nc,2)
   m1=numeric(nc)
   #dim(Mean)
   for(i in 1:nc){
      #i=1
      v=Mean[,i]
      CI[i,]=quantile(v,p=c(p1,p2),na.rm=T)
      m1[i]=mean(v)
   }
   mean=data.frame(m0,m1,CI)
   colnames(mean)=c("Orig",resample,"LL","UL")
   rownames(mean)=classnames

   order0=order(m0)
   rank0=numeric(nc)
   rank0[order0]=1:nc
   Rank1=matrix(0,times,nc)
   v=numeric(nc)
   for(i in 1:times){
      id=order(Mean[i,])
      v[id]=1:nc
      Rank1[i,]=v
   }
   CI=matrix(0,nc,2)
   order1=numeric(nc)
   for(i in 1:nc){
      v=Rank1[,i]
      CI[i,]=quantile(v,p=c(p1,p2),na.rm=T)
      order1[i]=mean(v)
   }

   rnk=data.frame(rank0,order1,CI)
   colnames(rnk)=c("Orig",resample,"LL","UL")
   rownames(rnk)=classnames
   result=list(mean=mean,rank=rnk,alpha=alpha)
   return(result)
}

Var2=function(y,class,cls2=NULL,resample,times=NULL,alpha=NULL,...){
   if(is.null(alpha))alpha=0.05
   if(is.null(times))tB=1000
   p1=alpha/2
   p2=1-p1
   m0=tapply(y,class,var,na.rm=T)
   nc=length(m0)
   Mean=matrix(0,times,nc)
   n=length(y)
   if(is.null(cls2))cls2=rep(1,n)
   cnames=unique(cls2)
   cn=length(cnames)
   for(i in 1:times){
      #i=1
      if(resample=="Boot"){
         y1=NULL
         class1=NULL
         for(j in 1:cn){
            id=which(cls2==cnames[j])
            y0=y[id]
            c0=class[id]
            n0=length(y0)
            id0=sample(n0,replace=T)
            y1=c(y1,y0[id0])
            class1=c(class1,c0[id0])
         }
         # 
         # 
         #id=sample(n,replace=T)
         #y1=y[id]
         #class1=class[id]
      }
      else{
        id=sample(n,replace=F)
        y1=y[id]
        class1=class
      }
      m1=tapply(y1,class1,var,na.rm=T)
      if(length(m1)<nc)m1=vectcomp(m0,m1)
      Mean[i,]=m1  ##tapply(y1,class1,mean,na.rm=T)
   }
   classnames=names(m0)
   #cn=length(m)
   CI=matrix(0,nc,2)
   m1=numeric(nc)
   #dim(Mean)
   for(i in 1:nc){
      #i=1
      v=Mean[,i]
      CI[i,]=quantile(v,p=c(p1,p2))
      m1[i]=mean(v)
   }
   mean=data.frame(m0,m1,CI)
   colnames(mean)=c("Orig",resample,"LL","UL")
   rownames(mean)=classnames
   result=list(Var=mean,alpha=alpha)
   return(result)
}



genmod.var1=function(y,class,resample,times,alpha=NULL,...){
   if(is.null(alpha))alpha=0.05
   a1=variation(y,class,resample,times,alpha)
   a2=cvar(y,class,resample,times,alpha)
   ok=complete.cases(a1)
   a1=a1[ok,]
   ok=complete.cases(a2)
   a2=a2[ok,]
   result=list(variation=a1,coefvar=a2)
   return(result)

}

cv=function (x, square = FALSE) 
{
   n <- length(x)
   V <- sqrt((n - 1) * var(x,na.rm=T)/n)/mean(x,na.rm=T)
   if (square) 
   V <- V^2
   V
}

variation=function(y,class,resample,times,alpha=NULL,...){
   if(is.null(alpha))alpha=0.05
   p1=alpha/2
   p2=1-p1
   v0=tapply(y,class,var,na.rm=T)
   nc=length(v0)
   Var=matrix(0,times,nc)
   n=length(y)
   for(i in 1:times){
      #i=1
      if(resample=="Boot"){
         id=sample(n,replace=T)
         y1=y[id]
         class1=class[id]
      }
      else{
        id=sample(n,replace=F)
        y1=y[id]
        class1=class
      }
      v1=tapply(y1,class1,var,na.rm=T)
      if(length(v1)<nc)v1=vectcomp(v0,v1)
      Var[i,]=v1  ##tapply(y1,class1,mean,na.rm=T)
   }
   classnames=names(v0)
   #cn=length(m)
   CI=matrix(0,nc,2)
   v1=numeric(nc)
   #dim(Mean)
   for(i in 1:nc){
      #i=1
      v=Var[,i]
      CI[i,]=quantile(v,p=c(p1,p2),na.rm=T)
      v1[i]=mean(v)
   }
   var1=data.frame(v0,v1,CI)
   colnames(var1)=c("Orig",resample,"LL","UL")
   rownames(var1)=classnames
   return(var1)
}

cvar=function(y,class,resample,times,alpha=NULL,...){
   if(is.null(alpha))alpha=0.05
   p1=alpha/2
   p2=1-p1
   v0=tapply(y,class,cv)
   nc=length(v0)
   Var=matrix(0,times,nc)
   n=length(y)
   for(i in 1:times){
      #i=1
      if(resample=="Boot"){
         id=sample(n,replace=T)
         y1=y[id]
         class1=class[id]
      }
      else{
        id=sample(n,replace=F)
        y1=y[id]
        class1=class
      }
      v1=tapply(y1,class1,cv)
      if(length(v1)<nc)v1=vectcomp(v0,v1)
      Var[i,]=v1  ##tapply(y1,class1,mean,na.rm=T)
   }
   classnames=names(v0)
   #cn=length(m)
   CI=matrix(0,nc,2)
   v1=numeric(nc)
   #dim(Mean)
   for(i in 1:nc){
      #i=1
      v=Var[,i]
      CI[i,]=quantile(v,p=c(p1,p2),na.rm=T)
      v1[i]=mean(v)
   }
   var1=data.frame(v0,v1,CI)
   colnames(var1)=c("Orig",resample,"LL","UL")
   rownames(var1)=classnames
   return(var1)
}


#B=100
FWStabilityR=function(y,Gen,Env,times=NULL,alpha=NULL,...){
   if(is.null(alpha))alpha=0.05
   if(is.null(times))times=1000
   p1=alpha/2
   p2=1-p1
   dat1=data.frame(y,Gen,Env)
   dat1=dat1[order(dat1$Gen,dat1$Env),]
   EI=tapply(dat1$y,dat1$Env,mean)
   EI.data<-as.data.frame.table(EI)
   colnames(EI.data)=c("Env","EI")
   EI.data<-as.data.frame.table(EI)
   colnames(EI.data)=c("Env","EI")
   gnames=sort(unique(Gen))
   gn=length(unique(Gen))
   r.square=numeric(gn)
   b=numeric(gn)
   name0=names(EI)
   DAT=list()
   for(i in 1:gn){
      #i=28
      cat("Interation=",i,"\n")
      id=which(dat1$Gen==gnames[i])
      d2=dat1[id,]
      ei=tapply(d2$y,d2$Env,mean)
      name1=names(ei)
      intername=intersect(name1,name0)
      id1=numeric(length(intername))
      id2=numeric(length(intername))
      for(j in 1:length(id1)){
         id0=which(name1==intername[j])
         id1[j]=id0
         id0=which(name0==intername[j])
         id2[j]=id0
      }
      ei1=ei[id1]
      ei2=EI[id2]
      dat2=data.frame(ei1,ei2)
      colnames(dat2)=c("y","EI")
      DAT[[i]]=dat2
      reg<-lm(y~EI,data=dat2)
      r.square[i]<-summary(reg)$r.squared
      b[i]<-coef(reg)[2]
   }
   names(r.square)=gnames
   names(b)=gnames
   result=list(r.square=r.square,slope=b)

   b0=b
   r0=r.square
   ##By bootstrapping test ######
    R2=matrix(0,times,gn)
    B1=matrix(0,times,gn)
    for(i in 1:gn){
       dat0=DAT[[i]]
       rc=nrow(dat0)
       #id=which(dat1$Gen==gnames[i])
       #dat0=data.frame(dat1[id,],EI)
       for(k in 1:times){
          id1=sample(rc,replace=T)
          y1=dat0$y[id1]
          x1=dat0$EI[id1]
          dat2=data.frame(y1,x1)
          reg<-lm(y1~x1,data=dat2)
          R2[k,i]=summary(reg)$r.squared
          B1[k,i]<-coef(reg)[2]
      }
    }
    
    CI1=matrix(0,gn,2)
    m1=numeric(gn)
    for(i in 1:gn){
       v=R2[,i]
       CI1[i,]=quantile(v,p=c(p1,p2))
       m1[i]=mean(v)
    }
    CI2=matrix(0,gn,2)
    m2=numeric(gn)
    for(i in 1:gn){
       v=B1[,i]
       CI2[i,]=quantile(v,p=c(p1,p2))
       m2[i]=mean(v)
    }
    rboot=data.frame(r0,m1,CI1)
    colnames(rboot)=c("Orig","Boot","LL","UL")
    rownames(rboot)=gnames
    bboot=data.frame(b0,m2,CI2)
    colnames(bboot)=c("Orig","Boot","LL","UL")
    rownames(bboot)=gnames
    result=list(rsquare=rboot,slope=bboot)
    return(result)
}



#Y=Yld
#ENV=L
#GEN=G
#Y=PWT
#REP=Rep
genmod.rank1=function(y,class,resample,B,alpha=NULL,...){
   if(is.null(alpha))alpha=0.05
   class=as.vector(class)
   result=Rank1(y,class,resample,B,alpha)
   return(result)
}


#y=y1
#X=X1
#data.frame(y1,X1)
reg.boot=function(y,X,times=NULL,boot=NULL,ALPHA=NULL,...){
   if(is.null(boot))boot="original"
   if(is.null(ALPHA))ALPHA=0.05
   if(is.null(times))times=1000
   p1=ALPHA/2
   p2=1-p1

   dat=data.frame(y,X)
   ok=complete.cases(dat)
   dat=dat[ok,]
   n=length(dat$y)
   mod0=lm(y~.,data=dat)
   #attributes(mod0)
   r0=summary(mod0)$r.squared
   res0=mod0$residual
   yhat0=mod0$fitted.values
   b0=mod0$coefficients
   nc=length(b0)
   B=matrix(0,times,nc)
   r1=numeric(times)
   if(boot=="original"){
      bootmethod="boot_original"
      for(i in 1:times){
         id=sample(1:n,replace=T)
         dat1=dat[id,]
         mod1=lm(y~.,data=dat1)
         B[i,]=mod1$coefficients
         r1[i]=summary(mod1)$r.squared
      }
   }

   else if(boot=="residual"){
      bootmethod="boot_residual"
      X1=dat[,-1]
      for(i in 1:times){
         id=sample(1:n,replace=T)
         y1=yhat0+res0[id]
         dat1=data.frame(y1,X1)
         #dat1=dat[id,]
         mod1=lm(y1~.,data=dat1)
         B[i,]=mod1$coefficients
         r1[i]=summary(mod1)$r.squared
      }
   }
   B=data.frame(B,r1)

   p=numeric(nc)
   classnames=c(names(b0),"r.squared")
   CI=matrix(0,nc+1,2)
   b1=numeric(nc+1)
   b0=c(b0,r0)
   for(i in 1:(nc+1)){
      #i=1
      v=B[,i]
      ok=complete.cases(v)
      v=v[ok]
      tm=length(v)
      CI[i,]=quantile(v,p=c(p1,p2),na.rm=T)
      b1[i]=mean(v)
      if(b0[i]==0)p[i]=1.000
      else if(b0[i]>0)p[i]=length(which(v<0))/tm
      else if(b0[i]<0)p[i]=length(which(v>0))/tm
   }
   p=p*2
   mean=data.frame(b0,b1,p,CI)
   colnames(mean)=c("Orig","boot","Prob","LL","UL")
   rownames(mean)=classnames
   res=list(parameters=mean,alpha=ALPHA,bootmethod=bootmethod)
   return(res)
}

reg.perm=function(y,X,times=NULL,perm=NULL,ALPHA=NULL,...){
   if(is.null(perm))perm="original"
   if(is.null(ALPHA))ALPHA=0.05
   if(is.null(times))times=1000
   p1=ALPHA/2
   p2=1-p1

   dat=data.frame(y,X)
   ok=complete.cases(dat)
   dat=dat[ok,]
   n=length(dat$y)
   mod0=lm(y~.,data=dat)
   res0=mod0$residual
   yhat0=mod0$fitted.values
   b0=mod0$coefficients
   r0=summary(mod0)$r.squared
   nc=length(b0)
   B=matrix(0,times,nc)
   r1=numeric(times)
   X0=dat[,-1]
   y0=dat[,1]
   p=numeric(nc)
   classnames=c(names(b0),"r.squared")
   if(perm=="original"){
      permmethod="perm_original"
      for(i in 1:times){
         id=sample(1:n,replace=F)
         y1=y0[id]
         dat1=data.frame(y1,X0)
         mod1=lm(y1~.,data=dat1)
         B[i,]=mod1$coefficients
         r1[i]=summary(mod1)$r.squared
      }
      B=data.frame(B,r1)
      p=numeric(nc)
      #classnames=c(names(b0),"r.squared")
      #CI=matrix(0,nc+1,2)
      #b1=numeric(nc+1)
      b0=c(b0,r0)
      for(i in 1:(nc+1)){
        v=B[,i]
        if(i==1){
          if(b0[i]==0)p[i]=1.000
          else if(b0[i]<0)p[i]=length(which(v<0))/times
          else if(b0[i]>0)p[i]=length(which(v>0))/times
        }

        else if(i>1){
          if(b0[i]==0)p[i]=1.000
          else if(b0[i]>0)p[i]=length(which(v>b0[i]))/times
          else if(b0[i]<0)p[i]=length(which(v<b0[i]))/times
        }
      }
   }
   else if(perm=="residual"){
      permmethod="perm_residual"
      X1=dat[,-1]
      for(i in 1:times){
         id=sample(1:n,replace=F)
         y1=yhat0+res0[id]
         dat1=data.frame(y1,X1)
         #dat1=dat[id,]
         mod1=lm(y1~.,data=dat1)
         B[i,]=mod1$coefficients
         r1[i]=summary(mod1)$r.squared
      }
      B=data.frame(B,r1)
      p=numeric(nc)
      #classnames=c(names(b0),"r.squared")
      #CI=matrix(0,nc+1,2)
      #b1=numeric(nc+1)
      b0=c(b0,r0)

      for(i in 1:(nc+1)){
         v=B[,i]
         if(b0[i]==0)p[i]=1.000
         else if(b0[i]>0)p[i]=length(which(v<0))/times
         else if(b0[i]<0)p[i]=length(which(v>0))/times
      }
   }
   p=p*2
   #classnames=names(b0)
   CI=matrix(0,(nc+1),2)
   b1=numeric(nc+1)
   for(i in 1:(nc+1)){
      #i=1
      v=B[,i]
      CI[i,]=quantile(v,p=c(p1,p2))
      b1[i]=mean(v)
   }
   mean=data.frame(b0,b1,p,CI)
   colnames(mean)=c("Orig","perm","Prob","LL","UL")
   rownames(mean)=classnames
   res=list(parameters=mean,alpha=ALPHA,permmethod=permmethod)
   return(res)
}

Rank1=function(y,class,resample,B=NULL,alpha=NULL,...){
   if(is.null(alpha))alpha=0.05
   if(is.null(B))B=1000
   p1=alpha/2
   p2=1-p1
   m0=tapply(y,class,mean,na.rm=T)
   nc=length(m0)
   Mean=matrix(0,B,nc)
   n=length(y)
   for(i in 1:B){
      #i=1
      if(resample=="Boot"){
         id=sample(n,replace=T)
         y1=y[id]
         class1=class[id]
      }
      else{
        id=sample(n,replace=F)
        y1=y[id]
        class1=class
      }
      m1=tapply(y1,class1,mean,na.rm=T)
      if(length(m1)<nc)m1=vectcomp(m0,m1)
      Mean[i,]=m1  ##tapply(y1,class1,mean,na.rm=T)
   }
   classnames=names(m0)
   #cn=length(m)
   CI=matrix(0,nc,2)
   m1=numeric(nc)
   #dim(Mean)
   for(i in 1:nc){
      #i=1
      v=Mean[,i]
      CI[i,]=quantile(v,p=c(p1,p2))
      m1[i]=mean(v)
   }
   mean=data.frame(m0,m1,CI)
   colnames(mean)=c("Orig",resample,"LL","UL")
   rownames(mean)=classnames

   order0=order(m0)
   rank0=numeric(nc)
   rank0[order0]=1:nc
   Rank1=matrix(0,B,nc)
   v=numeric(nc)
   for(i in 1:B){
      id=order(Mean[i,])
      v[id]=1:nc
      Rank1[i,]=v
   }
   CI=matrix(0,nc,2)
   order1=numeric(nc)
   for(i in 1:nc){
      v=Rank1[,i]
      CI[i,]=quantile(v,p=c(p1,p2))
      order1[i]=mean(v)
   }

   rnk=data.frame(rank0,order1,CI)
   colnames(rnk)=c("Orig",resample,"LL","UL")
   rownames(rnk)=classnames
   result=list(mean=mean,rank=rnk,alpha=alpha)
   return(result)
}
genmod.rank2=function(y,class,cls2,resample,B,alpha=NULL,...){
   if(is.null(alpha))alpha=0.05
   if(is.null(B))B=1000
   class=as.vector(class)
   result=Rank2(y,class,cls2,resample,B,alpha)
   return(result)
}


boot=function(fun,datatimes=NULL,...){
}

get.check.data=function(data,Genotype,GenComNames,...){
   gcn=length(GenComNames)
   for(i in 1:gcn){
      if(i==1)id=which(Genotype==GenComNames[i])
      else id=c(id,which(Genotype==GenComNames[i]))
   }
   return(id)
}

