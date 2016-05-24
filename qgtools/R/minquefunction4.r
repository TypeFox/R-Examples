
##Functions used for linear mixed model analysis

lmm.check=function (formula, data = list()) 
{
    if (is.null(data))gdata = mixed.data(formula)
    else gdata = mixed.data(formula, data)
    mk=length(gdata$U)
    xk=length(gdata$X)
    nx=numeric(xk)
    for(i in 1:xk){
       if(i==1)FixedNames=colnames(gdata$X[[i]])
       else FixedNames=c(FixedNames,colnames(gdata$X[[i]]))
       nx[i]=ncol(gdata$X[[i]])
    }
    names(nx)=names(gdata$X)
    res=list(VarCompNum=(mk+1),VarCompNames=names(gdata$U),FixedCompNum=nx,FixedEffectNames=FixedNames)
    return(res)
}
lmm.mq.simu=function(formula, data = list(),v0,b0,SimuNum=NULL,JacNum=NULL,JacRep=NULL,ALPHA=NULL){
    if (is.null(data))gdata = mixed.data(formula)
    else gdata = mixed.data(formula, data)
    if(is.null(SimuNum))SimuNum=200
    if(is.null(JacNum))JacNum=10
    if(is.null(JacRep))JacRep=1
    if(is.null(ALPHA))ALPHA=0.05
    result=genmod.simuold(gdata,v=v0,b=b0,SimuNum=SimuNum,JacNum=JacNum,JacRep=JacRep,ALPHA=ALPHA)
    return(result)
}    

##Minque(1) approach used for a general mixed linear model
mixed.minque=function(formula,data=list()){
   if(is.null(data))gdata=mixed.data(formula)
   else gdata=mixed.data(formula,data)
   res=genmod(gdata)
   res$FixedEffect=res$FixEffect
   return(res)
}


##Minque(1) approach with jackknife techniques used for a general mixed linear model
mixed.minque.jack=function(formula,data=list(),JacNum=NULL,JacRep=NULL){
   if(is.null(JacNum))JacNum=10
   if(is.null(JacRep))JacRep=1
   if(is.null(data))gdata=mixed.data(formula)
   else gdata=mixed.data(formula,data=data)
   res=genmod.jack(gdata,JacNum=JacNum,JacRep=JacRep)
   return(res)
}


mixed.genmod.minque=function(Y,Ped,Model,Cross){
   gdata=genmod.data(Y,Ped,Model,Cross)
   res=genmod(gdata)
   res$FixedEffect=res$FixEffect
   return(res)
}

mixed.genmod.minque.jack=function(Y,Ped,Model,Cross,JacNum=NULL,JacRep=NULL){
   if(is.null(JacNum))JacNum=10
   if(is.null(JacRep))JacRep=1
   gdata=genmod.data(Y,Ped,Model,Cross)
   res=genmod.jack(gdata,JacNum=JacNum,JacRep=JacRep)
   return(res)
}

#step_subset_simuq=function(G0,M,h2,REP,simunum,alpha=NULL){
#   if(is.null(alpha))alpha=0.01
#   mn=ncol(M)
#   n=length(G0)
#
#}

ginv=function (X, tol = sqrt(.Machine$double.eps)){ 

    if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X))) 
        stop("'X' must be a numeric or complex matrix")
    if (!is.matrix(X)) 
        X <- as.matrix(X)
    Xsvd <- svd(X)
    if (is.complex(X)) 
        Xsvd$u <- Conj(Xsvd$u)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
    if (all(Positive)) 
        Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
    else if (!any(Positive)) 
        array(0, dim(X)[2L:1L])
    else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
        t(Xsvd$u[, Positive, drop = FALSE]))
}

#step_subset_simu=function(G0,M,h2,REP,simunum,alpha=NULL){
#   if(is.null(alpha))alpha=0.01
#   mn=ncol(M)
#   n=length(G0)
#   vg=var(G0)
#   ce=sqrt((1/h2-1)*vg)
#   P=numeric(mn)
#   ID=list()
#   e=0
#   for(i in 1:simunum){
#      e=e+1
#      if(e==10){
#        cat("Simulation ",i,"/", simunum,"\n")
#        e=0
#      }
#      y=G0+ce*rnorm(n)
#      ptm <- proc.time()
#      id0=MDRSelect2(y,M,k=1,alpha=alpha)$ID[[1]]
#      M1=M[,id0]
#      if(length(id0)==1){
#         
#         #M1=as.matrix(M1)
#         #colnames(M1)=colnames(M)[id0]
#         id=id0
#      }
#      if(length(id0)>1){
#         res=step_subset(y,M1,alpha=alpha)
#         #res=forward_mdr_subset(y,M,k=k,REP,ALPHA=0.0001)
#         proc.time() - ptm
#         id=id0[res$idmax]
#      }
#      if(length(id)>0){
#         ID[[i]]=id
#         P[id]=P[id]+1
#      }
#   }
#   a=list(ID=ID,P=P)
#   return(a)
#}


Intersect=function(x,f){
   cls=unique(f)
   cn=length(cls)
   id=which(f==cls[1])
   v1=x[id]
   for(i in 2:cn){
      id=which(f==cls[i])
      v2=x[id]
      v1=intersect(v1,v2)
   }
   return(v1)
}

crecode=function(v){
  cname=unique(v)
  gn=length(cname)
  n=length(v)
  v1=numeric(n)
  for(i in 1:gn){
    id=which(v==cname[i])
    v1[id]=i
  }
  result=list(nv=v1,vname=cname)
  return(result)
}


#### A correlation function with dealing missing data#####
#y=Y[,1]

Cor=function(y,X){
   mn=ncol(X)
   r0=numeric(mn)
   p0=numeric(mn)
   for(i in 1:mn){
     #i=76
     x=X[,i]
     xy=data.frame(x,y)
     ok=complete.cases(xy)
     xy=xy[ok,]
     x=as.vector(xy$x)
     xy$x=as.numeric(x)
     x1<-xy[,1]                          # get rid of na's
     y1<-xy[,2]
     a=cor.test(x1,y1)
     r0[i]=a$estimate
     p0[i]=a$p.value
   }
   res=list(r0=r0,p0=p0)
   return(res)
}
#################################################################
##~~~Correlation between two groups~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#X=X12
mpcor=function(Y,X,...){
   Y=matrix(Y)
   TraitNum=ncol(Y)
   TraitNames=colnames(Y)
   mn=ncol(X)
   P=matrix(0,TraitNum,mn)
   R=matrix(0,TraitNum,mn)
   for(i in 1:TraitNum){
      #i=1
      res=Cor(Y[,i],X)
      R[i,]=res$r0
      P[i,]=res$p0
   }
   rownames(R)=rownames(P)=TraitNames
   colnames(R)=colnames(P)=colnames(X)

   result=list(R=R,P=P)
   return(result)
}
###Selection significant correlations ######
SelectCor=function(R,P,alpha=NULL,...){
  if(is.null(alpha))alpha=0.05
  TraitNum=nrow(R)
  TraitNames=rownames(R)
  MarkNames=colnames(R)
  ID=list(TraitNum)
  MarkID=list(TraitNum)
  RID=list(TraitNum)
  PID=list(TraitNum)
  Result=list(TraitNum)
  for(i in 1:TraitNum){
     id=which(P[i,]<=alpha)
     ID[[i]]=id
     MarkID[[i]]=MarkNames[id]
     RID[[i]]=R[i,id]
     PID[[i]]=P[i,id]
     Result[[i]]=cbind(id,MarkNames[id],R[i,id],P[i,id])
     colnames(Result[[i]])=c("ID","Mark","Corr","P_value")
  }
  names(Result)=TraitNames
  return(list(Result=Result, ID=ID))
}


##########################################################################################
## A function to generate a general linear function with class variables ################# 
factorlmform=function(data){
   formu=FactorLM(data)
   return(formu)
}

FactorLM=function(data){
  a=paste(names(data)[1],paste(paste("factor(",names(data)[-1],")"),collapse="+"),sep=" ~ ")
  a
}

###################################################################
#######Develop conditional variable################################

GetConVariable=function(Y){
   convar=GetConVar(Y)
   return(convar)
}

ConVariable=function(y,X){
   X=as.matrix(X)
   Cyx=cov(y,X)
   Cx=cov(X)
   k=Cyx%*%ginv(Cx)
   mux=apply(X,2,mean)
   x=t(X)-mux
   CY=y-t(k%*%x)
   return(CY)
}

GetConVar=function(Y){
   #require(gtools)
   TraitNames=colnames(Y)
   ComNum=length(TraitNames)-1
   e=1
   y=Y[,-(1:ComNum)]
   for(i in 1:ComNum){
      #i=1
      if(i==ComNum)C=matrix(1:i)
      else C=combn(ComNum,i)
      nr=ncol(C)
      name=NULL
      for(k in 1:nr){
          #k=1
          for(j in 1:i){
             if(j==1){
                 name[k]=TraitNames[C[j,k]]
             }
             else if(j>1){
                 name[k]=paste(name[k],TraitNames[C[j,k]],sep="&")
             }
          }
          x=Y[,C[,k]]
          x=as.matrix(x)
          if(e==1)ConY=ConVariable(y,x)
          else if(e>1)ConY=cbind(ConY,ConVariable(y,x))
          e=e+1
      }
      if(i==1)ConNames=name
      else if(i>1)ConNames=c(ConNames,name)
   }
   ConNum=length(ConNames)

   colnames(ConY)=paste(TraitNames[ComNum+1],ConNames,sep="|")
   ConY=cbind(ConY,y)
   colnames(ConY)[ConNum+1]=TraitNames[ComNum+1]
   con=list(Y=Y,ConY=ConY)
   return(con)
}


ConVariable1=function(y,X){
   data=data.frame(y,X)
   mod=lm(y~X,data=data)
   return(mod$residuals)
}

GetConVar1=function(Y){
   #require(gtools)
   TraitNames=colnames(Y)
   ComNum=length(TraitNames)-1
   e=1
   y=Y[,-(1:ComNum)]
   for(i in 1:ComNum){
      C=combn(ComNum,i)
      nr=nrow(C)
      name=NULL
      for(k in 1:nr){
          for(j in 1:i){
             if(j==1){
                 name[k]=TraitNames[C[k,j]]
             }
             else if(j>1){
                 name[k]=paste(name[k],TraitNames[C[k,j]],sep="&")
             }
          }
          x=Y[,C[k,]]
          x=as.matrix(x)
          if(e==1)ConY=ConVariable1(y,x)
          else if(e>1)ConY=cbind(ConY,ConVariable1(y,x))
          e=e+1
      }
      if(i==1)ConNames=name
      else if(i>1)ConNames=c(ConNames,name)
   }
   ConNum=length(ConNames)

   colnames(ConY)=paste(TraitNames[ComNum+1],ConNames,sep="|")
   ConY=cbind(ConY,y)
   colnames(ConY)[ConNum+1]=TraitNames[ComNum+1]
   con=list(Y=Y,ConY=ConY)
   return(con)
}


MINQUE=function(au=NULL,U,X=NULL,y){
    mk=length(U)
    if(is.null(au))au=rep(1,(mk+1))
    n=length(y)
    if(is.null(X))X=matrix(1,nrow=n,ncol=1)
    res=NULL
    Qa=GetQa(au,U,X)
    ML=MINQUE_L(Qa,U)
    MR=MINQUE_R(Qa,U,y)
    a=MINQUEVarPre(Qa,U,ML,MR,y)
    res$Var=a[[1]]
    res$FixedEffect=a[[2]]
    res$RandomEffect=a[[3]]
    res$X=X
    res$U=U
    res$au=au
    return(res)
}
#Qa=Qa0
#ML=ML0
#y=Y[,2]
#result[[3]]$Var[,22]
MINQUEVarPre=function(Qa,X,U,ML,MR,y,C=NULL,VI=NULL){
   svdinv=LargeInv
   #ML1=as.matrix(ML)
   #V0=ginv(ML)%*%MR ##Estimate variance components
   V0=svdinv(ML)%*%MR ##Estimate variance components
   v=as.vector(V0)
   mk=length(U)
   n=nrow(X)
   if(is.null(VI)){
      if(mk==0)VI=diag(1/v[1],n)
      else if(mk>0){
         id=which(v<=0)
         sumv=sum(v[id])
         smallvalue=sumv/10^3
         id=which(v<smallvalue)
         v[id]=smallvalue
         VI=Woodbury(v,U)
      }
   } 
   #X0=X
   if(ncol(X)==1){
       HF=t(X[,1])%*%VI%*%X[,1]
       bhat=svdinv(HF)%*%(t(X[,1])%*%VI%*%y)  ##Estimate fixed effects
   }
   if(is.null(C)){
       HF=t(X)%*%VI%*%X
       bhat=svdinv(HF)%*%(t(X)%*%VI%*%y)
   }
   else{
      HF=t(X)%*%VI%*%X+C%*%t(C)
      bhat=svdinv(HF)%*%(t(X)%*%VI%*%y)
   }
   #dim(X)
   #n=ncol(U[[1]])
   #mk=length(U)
   if(mk>0){
     k=numeric()
     for(i in 1:mk){
        nu=ncol(U[[i]])-1
        if(v[i]<=0){
          k[i]=0
          ##v[i]=0
        }
        else if(v[i]>0)k[i]=sqrt(nu*v[i]/MR[i])
     }
     PRE=numeric()
     vc=numeric()

     for(i in 1:mk){
        if(k[i]==0){
          es=rep(0,ncol(U[[i]]))
          PRE=c(PRE,es)
        }
        else if(k[i]>0){
          es<-k[i]*t(U[[i]])%*%Qa%*%y
          #es<-t(U[[i]])%*%Qa%*%y
          PRE=c(PRE,es)
        }
        #vc[i]=var(es)
     }
   }
   if(mk==0)result=list(v=v,b=bhat,Pre=NULL,HF=HF)
   else result=list(v=v,b=bhat,Pre=PRE,HF=HF)
   return(result)
}

MINQUEVarLUP=function(Qa,VI,X,U,ML,MR,y,C=NULL){
   svdinv=LargeInv
   V0=svdinv(ML)%*%MR ##Estimate variance components
   v=as.vector(V0)
   mk=length(U)
   if(ncol(X)==1) bhat=svdinv(t(X[,1])%*%VI%*%X[,1])%*%(t(X[,1])%*%VI%*%y)  ##Estimate fixed effects
   if(is.null(C)) bhat=svdinv(t(X)%*%VI%*%X)%*%(t(X)%*%VI%*%y)
   else bhat=svdinv(t(X)%*%VI%*%X+C%*%t(C))%*%(t(X)%*%VI%*%y)
   
   if(mk>0){
     k=numeric()
     for(i in 1:mk){
        nu=ncol(U[[i]])-1
        if(v[i]<=0){
          k[i]=0
          ##v[i]=0
        }
        else if(v[i]>0)k[i]=sqrt(nu*v[i]/MR[i])
     }
     PRE=numeric()
     vc=numeric()

     for(i in 1:mk){
        if(k[i]==0){
          es=rep(0,ncol(U[[i]]))
          PRE=c(PRE,es)
        }
        else if(k[i]>0){
          es<-k[i]*t(U[[i]])%*%Qa%*%y
          #es<-t(U[[i]])%*%Qa%*%y
          PRE=c(PRE,es)
        }
        #vc[i]=var(es)
     }
   }
   if(mk==0)result=list(v=v,b=bhat,Pre=NULL)
   else result=list(v=v,b=bhat,Pre=PRE)
   #result=list(v=v,b=bhat,Pre=NULL)
   return(result)
}



SVDInv=function (x, eps = 1e-06) 
{
    if (any(is.na(x))) 
        stop("NA(s) encountered in matrix, cannot invert.")
    if (length(x) > 1) {
        savesvd <- svd(x, LINPACK = TRUE)
        U.svd <- savesvd$u
        V.svd <- savesvd$v
        d.svd <- savesvd$d
        maxd <- max(d.svd)
        w <- ifelse((d.svd/maxd) < eps, rep(0, length(d.svd)), 
            1/d.svd)
        rank <- sum(d.svd/maxd >= eps)
        Ginv <- V.svd %*% diag(w) %*% t(U.svd)
    }
    else {
        Ginv <- ifelse(x < eps, 0, 1/x)
        rank <- ifelse(x < eps, 0, 1)
    }
    #list(Ginv = Ginv, rank = rank)
    return(Ginv)
}

    
#SVDInv=function(A){
#  ##A=vi
#  ##s=svd(A,LINPACK = FALSE)
#  s=svd(A,LINPACK = T)
#  r=nrow(A)
#  SMALLVALUE=1e-10
#  D=matrix(0,nrow=r,ncol=r)
#  for(i in 1:r){
#    if(s$d[i]<=SMALLVALUE)D[i,i]=SMALLVALUE
#    else if(s$d[i]>SMALLVALUE)D[i,i]=1/s$d[i]
#  }
#  A1=s$v%*%D%*%t(s$u)
#  return(A1)
#}

   
GetQa=function(au,U,X,C=NULL){
  svdinv=LargeInv
  mk=length(U)
  n=nrow(X)
  if(mk>0)VI=Woodbury(au,U)

  else if(mk==0)VI=diag(1/au,n)
  if(ncol(X)==1)VX=VI%*%X[,1]
  else VX=VI%*%X
  if(is.null(C))Qa=VI-VX%*%ginv(t(X)%*%VX)%*%t(VX)
  else Qa=VI-VX%*%svdinv(t(X)%*%VX+C%*%t(C))%*%t(VX)
  res=list(Qa=Qa,VI=VI)
  return(res)

}

## Calculate the left side of MINQUE normal equations

MINQUE_L=function(Qa,U){
  mk=length(U)
  ML=matrix(0, nrow=mk+1,ncol=mk+1)
  if(mk>0){
  for(i in 1:mk){
     for(j in i:mk){
        a1=t(U[[i]])%*%Qa%*%U[[j]]
        a=a1%*%t(a1)
        b=sum(diag(a))
        ML[i,j]=ML[j,i]=b[1]
     }
     a1=t(U[[i]])%*%Qa
     a=a1%*%t(a1)
     b=sum(diag(a))
     ML[i,(mk+1)]=ML[(mk+1),i]=b[1]
  }
  }

  a1=Qa
  a=a1%*%t(a1)
  b=sum(diag(a))
  ML[(mk+1),(mk+1)]=b[1]
   
  return(ML)
}
#class(y)
## calculate the right side of MINQUE normal equations
#y=Y[,1]
MINQUE_R=function(Qa,U,y){
  mk=length(U)
  MR<-numeric()
  y=as.numeric(y)
  if(mk>0){
     for(i in 1:mk){
       #i=1
       t=t(y)%*%Qa%*%U[[i]]
       a=t%*%t(t)
       MR[i]=a[1,1]
     }
  }
  t=t(y)%*%Qa
  a=t%*%t(t)
  MR[mk+1]=a[1,1]
  return(MR)
}

##y=YD[,1]
GetVarPre=function(Qa,U,X,ML,y){
   svdinv=LargeInv
   MR=MINQUE_R(Qa,U,y)
   V0=svdinv(ML)%*%MR


   v=as.vector(V0)
   mk=length(U)
   n=nrow(X)
   if(mk==0)VI=diag(1,n)
   else VI=Woodbury(v,U) 
   X0=as.matrix(X)                      
   ##Estimate fixed effects by GLSE
   if(ncol(X)==1) bhat=ginv(t(X[,1])%*%VI%*%X[,1])%*%(t(X[,1])%*%VI%*%y)  ## I don't know why, but it works this way
   else if(ncol(X)>1)bhat=ginv(t(X0)%*%VI%*%X0)%*%(t(X0)%*%VI%*%y)
   
   #n=ncol(U[[1]])
   #mk=length(U)
   #v=as.vector(V0)
   if(mk>0){
     k=numeric()
     for(i in 1:mk){
       nu=ncol(U[[i]])-1
       if(v[i]<=0){
          k[i]=0
          ##v[i]=0
       }
       else if(v[i]>0)k[i]=sqrt(nu*v[i]/MR[i])
     }
     PRE=numeric()
     vc=numeric()

     for(i in 1:mk){
      if(k[i]==0){
          es=rep(0,ncol(U[[i]]))
          PRE=c(PRE,es)
      }
      else if(k[i]>0){
         es<-k[i]*t(U[[i]])%*%Qa%*%y
         #es<-t(U[[i]])%*%Qa%*%y
         PRE=c(PRE,es)
      }
      #vc[i]=var(es)
    }
   }
   if(mk==0)PRE=NULL
   result=list(v=v,b=bhat,Pre=PRE)
   return(result)
}

GetFit=function(U,X,b,Pre){
     fit=X%*%b
     #Pre=PRE0
     mk=length(U)
     BigU=U[[1]]
     for(i in 2:mk)BigU=cbind(BigU,U[[i]])
     fit=fit+BigU%*%Pre
     return (fit)
     #dim(BigU)
}

#######################################################################################
####### Function to Drop the columns of any matrix in which all elements are zero.#####
#######################################################################################

REML=function(U,X,y,ITMAX=NULL,...){
    #U=U1
    #y=yd
    if(is.null(ITMAX))ITMAX=5
    SMALLVALUE=0.001
    DIS=100
    it=1
    mk=length(U)
    v0=rep(1,mk+1)
    S=which(v0<=0)
    v0[S]=0

    while(DIS>SMALLVALUE&&it<=ITMAX){
   
        au=v0
        Qa=GetQa(au,U,X)$Qa
        ML=MINQUE_L(Qa,U)
        result=GetVarPre(Qa,U,ML,y)
        v1=result[[1]]
        S=which(v1<=0)
        v1[S]=0
        if(it==1)ITV=v1
        else if(it>1)ITV=cbind(ITV,v1)
        it=it+1
        d=v0-v1
        DIS=sqrt(mean(d^2)/sum(v0))
        v0=v1
    }
    Pre=result[[2]]
    reml=list(v1,Pre)
    return(reml)
}

Woodbury=function(au=NULL,U,...){
   svdinv=LargeInv
   mk=length(U)
   if(is.null(au))au=rep(1,mk+1)
   n=nrow(U[[1]])
   V=diag(1/(au[mk+1]),n)
   for(i in 1:mk){
     a=au[i]
     if(a>0){
        c=ncol(U[[i]])
        m1=V%*%U[[i]]
        a=au[i]
        d=diag(1,c)
        vi=d+a*t(U[[i]])%*%m1
        vi=svdinv(vi)
        VI=V-a*m1%*%vi%*%t(m1)
        V=VI
     }
   }
   VI=V
   return(VI)
}
#U=gdata$U
Woodbury0=function(au=NULL,U){
   svdinv=LargeInv
   mk=length(U)
   if(is.null(au))au=rep(1,mk+1)
   n=nrow(U[[1]])
   V=diag(au[mk+1],n)
   for(i in 1:mk){
     a=au[i]
     if(a>0)V=V+a*U[[i]]%*%t(U[[i]])
   }
   V=svdinv(V)
   return(V)
}

SimuData0=function(U,v,X=NULL,b,SimuNum){
  #v=v0
  mk=length(U)
  n=nrow(U[[1]])
  tc=numeric()
  if(is.null(X))X=matrix(1,n,1)
  for(i in 1:mk){
     tc[i]=ncol(U[[i]])
     if(i==1)BigU=U[[i]]
     else if(i>1)BigU=cbind(BigU,U[[i]])
  }
  mt=sum(tc)
  RE0=matrix(0,nrow=mt,ncol=SimuNum)
  YS=matrix(0,nrow=n,ncol=SimuNum)
  b=as.matrix(b)
  for(i in 1:SimuNum){
     for(j in 1:mk){
        if(v[j]>0)a=rnorm(c[j],0,sqrt(v[j]))
        else if(v[j]<=0)a=rep(0,c[j])
        if(j==1)RE=a
        else if(j>1)RE=c(RE,a)
     }
     RE0[,i]=RE
     YS[,i]=X%*%b+BigU%*%RE+rnorm(n,0,sqrt(v[mk+1]))
     
  }
  res=list(SimuY=YS,SimuEffect=RE0)
  return(YS)
}
SimuData1=function(gdata,v,b,SimuNum){
  #v=v0
  U=gdata$U
  nk=length(gdata$X)
  if(nk==1)X=gdata$X[[1]]
  else{
    for(i in 1:nk){
       if(i==1)X=gdata$X[[i]]
       else X=cbind(X,gdata$X[[i]])
    }
  }
  mk=length(U)
  n=nrow(U[[1]])
  c=numeric()
  for(i in 1:mk){
     c[i]=ncol(U[[i]])
     if(i==1)BigU=U[[i]]
     else if(i>1)BigU=cbind(BigU,U[[i]])
  }
  b=as.matrix(b)
  mt=sum(c)
  RE0=matrix(0,nrow=mt,ncol=SimuNum)
  YS=matrix(0,nrow=n,ncol=SimuNum)
  for(i in 1:SimuNum){
     for(j in 1:mk){
        if(v[j]>0)a=rnorm(c[j],0,sqrt(v[j]))
        if(v[j]<=0)a=rep(0,c[j])
        if(j==1)RE=a
        else if(j>1)RE=c(RE,a)
     }
     RE0[,i]=RE
     YS[,i]=X%*%b+BigU%*%RE+rnorm(n,0,sqrt(v[mk+1]))
     
  }
  res=list(SimuY=YS,SimuEffect=RE0)
  gdata$Y=YS
  gdata$VC=rep(1,mk+1)
  return(gdata)
}

BigM=function(gdata,...){

  U=gdata$U
  nk=length(gdata$X)
  if(nk==1)X=gdata$X[[1]]
  else{
    for(i in 1:nk){
       if(i==1)X=gdata$X[[i]]
       else X=cbind(X,gdata$X[[i]])
    }
  }
  mk=length(U)
  n=nrow(U[[1]])
  tc=numeric()
  for(i in 1:mk){
     tc[i]=ncol(U[[i]])
     if(i==1)BigU=U[[i]]
     else if(i>1)BigU=cbind(BigU,U[[i]])
  }
  mt=sum(tc)
  mat=list(BigU=U,X=X,mt=mt)
  return(mat)

}

simudata=function(gdata,v,b,SimuNum=NULL){
   if(is.null(SimuNum))SimuNum=200
   gdata=SimuData_Old(gdata,v,b,SimuNum)
   return(gdata)
}
   
SimuData=function(gdata,v,b,SimuNum=NULL){
   if(is.null(SimuNum))SimuNum=200
   gdata=SimuData_Old(gdata,v,b,SimuNum)
   return(gdata)
}

SimuData_Old=function(gdata,v,b,SimuNum=NULL){
  if(is.null(SimuNum))SimuNum=200
  U=gdata$U
  nk=length(gdata$X)
  if(nk==1)X=gdata$X[[1]]
  else{
    for(i in 1:nk){
       if(i==1)X=gdata$X[[i]]
       else X=cbind(X,gdata$X[[i]])
    }
  }
  mk=length(U)
  n=nrow(U[[1]])
  c=numeric()
  for(i in 1:mk){
     c[i]=ncol(U[[i]])
     if(i==1)BigU=U[[i]]
     else if(i>1)BigU=cbind(BigU,U[[i]])
  }
  dim(X)
  mt=sum(c)
  #b=as.matrix(b)
  RE0=matrix(0,nrow=mt,ncol=SimuNum)
  YS=matrix(0,nrow=n,ncol=SimuNum)
  for(i in 1:SimuNum){
     #i=1
     for(j in 1:mk){
        #j=1
        if(v[j]>0)a=rnorm(c[j],0,sqrt(v[j]))
        if(v[j]<=0)a=rep(0,c[j])
        if(j==1)RE=a
        else if(j>1)RE=c(RE,a)
     }
     RE0[,i]=RE
     #RE=as.matrix(RE)
     YS[,i]=X%*%b+BigU%*%RE+rnorm(n,0,sqrt(v[mk+1]))
     
  }
  #res=list(SimuY=YS,SimuEffect=RE0)
  gdata$Y=YS
  gdata$SimuEffect=RE0
  gdata$VC=rep(1,mk+1)
  return(gdata)
}


RandJackCond=function(U,X,Y,VC,VNames,JacNum,JacRep,Method){
   n=nrow(Y)
   TraitNum=ncol(Y)
   TraitNames=colnames(Y)
   JB=rep(1:JacNum,length=n)
   #JacNum=length(JB)
   a=1-(1-0.05)^(1/(JacNum-2))
   t=qt(1-a/2,JacNum)
   mk=length(U)
   
   ml=mk+1
   mt=0
   for(i in 1:mk)mt=mt+ncol(U[[i]])
   Str=NULL
   for(i in 1:mk)Str=c(Str,colnames(U[[i]]))
   RES0=NULL
   VPower=numeric(ml)
   V0=matrix(0,TraitNum,ml)
   P0=matrix(0,TraitNum,mt)
   if(Method==1){
      au=rep(1,ml)
      Qa0=GetQa(au,U,X)$Qa
      ML0=MINQUE_L(Qa0,U)
      for(i in 1:TraitNum){
         res=GetVarPre(Qa0,U,ML0,Y[,i])
         V0[i,]=res[[1]]*VC
         P0[i,]=res[[2]]
      }
   }
   else if(Method==0){
       for(i in 1:TraitNum){
          res=REML(U,X,Y[,i],5)
          V0[i,]=res[[1]]*VC
          P0[i,]=res[[2]]
       }
   }
   #list(V0,P0)
   
   Vp=matrix(0,TraitNum,ml)
   Cv=matrix(0,TraitNum,ml)
   Cv1=matrix(0,TraitNum,ml)
   Cv2=matrix(0,TraitNum,ml)
   Vm=matrix(0,TraitNum,ml)

   CRp=matrix(0,TraitNum-1,ml)
   CRv=matrix(0,TraitNum-1,ml)
   CRv1=matrix(0,TraitNum-1,ml)
   CRv2=matrix(0,TraitNum-1,ml)
   CRm=matrix(0,TraitNum-1,ml)


   VCp=matrix(0,TraitNum,ml)
   VCv=matrix(0,TraitNum,ml)
   VCv1=matrix(0,TraitNum,ml)
   VCv2=matrix(0,TraitNum,ml)
   VCm=matrix(0,TraitNum,ml)


   Pp=matrix(0,TraitNum,mt)
   Cp=matrix(0,TraitNum,mt)
   Cp1=matrix(0,TraitNum,mt)
   Cp2=matrix(0,TraitNum,mt)
   Pm=matrix(0,TraitNum,mt)
   
   CRPp=matrix(0,TraitNum-1,mt)
   CRCp=matrix(0,TraitNum-1,mt)
   CRCp1=matrix(0,TraitNum-1,mt)
   CRCp2=matrix(0,TraitNum-1,mt)
   CRPm=matrix(0,TraitNum-1,mt)


   
   VJR=matrix(0,nrow=ml,ncol=JacRep)
   JacVar=matrix(0,nrow=JacNum,ncol=ml)
   JacPre=matrix(0,nrow=JacNum,ncol=mt)
   
   JV=matrix(0,nrow=ml,ncol=JacNum)
   JP=matrix(0,nrow=mt,ncol=JacNum)

   df=JacNum-1

   JACVAR=matrix(0,TraitNum,ml)
   JACPRE=matrix(0,TraitNum,mt)

   for(k in 1:JacRep){
      #k=1
      #e1=1
      U1=NULL
      JacRes=NULL
      
      JAC1=NULL
      JAC2=NULL
      #N1=sample(1:n,replace=F)
      JB=sample(JB)
      for(i in 1:JacNum){
         #e2=e1+JB[i]-1
         #index=N1[e1:e2]
         index=which(JB==i)
         for(j in 1:mk){
            m=U[[j]]
            m=m[-index,]
            U1[[j]]=m
         }
         YD=Y[-index,]
         X1=X[-index,]
         if(Method==1){
             Qa=GetQa(au,U1,X1)$Qa
             ML=MINQUE_L(Qa,U1)
             for(j in 1: TraitNum){
                 res=GetVarPre(Qa,U1,ML,YD[,j])
                 JACVAR[j,]=res[[1]]*VC
                 JACPRE[j,]=res[[2]]
             }
         }
         else if(Method==0){
            for(j in 1: TraitNum){
                res=REML(U1,X1,YD[,j],5)
                JACVAR[j,]=res[[1]]*VC
                JACPRE[j,]=res[[2]]
            }

         }
         JAC1[[i]]=JACVAR
         JAC2[[i]]=JACPRE
         #e1=e2+1
      }
   
      #for(i in 1: JacNum){
      #   JV[,i]=JacNum*v0-df*JacVar[,i]
      #   JP[,i]=JacNum*pre0-df*JacPre[,i]
      #}
      #########################################################
      ## non-psudeo value jackknife estimate###################
      for(s in 1:TraitNum){
          #s=1
          for(j in 1:JacNum){
             JacVar[j,]=JAC1[[j]][s,]
             JacPre[j,]=JAC2[[j]][s,]
          }
          for(j in 1:JacNum){
             index=which(JacVar[j,]<0)
             JacVar[j,index]=0
          }
          if(s==TraitNum){
              JVY=JacVar
              JPY=JacPre
          }
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
                
          for(i in 1:mt){
             m=mean(JacPre[,i])
             se=sqrt(var(JacPre[,i]))
             t1=m/se
             p=pt(-abs(t1),df)
             #p=a$p.value
             p=1-(1-p)^(JacNum-2)
             Pp[s,i]=Pp[s,i]+p
             Pm[s,i]=Pm[s,i]+m
             Cp[s,i]=Cp[s,i]+se
          }
      }
      
      for(s in 1:(TraitNum-1)){
          
          for(j in 1:JacNum){
             JacVar[j,]=JAC1[[j]][s,]
             JacPre[j,]=JAC2[[j]][s,]
          }
          for(j in 1:JacNum){
             index=which(JacVar[j,]<0)
             JacVar[j,index]=0
          }
          for(i in 1:ml){
             pr=prop(JacVar[,i],JVY[,i])
             m=mean(pr)
             se=sqrt(var(pr))
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
             CRp[s,i]=CRp[s,i]+p
             CRm[s,i]=CRm[s,i]+m
             CRv[s,i]=CRv[s,i]+se
          } 
                
          for(i in 1:mt){
             cr=JPY[,i]-JacPre[,i]
             m=mean(cr)
             se=sqrt(var(cr))
             t1=m/se
             p=pt(-abs(t1),df)
             #p=a$p.value
             p=1-(1-p)^(JacNum-2)
             CRPp[s,i]=Pp[s,i]+p
             CRPm[s,i]=Pm[s,i]+m
             CRCp[s,i]=Cp[s,i]+se
          }
      }


   }
   
   Vm=Vm/JacRep
   Vp=Vp/JacRep
   Cv=Cv/JacRep
   Cv1=Vm-t*Cv
   Cv2=Vm+t*Cv

   CRm=CRm/JacRep
   CRp=CRp/JacRep
   CRv=CRv/JacRep
   CRv1=CRm-t*CRv
   CRv2=CRm+t*CRv


   VCm=VCm/JacRep
   VCp=VCp/JacRep
   VCv=VCv/JacRep
   VCv1=VCm-t*VCv
   VCv2=VCm+t*VCv

   
   Pm=Pm/JacRep
   Pp=Pp/JacRep
   Cp=Cp/JacRep
   Cp1=Pm-t*Cp
   Cp2=Pm+t*Cp
   
   CRPm=CRPm/JacRep
   CRPp=CRPp/JacRep
   CRCp=CRCp/JacRep
   CRCp1=CRPm-t*CRCp
   CRCp2=CRPm+t*CRCp


   VAR=NULL
   PRE=NULL
   VC=NULL
   CR=NULL
   CP=NULL
   CNames=c("Estimate","SE","PValue","CL","CU")
   for(i in 1: TraitNum){
      VAR[[i]]=data.frame(Vm[i,],Cv[i,],Vp[i,],Cv1[i,],Cv2[i,])
      VC[[i]]=data.frame(VCm[i,],VCv[i,],VCp[i,],VCv1[i,],VCv2[i,])
      PRE[[i]]=data.frame(Pm[i,],Cp[i,],Pp[i,],Cp1[i,],Cp2[i,])
      if(i<TraitNum){
         CR[[i]]=data.frame(CRm[i,],CRv[i,],CRp[i,],CRv1[i,],CRv2[i,])
         CP[[i]]=data.frame(CRPm[i,],CRCp[i,],CRPp[i,],CRCp1[i,],CRCp2[i,])
         rownames(CR[[i]])=VNames
         colnames(CR[[i]])=CNames
         rownames(CP[[i]])=Str
         colnames(CP[[i]])=CNames

      }
      rownames(VAR[[i]])=VNames
      colnames(VAR[[i]])=CNames

      

      rownames(VC[[i]])=paste(VNames,"/VP",sep="")
      colnames(VC[[i]])=CNames
      rownames(PRE[[i]])=Str
      colnames(PRE[[i]])=CNames
      
      

   }
   
   names(VAR)=TraitNames
   names(VC)=TraitNames
   names(PRE)=TraitNames
   names(CR)=TraitNames[1:(TraitNum-1)]
   names(CR)=TraitNames[1:(TraitNum-1)]
   a=list(VAR=VAR,VARP=VC,CRaitio=CR,PRE=PRE,CRPre=CP)
   return(a)
   
  
} 

prop=function(x,y){
   n=length(x)
   c=numeric(n)
   for(i in 1:n){
     a=x[i]
     b=y[i]
     if(a<=0)a=0
     if(b<=0)c[i]=0
     else if(b>0)c[i]=1-a/b
  }
  return(c)
}

vectcomp=function(v0,v1){

  name0=names(v0)
  name1=names(v1)
  vn=length(v0)
  v2=numeric(vn)
  for(i in 1:vn){
     id=which(name0[i]==name1)
     if(length(id)>0)v2[i]=v1[id]
  }
  names(v2)=name0
  return(v2)
}
    
svdwrapper = function( x, nu, nv, verbose=T ){ 
  #   Caution: I have not tested this much. 
  #   It's here as an example for an R-Help discussion. 
    gotit = F 
    try( {svdx = svd(x,nu,nv); gotit=T}, silent = !verbose ) 
    if( gotit )return(svdx) 
    try( {svdtx = svd(t(x),nv,nu); gotit=T}, silent = !verbose ) 
    if( !gotit )stop("svd(x) and svd(t(x)) both failed.") 
    if( verbose )print("svd(x) failed but svd(t(x)) worked.") 
    temp    = svdtx$u 
    svdtx$u = svdtx$v 
    svdtx$v = temp 
    svdtx 
} 

svdinv=function(A,...){
   return(LargeInv(A))
}

svdinv0=function(x,...) 
{
    eps = 1e-06
    verbose=T
    if (any(is.na(x))) 
        stop("NA(s) encountered in matrix, cannot invert.")
    if (length(x) > 1) {
        gotit=F
        try( {svdx = svd(x); gotit=TRUE}, silent = !verbose ) 
        if(gotit)savesvd=svdx
        if(!gotit){
           try( {svdtx=svd(t(x));gotit=T},silent = !verbose ) 
           if( !gotit)stop("svd(x) and svd(t(x)) both failed.") 
           if( verbose)print("svd(x) failed but svd(t(x)) worked.") 
           temp    = svdtx$u 
           svdtx$u = svdtx$v 
           svdtx$v = temp 
           savesvd=svdtx
        }
        #savesvd <- svd(x, LINPACK = FALSE)
        U.svd <- savesvd$u
        V.svd <- savesvd$v
        d.svd <- savesvd$d
        maxd <- max(d.svd)
        w <- ifelse((d.svd/maxd) < eps, rep(0, length(d.svd)), 
            1/d.svd)
        rank <- sum(d.svd/maxd >= eps)
        Ginv <- V.svd %*% diag(w,) %*% t(U.svd)
    }
    else {
        Ginv <- ifelse(x < eps, 0, 1/x)
        rank <- ifelse(x < eps, 0, 1)
    }
    #list(Ginv = Ginv, rank = rank)
    return(Ginv)
}

##A function used for Monte Carlo simulations
genmod.simu=function(gdata,v,b,SimuNum,JacNum=NULL,JacRep=NULL,ALPHA=NULL,...){
   if(is.null(JacNum))JacNum=10
   if(is.null(JacRep))JacRep=1
   if(is.null(ALPHA))ALPHA=0.05
   
   result=genmod.simuold(gdata,v,b,SimuNum,JacNum,JacRep,ALPHA)
   return(result)
}

##A function used for Monte Carlo simulations

genmod.simuold=function(gdata,v,b,SimuNum=NULL,JacNum=NULL,JacRep=NULL,ALPHA=NULL,...){
   if(is.null(SimuNum))SimuNum=200
   if(is.null(ALPHA))ALPHA=0.05
   gdata=SimuData(gdata,v=v,b=b,SimuNum=SimuNum)
   if(is.null(JacNum))JacNum=10
   if(is.null(JacRep))JacRep=1
   
   jac=genmod.jack(gdata,JacNum=JacNum,JacRep=JacRep,LUP=1)
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

mixed.data1=function(formula,data=list(),...){
   
   if(is.null(data))result=mixed.data1(formula)
   else{
       nc=ncol(data)
       for(i in 1:nc){
           v=data[,i]
           if(class(v)=="factor")v=as.character(v)
           data[,i]=v
       }
       result=mixed.data1(formula,data)
   }
   return(result)
}

mixed.data2=function(formula,data=list(),...){
   
   if(is.null(data))result=mixed.data1(formula)
   else {
       nc=ncol(data)
       for(i in 1:nc){
          v=as.vector(data[,i])
          if(i==1)dat1=v
          else dat1=cbind(dat1,v)
       }
       colnames(dat1)=colnames(data)
       dat1=as.data.frame(dat1)
       result=mixed.data1(formula,dat1)
   }
   return(result)
}


mixed.model=function(formula,data=list(),...){
   if(is.null(data))return(mixed.data(formula=formula))
   else return(mixed.data(formula=formula,data=data))
}

mixed.data=function(formula,data=list(),...){
   if(is.null(data))mq=lm.data(formula=formula)
   else mq=lm.data(formula=formula,data=data)
   Y=mq$Y
   X=mq$X
   U=mq$U
   xk=length(X)
   nx=0
   for(i in 1:xk)nx=nx+ncol(X[[i]])
   if(xk==1)C=NULL
   if(xk>1){
      C=matrix(0,nx,xk-1)
      e1=ncol(X[[1]])+1
      for(i in 2:xk){
         e2=e1+ncol(X[[i]])-1
         C[(e1:e2),(i-1)]=1
         e1=e2+1
      }
   }
      
   #for(i in 1:xk)X=cbind(X,mq$X[[i]])
   mk=length(U)
   VC=rep(1,(mk+1))
   if(class(Y)=="numeric")Y=data.frame(Y)
   result=list(Y=Y,X=X,U=U,VC=VC,C=C,Model=formula)
   result$call=match.call()
   class(result)="genmod.data"
   return(result)
}

##genmod(gdata)
##SIMUYES=1

genmod=function(gdata,au=NULL,LUP=NULL,...){
   Y=gdata$Y
   U=gdata$U
   C=gdata$C
   ml=length(gdata$U)+1
   if(is.null(au))au=rep(1,ml)
   xk=length(gdata$X)
   if(xk==1)X=gdata$X[[1]]
   if(xk>1)X=BigX(gdata$X)
   nx=ncol(X)
   #for(i in 1:nx)X[,i]=as.numeric(X[,i])
   VC=gdata$VC
   #Model=gdata$model
   #BigX=NULL
   #class(X[,3])
   #X=as.matrix(X)
   if(is.null(LUP))LUP=NULL
   #result=minque.default(X,U,Y,C,au)
   result=minque.default(X,U,Y,C,au,LUP)
   result$Var=result$Var*VC
   #result$call=match.call()
   class(result)="minque"
   return(result)
}

genmod.jack2=function(gdata,JacNum,JacRep,ALPHA=NULL,LUP=NULL,...){
   Y=gdata$Y
   dim(Y)
   U=gdata$U
   xk=length(gdata$X)
   if(xk==1)X=gdata$X[[1]]
   else X=BigX(gdata$X)
   #X=gdata$X
   VC=gdata$VC
   #if(is.null(gdata$C))
   C=gdata$C
   #Model=gdata$model
   if(is.null(ALPHA))ALPHA=0.05
   if(is.null(LUP))LUP=NULL
   #result=minque.jack1(Y,X,U,JacNum,JacRep,VC,C,ALPHA)
   #else 
   result=minque.jack1(Y,X,U,JacNum,JacRep,VC,C,ALPHA,LUP)

   #result=minque.jack(X,U,Y)
   #result$Var=result$Var*VC
   result$call=match.call()
   class(result)="minque.jack"
   return(result)
}

#gdata=mod
print.genmod.data=function(gdata,...){
   cat("\nCall:\n")
   print(gdata$call)

   cat("\nThe information data are very complicated. We only provide summarized information as follows:\n")
   
   cat("\nSample size=", nrow(gdata$Y),"\n")
   cat("\nNo of traits=", ncol(gdata$Y),"\n")
  
   cat("\nNo of information matrices X for fixed effects is: ", length(gdata$X), "\n")
   mk=length(gdata$U)
   if(mk==0) cat("\nNo matrix for ramdom effects available\n")
   else if (mk>0) cat("\nNo of information matrices for random effects is:", length(gdata$U), "\n")
   
}



#minque=function(Y,formula,...) UseMethod("minque")

minque.resample=function(Y,formula,JacNum,JacRep,data=list(),method=NULL,...){
   if(is.null(data))mq=lm.data(formula=formula)
   else mq=mixed.data(formula=formula,data=data)
   X=mq$Xmatrix
   U=mq$Umatrix
   mq$fid
   #res0=minque.default(X,U,Y)
   result=minque.jack(Y,X,U,JacNum,JacRep)
   #result=list(res0=res0,res.jack=jack)
   #result=list(X=X,U=U,Y=Y)
   result$call=match.call()
   class(result)="minque.jack"
   return(result)
}


minque=function(Y,formula,data=list(),method=NULL){
   if(is.null(data))mq=mixed.data(formula=formula)
   else mq=mixed.data(formula=formula,data=data)
   X=mq$Xmatrix
   U=mq$Umatrix
   mq$fid
   res0=minque.default(X,U,Y)
   #jack=minque.jack(Y,X,U,10,5)
   #result=list(res0=res0,res.jack=jack)
   #result=list(X=X,U=U,Y=Y)
   res0#call=match.call()
   class(res0)="minque"
   return(res0)
}


minque.default=function(X,U,Y,C=NULL,au=NULL,LUP=NULL,...){
    mk=length(U)
    n=nrow(U[[1]])
    if(is.null(au))au=rep(1,(mk+1))
    if(is.null(X))X=matrix(1,nrow=n,ncol=1)
    if(is.null(C))C=NULL
    xnames=colnames(X)
    if(mk==0)vnames="V(e)"
    else if(mk>0){
       rnames=NULL
       vnames=c(paste("V(",names(U),")",sep=""),"V(e)")
       for(i in 1:mk)rnames=c(rnames,colnames(U[[i]]))
    }

    aa=GetQa(au,U,X,C)
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
    if(is.null(LUP))VI=NULL
    for(i in 1:TraitNum){
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


minque.jack=function(Y,X,U,JacNum,JacRep,VC=NULL,C=NULL,VNames=NULL,Method=NULL,ALPHA=NULL,LUP=NULL,...){
   if(is.null(LUP))LUP=NULL
   #result=minque.jack1(Y,X,U,JacNum,JacRep,VC,C,VNames,Method,ALPHA)
   #else 
   result=minque.jack1(Y,X,U,JacNum,JacRep,VC,C,VNames,Method,ALPHA,LUP)
   return(result)
}

minque.jack1=function(Y,X,U,JacNum,JacRep,VC=NULL,C=NULL,VNames=NULL,Method=NULL,ALPHA=NULL,LUP=NULL,...){
   n=nrow(Y)
   mk=length(U)
   if(is.null(VC))VC=rep(1,mk+1)
   if(is.null(VNames)){
      if(mk>0)VNames=c(paste("Var(",names(U),")",sep=""),"Ve")
      else VNames="Ve"
   }
   if(is.null(C))C=NULL
   SIMUYES=NULL
   
   if(is.null(Method))Method=1
   TraitNum=ncol(Y)
   TraitNames=colnames(Y)
   nx=ncol(X)
   XNames=colnames(X)
   if(length(XNames)==0)XNames=paste("B<",(1:nx),">",sep="")
   JB=rep(1:JacNum,length=n)
   if(is.null(ALPHA))ALPHA=0.05
   a=1-(1-ALPHA)^(1/(JacNum-2))
   t=qt(1-a/2,JacNum-1)
   
   #VC=vc
   ml=mk+1
   mt=0
   if(mk>0){
     for(i in 1:mk)mt=mt+ncol(U[[i]])
     Str=NULL
     for(i in 1:mk)Str=c(Str,colnames(U[[i]]))
   }

   RES0=NULL
   VPower=numeric(ml)
   V0=matrix(0,TraitNum,ml)
   if(mk>0)P0=matrix(0,TraitNum,mt)
   B0=matrix(0,TraitNum,ncol(X))
   #dim(U[[5]])

   if(Method==1){
      au=rep(1,ml)
      aa=GetQa(au,U,X,C)
      Qa0=aa$Qa
      if(is.null(SIMUYES)==F)VI0=aa$VI
      ML0=MINQUE_L(Qa0,U)
      if(is.null(SIMUYES)){  ##LUP method
         for(i in 1:TraitNum){
            #i=2 sum(Y[,i])
            MR=MINQUE_R(Qa0,U,Y[,i])
            res=MINQUEVarPre(Qa0,X,U,ML0,MR,Y[,i],C)
         
            #(Qa,X,U,ML,MR,y,C=NULL)
            #res=MINQUEVarPre(Qa0,X,U,ML0,MY[,i],C)
            V0[i,]=res$v*VC
            B0[i,]=res$b
            if(mk>0)P0[i,]=res$Pre
         }
     }
     else{                   ##AUP method
         for(i in 1:TraitNum){
            MR=MINQUE_R(Qa0,U,Y[,i])
            res=MINQUEVarPre(Qa0,VI0,X,U,ML0,MR,Y[,i],C)
            V0[i,]=res$v*VC
            B0[i,]=res$b
            #if(mk>0)P0[i,]=res$Pre
         }
     }
     
   }
   else if(Method==0){
       for(i in 1:TraitNum){
          res=REML(U,X,Y[,i],5)
          V0[i,]=res[[1]]*VC
          if(mk>0)P0[i,]=res[[3]]
          B0[i,]=res[[2]]
       }
   }
   #list(V0,P0)
   
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

   for(k in 1:JacRep){
      #k=1
      U1=NULL
      JacRes=NULL
      JAC1=NULL
      JAC2=NULL
      JAC3=NULL
      JB=sample(JB)
      for(i in 1:JacNum){
         #i=1
         index=which(JB==i)
         if(mk==0)U1=NULL
         else if (mk>0){
            for(j in 1:mk){
              m=U[[j]]
              m=m[-index,]
              U1[[j]]=m
            }
         }
         YD=Y[-index,]
         if(TraitNum==1)YD=matrix(YD,nrow=length(YD),ncol=1) ##or YD=data.frame(YD)
         #else if(TraitNum>1)YD=Y[-index,]
         if(nx==1)X1=as.matrix(X[-index,])
         else X1=X[-index,]
         if(Method==1){
             aa=GetQa(au,U1,X1,C)
             Qa1=aa$Qa
             if(is.null(SIMUYES)==F)VI1=aa$VI
             ML1=MINQUE_L(Qa1,U1)
             if(is.null(SIMUYES)){
                for(j in 1: TraitNum){
                    
                    MR1=MINQUE_R(Qa1,U1,YD[,j])
                    res=MINQUEVarPre(Qa1,X1,U1,ML1,MR1,YD[,j],C)
                    JACVAR[j,]=res$v*VC
                    if(mk>0)JACPRE[j,]=res$Pre
                    JACB[j,]=res$b
               }
             }
             else{
                for(j in 1: TraitNum){
                    MR1=MINQUE_R(Qa1,U1,YD[,j])
                    res=MINQUEVarPre(Qa1,VI1,X1,U1,ML1,MR1,YD[,j],C)
                    JACVAR[j,]=res$v*VC
                    #if(mk>0)JACPRE[j,]=res$Pre
                    JACB[j,]=res$b
               }
             }

         }
         else if(Method==0){
            for(j in 1: TraitNum){
                res=REML(U1,X1,YD[,j],5)
                JACVAR[j,]=res[[1]]*VC
                if(mk>0)JACPRE[j,]=res[[3]]
                JACB[j,]=res[[2]]
            }

         }
         #ncol(JACVAR)
         JAC1[[i]]=JACVAR
         if(is.null(SIMUYES))if(mk>0)JAC2[[i]]=JACPRE
         #if(nx==1)JACB
         JAC3[[i]]=JACB
         #e1=e2+1
      } ## end of loop i
   
      #for(i in 1: JacNum){
      #   JV[,i]=JacNum*v0-df*JacVar[,i]
      #   JP[,i]=JacNum*pre0-df*JacPre[,i]
      #}
      for(s in 1:TraitNum){
          #s=1
          for(j in 1:JacNum){
             #j=1
             
             if(mk==0&&TraitNum==1)JacVar[j,]=JAC1[[j]]
             else JacVar[j,]=JAC1[[j]][s,]
             if(is.null(SIMUYES)){if(mk>0)JacPre[j,]=JAC2[[j]][s,]}
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
          if(is.null(SIMUYES)){ 
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
          }
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
   if(is.null(SIMUYES)){
   if(mk>0){
     ##Random effects   
     Pm=Pm/JacRep
     Pp=Pp/JacRep
     Cp=Cp/JacRep
     Cp1=Pm-t*Cp
     Cp2=Pm+t*Cp
   }
   }
   ##Fixed effects
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
   CNames=c("Estimate","SE","PValue","2.5%CL","97.5%CU")
   for(i in 1: TraitNum){
      VAR[[i]]=data.frame(Vm[i,],Cv[i,],Vp[i,],Cv1[i,],Cv2[i,])
      VC[[i]]=data.frame(VCm[i,],VCv[i,],VCp[i,],VCv1[i,],VCv2[i,])
      if(is.null(SIMUYES)){if(mk>0)PRE[[i]]=data.frame(Pm[i,],Cp[i,],Pp[i,],Cp1[i,],Cp2[i,])}
      B[[i]]=data.frame(Bm[i,],Cb[i,],Bp[i,],Cb1[i,],Cb2[i,])

      rownames(VAR[[i]])=VNames
      colnames(VAR[[i]])=CNames
      rownames(VC[[i]])=paste(VNames,"/VP",sep="")
      colnames(VC[[i]])=CNames
      if(is.null(SIMUYES)){
      if(mk>0){
         rownames(PRE[[i]])=Str
         colnames(PRE[[i]])=CNames
      }
      }
      colnames(B[[i]])=CNames
      rownames(B[[i]])=XNames

   }
   #VM=numeric(ml)
   #MSE=numeric(ml)
   #if(SimuID==1){
   #   for(i in 1:ml){
   #      VPower[i]=length(which(Vp[,i]<=0.05))/TraitNum
   #      VM[i]=mean(Vm[,i])
   #      MSE[i]=svar
   #   }
   #}
   names(VAR)=TraitNames
   names(VC)=TraitNames
   if(is.null(SIMUYES)){if(mk>0)names(PRE)=TraitNames}
   names(B)=TraitNames
   ##a=list(Var=VAR,PVar=VC,RandomEffect=PRE,FixedEffect=B,Method=Method,JacNum=JacNum,JacRep=JacRep,VNames=VNames,Y=Y,X=X,U=U)
   if(is.null(SIMUYES))a=list(Var=VAR,PVar=VC,RandomEffect=PRE,FixedEffect=B)##,Method=Method,JacNum=JacNum,JacRep=JacRep,VNames=VNames,Y=Y,X=X,U=U)
   else a=list(Var=VAR,PVar=VC,FixedEffect=B)
   #class(a)="jackknife.result"
   a$call=match.call()
   a$Var
   return(a)
} 



print.minque=function(mq,...){
   cat("Call: \n")
   print(mq$call)
   cat("\n")
   
   cat("\nMinque approach for variance components and for fixed and random effects\n")

   cat("\nVariance components:\n")
   print(mq$Var)
   cat("\n")

   cat("Fixed effects:\n")
   print(mq$FixEffect)
   cat("\n")
   
   cat("Random effects:\n")
   print(mq$RandomEffect)
   cat("\n")
}


print.minque.jack=function(mq,...){
   cat("Call: \n")
   print(mq$call)
   cat("\n")
   
   cat("\nJackknife technique with randomized grouping\n")

   cat("\nVariance components:\n")
   print(mq$Var)
   cat("\n")

   cat("Proportional variance components:\n")
   print(mq$PVar)
   cat("\n")
   
   cat("Fixed effects:\n")
   print(mq$FixedEffect)
   cat("\n")
   
   cat("Random effects:\n")
   print(mq$RandomEffect)
   cat("\n")
}

genmod.perm=function(gdata,au=NULL,PermNum=NULL,LUP=NULL,...){
    if(is.null(au))au=NULL
    if(is.null(PermNum))PermNum=200
    if(is.null(LUP))LUP=NULL
    #if(is.null(SIMUYES))result=genmod.oldperm(gdata,PermNum)
    #else
    result=genmod.oldperm(gdata,au,PermNum,LUP)
    return(result)

}

mq.perm=function(gdata,au=NULL,PermNum=NULL,LUP=NULL,...){
    if(is.null(au))au=NULL
    if(is.null(PermNum))PermNum=200
    if(is.null(LUP))LUP=NULL
    #if(is.null(SIMUYES))result=genmod.oldperm(gdata,PermNum)
    #else
    result=genmod.oldperm(gdata,au,PermNum,LUP)
    return(result)

}

genmod.oldperm=function(gdata,au=NULL,PermNum,LUP=NULL,...){
   if(is.null(LUP))LUP=NULL
   if(is.null(au))au=NULL
   result0=mq(gdata)
   tn=ncol(gdata$Y)
   n=nrow(gdata$Y)
   Y0=gdata$Y
   gc=length(gdata$VC)
   mc=nrow(result0$RandomEffect)
   fc=nrow(result0$FixedEffect)
   VP=matrix(0,gc,tn)
   VT=list()
   ET=list()
   FT=list()
   for(i in 1:tn){
      vp=numeric(gc)
      ep=numeric(mc)
      fp=numeric(fc)
      y=Y0[,i]
      Y=matrix(0,n,PermNum)
      for(j in 1:PermNum)Y[,j]=sample(y)
      colnames(Y)=paste("R",1:PermNum,sep="")
      gdata$Y=Y
      result=mq(gdata)
      V=result$Var
      v0=result0$Var[,i]
      e0=result0$RandomEffect[,i]
      f0=result0$FixedEffect[,i]
      PE=result$RandomEffect
      FE=result$FixedEffect
      for(j in 1:gc){
         id=which(V[j,]>v0[j])
         vp[j]=length(id)/PermNum
      }
      for(j in 1:mc){
         if(e0[j]<0){
             id=which(PE[j,]<e0[j])
             ep[j]=length(id)/PermNum
         }
         if(e0[j]>=0){
             id=which(PE[j,]>e0[j])
             ep[j]=length(id)/PermNum
         }
      }
      for(j in 1:fc){
         if(f0[j]<0){
             id=which(FE[j,]<f0[j])
             fp[j]=length(id)/PermNum
         }
         if(f0[j]>=0){
             id=which(FE[j,]>f0[j])
             fp[j]=length(id)/PermNum
         }
      }

      vp[gc]=1-vp[gc]
      VT[[i]]=cbind(v0,vp)
      ET[[i]]=cbind(e0,ep)
      FT[[i]]=cbind(f0,fp)
      colnames(VT[[i]])=c("Est","Pvalue")
      rownames(VT[[i]])=rownames(result0$Var)
      
      colnames(ET[[i]])=c("Pre","Pvalue")
      rownames(ET[[i]])=rownames(result0$RandomEffect)
      
      colnames(FT[[i]])=c("Est","Pvalue")
      rownames(FT[[i]])=rownames(result0$FixedEffect)

  }
  names(VT)=colnames(result0$Var)
  names(ET)=colnames(result0$Var)
  
  res=list(Var=VT,RandomEffect=ET,FixedEffect=FT)
  return(res)
}

reml.perm=function(gdata,au=NULL,PermNum,LUP=NULL,...){
   if(is.null(LUP))LUP=NULL
   if(is.null(au))au=NULL
   result0=reml(gdata)
   tn=ncol(gdata$Y)
   n=nrow(gdata$Y)
   Y0=gdata$Y
   gc=length(gdata$VC)
   mc=nrow(result0$RandomEffect)
   fc=nrow(result0$FixedEffect)
   VP=matrix(0,gc,tn)
   VT=list()
   ET=list()
   FT=list()
   for(i in 1:tn){
      vp=numeric(gc)
      ep=numeric(mc)
      fp=numeric(fc)
      y=Y0[,i]
      Y=matrix(0,n,PermNum)
      for(j in 1:PermNum)Y[,j]=sample(y)
      colnames(Y)=paste("R",1:PermNum,sep="")
      gdata$Y=Y
      result=reml(gdata)
      V=result$Var
      v0=result0$Var[,i]
      e0=result0$RandomEffect[,i]
      f0=result0$FixedEffect[,i]
      PE=result$RandomEffect
      FE=result$FixedEffect
      for(j in 1:gc){
         id=which(V[j,]>v0[j])
         vp[j]=length(id)/PermNum
      }
      for(j in 1:mc){
         if(e0[j]<0){
             id=which(PE[j,]<e0[j])
             ep[j]=length(id)/PermNum
         }
         if(e0[j]>=0){
             id=which(PE[j,]>e0[j])
             ep[j]=length(id)/PermNum
         }
      }
      for(j in 1:fc){
         if(f0[j]<0){
             id=which(FE[j,]<f0[j])
             fp[j]=length(id)/PermNum
         }
         if(f0[j]>=0){
             id=which(FE[j,]>f0[j])
             fp[j]=length(id)/PermNum
         }
      }

      vp[gc]=1-vp[gc]
      VT[[i]]=cbind(v0,vp)
      ET[[i]]=cbind(e0,ep)
      FT[[i]]=cbind(f0,fp)
      colnames(VT[[i]])=c("Est","Pvalue")
      rownames(VT[[i]])=rownames(result0$Var)
      
      colnames(ET[[i]])=c("Pre","Pvalue")
      rownames(ET[[i]])=rownames(result0$RandomEffect)
      
      colnames(FT[[i]])=c("Est","Pvalue")
      rownames(FT[[i]])=rownames(result0$FixedEffect)

  }
  names(VT)=colnames(result0$Var)
  names(ET)=colnames(result0$Var)
  
  res=list(Var=VT,RandomEffect=ET,FixedEffect=FT)
  return(res)
}

reml1.perm=function(gdata,au=NULL,PermNum=NULL,LUP=NULL,...){
   if(is.null(LUP))LUP=NULL
   if(is.null(au))au=NULL
   if(is.null(PermNum))PermNum=100
   result0=reml(gdata)
   tn=ncol(gdata$Y)
   n=nrow(gdata$Y)
   Y0=gdata$Y
   gc=length(gdata$VC)
   mc=nrow(result0$RandomEffect)
   VP=matrix(0,gc,tn)
   VT=list()
   ET=list()
   for(i in 1:tn){
      vp=numeric(gc)
      ep=numeric(mc)
      y=Y0[,i]
      Y=matrix(0,n,PermNum)
      for(j in 1:PermNum)Y[,j]=sample(y)
      colnames(Y)=paste("R",1:PermNum,sep="")
      gdata$Y=Y
      result=reml(gdata)
      V=result$Var
      v0=result0$Var[,i]
      e0=result0$RandomEffect[,i]
      PE=result$RandomEffect
      for(j in 1:gc){
         id=which(V[j,]>v0[j])
         vp[j]=length(id)/PermNum
      }
      for(j in 1:mc){
         if(e0[j]<0){
             id=which(PE[j,]<e0[j])
             ep[j]=length(id)/PermNum
         }
         if(e0[j]>=0){
             id=which(PE[j,]>e0[j])
             ep[j]=length(id)/PermNum
         }
      }
      vp[gc]=1-vp[gc]
      VT[[i]]=cbind(v0,vp)
      ET[[i]]=cbind(e0,ep)
      colnames(VT[[i]])=c("Est","Pvalue")
      rownames(VT[[i]])=rownames(result0$Var)
      
      colnames(ET[[i]])=c("Pre","Pvalue")
      rownames(ET[[i]])=rownames(result0$RandomEffect)
  }
  names(VT)=colnames(result0$Var)
  names(ET)=colnames(result0$Var)
  
  res=list(Var=VT,Pre=ET)
  return(res)
}
#SIMUYES=1
genmod.jack=function(gdata,JacNum=NULL,JacRep=NULL,ALPHA=NULL,au=NULL,LUP=NULL,...){
  if(is.null(JacNum))JacNum=10
  if(is.null(JacRep))JacRep=1
  if(is.null(ALPHA))ALPHA=0.05
  if(is.null(au))au=NULL
  if(is.null(LUP))LUP=NULL
  #if(is.null(SIMUYES))result=genmod.jack1(gdata=gdata,JacNum=JacNum,JacRep=JacRep,au=au,ALPHA=ALPHA)
  result=genmod.jack1(gdata,JacNum,JacRep,ALPHA,au,LUP)
}

genmod0=function(gdata,au=NULL,LUP=NULL){
   gdata0=gdata
   if(is.null(au))au=NULL
   if(is.null(LUP))LUP=NULL
   #return(list(au=au,SIMUYES=SIMUYES))

   res0=genmod(gdata0,au,LUP)
}
#JacNum=10
#JacRep=1
genmod.jack1=function(gdata,JacNum,JacRep,ALPHA,au=NULL,LUP=NULL,...){
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

   res0=genmod0(gdata0,au=au,LUP=LUP)
   V0=res0$Var
   B0=res0$FixEffect
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
         res=genmod(gdata,au=au,LUP=LUP)
         JACVAR=t(res$Var)
         JACB=t(res$FixEffect)
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
   a=list(Var=VAR,PVar=VC,RandomEffect=PRE,FixedEffect=B,Prediction=Prediction)

   #class(a)="jackknife.result"
   #a$call=match.call()
   #a$Var
   return(a)
}

crecode1=function(v){
  cname=unique(v)
  gn=length(cname)
  n=length(v)
  v1=numeric(n)
  for(i in 1:gn){
    id=which(v==cname[i])
    v1[id]=i
  }
  result=list(nv=v1,vname=cname)
  return(result)
}

comid=function(name0,name1){
   id1=NULL
   cn=length(name1)
   for(i in 1:cn)id1=c(id1,which(name0==name1[i]))
   return(id1)
}
