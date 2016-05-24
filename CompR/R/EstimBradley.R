EstimBradley<-function(Data,Constraint=0,Tcla=1,eps=0.0001,eps1=0.0001,TestPi=TRUE)
{


  
  nsujet<-length(Data@Cons)
  ncrit<-length(Data@Crit)
  Appart<-ClassifPaired(Data,Tcla)
  pcla<-matrix(rep(0,nsujet*Tcla), ncol=Tcla)
  namePi<-matrix(0,nrow=2,ncol=Tcla)
  colnames(namePi)<-colnames(namePi,do.NULL = FALSE, prefix = "Class")
  for (t in 1:Tcla)
  {
    pcla[which(Appart==t),t]<-1
  }
  Zht<-matrix(0,nrow=nsujet,ncol=Tcla)
  for (t in 1:Tcla)
  {
    Zht[,t]<-pcla[,t]
  }
  lambda=rep(0,Tcla)
  
  for (t in 1:Tcla)
  {
    lambda[t]<-length(Data@Cons[Appart==t])/length(Data@Cons)
  }
  
  Mat<-vector("list")
  Matcla<-vector("list")
  for (t in 1:Tcla)
  {
    for (k in 1:ncrit)
    {
      Mat[[k]]<-apply(simplify2array(Data@Paircomp[[k]]),c(2,1),somme,Zht[,t])
    }
    Matcla[[t]]<-Mat
  }
  Matclaold<-Matcla
  Pinouv<-vector("list")
  lnLnouv<-vector("list")
  
  for (t in 1:Tcla)
  {
    Pinouvt<-NULL
    lnLnouvt<-NULL
    for (k in 1:ncrit)
    {
      nprod<-length(Data@Prod)
      Matpair<-matrix(simplify2array(Matcla[[t]][k]),nrow=nprod,ncol=nprod,byrow=TRUE)
      Estim<-C_piBTL(Matpair,Constraint=Constraint,eps1=eps1,Zht=Zht)
      Pinouvt<-cbind(Pinouvt,Estim$Pi)
      lnLnouvt<-cbind(lnLnouvt,Estim$lnL)
      rownames(Pinouvt)<-Data@Prod
    }
    Pinouv[[t]]<-Pinouvt
    lnLnouv[[t]]<-lnLnouvt
  }
  Piold<-Pinouv
  lnLold<-lnLnouv
  Zhtold<-Zht
  
  lambdanouv<-estimlambda(Piold,Data,lambda)
  maxiter<-1000
  iter<-1
  Lvrold<--100000000
  Lvrnouv<-lambdanouv$lvrnouv
  difLvr<-Lvrnouv-Lvrold
  Lvr<-t(c(iter,Lvrold,Lvrnouv,difLvr))
  
  
  while ((difLvr>eps)&(iter<maxiter))
    
  {
    Zht<-lambdanouv$Zhtnouv
    Piold<-Pinouv
    lambda<-lambdanouv$lambdanouv
    Lvrold<-lambdanouv$lvrnouv
    Mat<-vector("list")
    Matcla<-vector("list")
    Pinouv<-vector("list")
    lnLnouv<-vector("list")
    
    Mat<-vector("list")
    Matcla<-vector("list")
    for (t in 1:Tcla)
    {
      for (k in 1:ncrit)
      {
        Mat[[k]]<-apply(simplify2array(Data@Paircomp[[k]]),c(2,1),somme,Zht[,t])
      }
      Matcla[[t]]<-Mat
    }
    
    for (t in 1:Tcla)
    {
      Pinouvt<-NULL
      lnLnouvt<-NULL
      for (k in 1:ncrit)
      {
        
        Matpair<-matrix(simplify2array(Matcla[[t]][k]),nrow=nprod,ncol=nprod,byrow=TRUE)
        Estim<-C_piBTL(Matpair,Pi=Piold[[t]][,k],Constraint=Constraint,eps1=eps1,Zht=Zht)
        Pinouvt<-cbind(Pinouvt,Estim$Pi)
        lnLnouvt<-cbind(lnLnouvt,Estim$lnL)
      }
      rownames(Pinouvt)<-rownames(Pinouvt,do.NULL=FALSE,prefix="P")
      Pinouv[[t]]<-Pinouvt
      colnames(Pinouv[[t]])<-Data@Crit
      lnLnouv[[t]]<-lnLnouvt
    }
    lambdanouv<-estimlambda(Pinouv,Data,lambda)
    Lvrnouv<-lambdanouv$lvrnouv
    difLvr<-Lvrnouv-Lvrold
    
    iter<-iter+1
    Lvr<-rbind(Lvr,t(c(iter,Lvrold,Lvrnouv,difLvr)))
  }
  
  Zhtclass<-apply(lambdanouv$Zhtnouv,1,which.max)
  
  
  ddl<-ncrit*(Tcla*(nprod-1))+Tcla-1
  aic<--2*max(Lvr[,2])+2*ddl
  n<-0
  for (k in 1:ncrit)
  {
    n<-n+sum(apply(simplify2array(Data@Paircomp[[k]]),c(1,2),sum))
  }
  bic<--2*max(Lvr[,2])+ddl*log(n)
  caic<--2*max(Lvr[,2])+ddl*(log(n)+1)
  prefix<-rep("Class",Tcla)
  Class<-paste(prefix,c(1:Tcla),sep="")
  Zht<-lambdanouv$Zhtnouv
  colnames(Zht)<-Class
  Zh<-cbind(Zht,Zhtclass)
  Ic<-cbind(Tcla,aic,bic,caic)
  Lvriter<-Lvr
  colnames(Lvriter)<-c("iter","LvrOld","LvrNew","DiffLvr")
  Lvr<-lambdanouv$lvrnouv
  Pi<-Piold
  names(Pi)<-Class
  Lambda<-t(lambdanouv$lambdanouv)
  colnames(Lambda)<-Class
  if (TestPi==FALSE)
  {
    new(Class="BradleyEstim",Lvriter=Lvriter,Lvr=Lvr,Lambda=Lambda,Pi=Pi,Zh=Zh,Ic=Ic)
  }
  else
  {
    
    
    lvrH0<-matrix(0,nrow=Tcla,ncol=ncrit)
    lvrH1<-matrix(0,nrow=Tcla,ncol=ncrit)
    
    for (t in 1:Tcla)
    {
      for (k in 1:ncrit)
      {
        Pitk<-Piold[[t]][,k]
        for (h in 1:nsujet)
        {
          Mcomp<-Data@Paircomp[[k]][[h]]
          for (i in 1:(nprod-1))
          {
            for (j in (i+1):nprod)
            {
              lvrH0[t,k]<-lvrH0[t,k]+lambdanouv$Zhtnouv[h,t]*log(0.5)*(Mcomp[i,j]+Mcomp[j,i])
              
              lvrH1[t,k]<-lvrH1[t,k]+lambdanouv$Zhtnouv[h,t]*(log(Pitk[i]/(Pitk[i]+Pitk[j]))*Mcomp[i,j]+log(Pitk[j]/(Pitk[i]+Pitk[j]))*Mcomp[j,i])
            }
          }                                              
        }
      }
    }
    colnames(lvrH0)<-Data@Crit
    rownames(lvrH0)<-rownames(lvrH0,do.NULL = FALSE, prefix = "Class")
    colnames(lvrH1)<-Data@Crit
    rownames(lvrH1)<-rownames(lvrH1,do.NULL = FALSE, prefix = "Class")
    diff<-2*(lvrH1-lvrH0)
    PValue<-(1-pchisq(diff,(nprod-1)))
    colnames(PValue)<-Data@Crit
    rownames(PValue)<-rownames(PValue,do.NULL = FALSE, prefix = "Class")
    H1<-(PValue<0.05)
    restestglob<-list(lvrH0=lvrH0,lvrH1=lvrH1,lRatio=diff,PValue=PValue,H1=H1)
    
    Varcov<-vector("list")
    Restest<-vector("list")
    for (t in 1:Tcla)
    {
      Varcovt<-vector("list")
      Restestk<-vector("list")
      for (k in 1:ncrit)
      {
          Pitk<-Piold[[t]][,k]
          Fishtk<-matrix(0,nrow=nprod,ncol=nprod)
          Matsegtk<-matrix(0,nrow=nprod,ncol=nprod)
          for (h in 1:nsujet)
          {
            Matsegtk<-Matsegtk+lambdanouv$Zhtnouv[h,t]*Data@Paircomp[[k]][[h]]
          }
          for (i in 1:nprod)
          {
            for (j in 1:nprod)
            {
              if (j==i)
              {
                termeii<-0
                for (l in 1:nprod)
                {
                  if (l!=i)
                  {
                    termeii<-termeii+(Matsegtk[i,l]+Matsegtk[l,i])*Pitk[l]/((Pitk[i]+Pitk[l])^2)
                  }
                }
                Fishtk[i,i]<-termeii/Pitk[i]
              }
              else
              {
                Fishtk[i,j]<--(Matsegtk[i,j]+Matsegtk[j,i])/((Pitk[i]+Pitk[j])^2)
                Fishtk[j,i]<-Fishtk[i,j]
              }
            }
          }
          Varcovtk<-ginv(Fishtk)
          Varcovt[[k]]<-Varcovtk
          Sigma2tk<-diag(Varcovtk)
          restesttkij<-NULL
          for (i in 1:(nprod-1))
          {
            for (j in (i+1):nprod)
            {
              testtkij<-(Pitk[i]-Pitk[j])/(sqrt(Sigma2tk[i]+Sigma2tk[j]-2*Varcovtk[i,j]))
              PValuetkij<-2*(pnorm(-abs(testtkij)))
              H1ij<-(PValuetkij<0.05)
              restesttkij<-rbind(restesttkij,c(t,k,i,j,testtkij,PValuetkij,H1ij))
            }
          }
          rownames(PValue)<-rownames(PValue,do.NULL = FALSE, prefix = "Class")       
          colnames(restesttkij)<-c("Class","Crit","i","j","test i=j","PValue","H1")
          Restestk[[k]]<-restesttkij
      }
      names(Varcovt)<-Data@Crit
      Varcov[[t]]<-Varcovt
      names(Restestk)<-Data@Crit
      Restest[[t]]<-Restestk
    }
    names(Varcov)<-paste("Class",c(1:Tcla),sep="")
    names(Restest)<-rownames(lvrH0)

    
    new(Class="BradleyEstim",Lvriter=Lvriter,Lvr=Lvr,Lambda=Lambda,Pi=Pi,Zh=Zh,Ic=Ic,Restestglob=restestglob,Restestprod=Restest,Varcov=Varcov)
  }

  
}