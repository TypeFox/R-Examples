C_piBTL<-function(Matpair,Constraint=0,eps1=0.0001,Pi=NULL,TestPi=FALSE,Zht=NULL)
{

  if(is(Matpair,"DataPairComp")==TRUE)
  {
    nsujet<-length(Matpair@Cons)
    nprod<-length(Matpair@Prod)
    Matcla<-vector("list")
    Mat<-vector("list")
    if (is.null(Zht)==TRUE)
    {
       Zht<-matrix(1,nsujet,1)
       t<-1
       k<-1
    }
    Mat[[k]]<-apply(simplify2array(Matpair@Paircomp[[k]]),c(2,1),somme,Zht[,t])
    Matcla[[t]]<-Mat
    Matpair<-matrix(simplify2array(Matcla[[t]][k]),nrow=nprod,ncol=nprod,byrow=TRUE)
  }
  
  
  nprod<-ncol(Matpair)
  if (is.null(Pi))
  {
    Pi<-rowSums(Matpair)
    if (Constraint==0)
    {
      Pi<-(t(t(Pi))+0.000000001)/sum(Pi)
    }    else
    {
      c<-(prod(Pi))^(1/nprod)
      Pi<-t(t(Pi))/c
    }
  }
  lnL<--100000000
  diff<-100000
  maxiteru<-50
  iteru<-1
  while ((iteru<maxiteru)&(diff>eps1))
  {
    lnLold<-lnL
    Piold<-Pi
    som1<-rowSums(Matpair)
    vpi<-rbind(t(Piold),matrix(Piold,nrow=nprod,ncol=nprod))
    u<-diag(1,nprod)
    u<-cbind(rep(1,nprod),u)

    denom<-u%*%vpi   
    w<-t(Matpair)+Matpair
    som2<-(w-diag(diag(w)))/denom

    som2<-colSums(som2)   
    Pi<-som1/som2
    if (Constraint==0)
    {
      Pi<-t(t(Pi+0.00000001))/sum(Pi)
    }   else
    {
      c<-(prod(Pi))^(1/nprod)
      Pi<-t(t(Pi))/c
    }
    lnL<-0
    for (i in 1:(nprod-1))
    {
      for (j in (i+1):nprod)
      {
        lnL<-lnL+Matpair[i,j]*log((Pi[i])/(Pi[j]))+(Matpair[i,j]+Matpair[j,i])*log((Pi[j])/((Pi[i])+(Pi[j])))
      }
    }
    diff<-lnL-lnLold
    iteru<-iteru+1
  }
  
  
  
  if (TestPi==TRUE)
  {
    lvrH0<-0
    lvrH1<-lnL
    for (i in 1:(nprod-1))
    {
      for (j in (i+1):nprod)
      {
        lvrH0<-lvrH0+(Matpair[i,j]+Matpair[j,i])*log(0.5)
      }
    }
    
    diff<-2*(lvrH1-lvrH0)
    PValue<-(1-pchisq(diff,(nprod-1)))
    H1<-(PValue<0.05)
    
    
    Restest<-vector("list")
    
    
    Fish<-matrix(0,nrow=nprod,ncol=nprod)
    
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
              termeii<-termeii+(Matpair[i,l]+Matpair[l,i])*Pi[l]/((Pi[i]+Pi[l])^2)
            }
          }
          Fish[i,i]<-termeii/Pi[i]
        }
        else
        {
          Fish[i,j]<--(Matpair[i,j]+Matpair[j,i])/((Pi[i]+Pi[j])^2)
          Fish[j,i]<-Fish[i,j]
        }
      }
    }
    
    Varcov<-ginv(Fish)
    Sigma2<-diag(Varcov)
    restestij<-NULL
    for (i in 1:(nprod-1))
    {
      for (j in (i+1):nprod)
      {
        testij<-(Pi[i]-Pi[j])/(sqrt(Sigma2[i]+Sigma2[j]-2*Varcov[i,j]))
        PValueij<-2*(pnorm(-abs(testij)))
        H1ij<-(PValueij<0.05)
        restestij<-rbind(restestij,c(i,j,testij,PValueij,H1ij))
      }
    }
    colnames(restestij)<-c("i","j","test i=j","PValue","H1")
    res<-list(Pi=Pi,lnL=lnL,lvrH0=lvrH0,lvrH1=lvrH1,lRatio=diff,PValue=PValue,H1=H1,VarcovPi=Varcov,restestij=restestij)
  }
  else
  {
    res<-list(Pi=Pi,lnL=lnL)
  }
  
  return(res)
}
