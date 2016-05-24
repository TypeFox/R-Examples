#######################################################################
# this code is for calculating the asymptotic covariance 
# matrix of CL1 estimates
###########################################################################

bvn<-function(lb, ub, rh)
{ if(rh>0)
  { exchmvn(lb, ub, rh) }
  else
  { rhmat=matrix(c(1,rh,rh,1),2,2)
    pmvnorm(lb,ub,c(0,0),rhmat)[1]
  }
}

bvn.deriv.margin<-function(lb, ub, rh, k, ksign)
{ if(rh>0)
  { exchmvn.deriv.margin(lb, ub, rh, k, ksign) }
  else
  { rhmat=matrix(c(1,rh,rh,1),2,2)
    mvn.deriv.margin(lb,ub,c(0,0),rhmat,k,ksign)$deriv
  }
}

bvn.deriv.rho<-function(lb, ub, rh)
{ if(rh>0)
  { exchmvn.deriv.rho(lb, ub, rh) }
  else
  { rhmat=matrix(c(1,rh,rh,1),2,2)
    mu<-c(0,0)
    dmvnorm(ub,mu,rhmat)-dmvnorm(c(ub[1],lb[2]),mu,
    rhmat)-dmvnorm(c(lb[1],ub[2]),mu,rhmat)+dmvnorm(lb,mu,rhmat)
  }
}



# derivative of the univariate cdf wrt gam
# input:
# y: the response cector
# nu: the nu's
# margmodel: see weightesScores
# link: see weightesScores
# output: the derivative
der.cdf.gam<-function(y,mu,gam,u,link)
{ if(y<1) return(0)
  k<-1:y
  cuts<-c(-10,gam,10)
  lb=cuts[k]+mu
  ub=cuts[k+1]+mu
  if(link=="probit") { dlatent=dnorm } else { dlatent=dlogis }
  dlatentub<-dlatent(ub)
  dlatentlb<-dlatent(lb)
  res<-rep(NA,y)
  for(i in 1:y)
  { if(u==i)
    { res[i]=dlatentub[i] } 
    else
    { if(u==i-1)
      { res[i]=-dlatentlb[i] } 
    else {res[i]=0}}
  }
  sum(res)
}




# derivative with respect to ga_1,ga_2,...,ga_d of the BCL likelhood
# input:
# y: the d-variate vector response vector
# ...: see weightedScores
# output: a d-variate vector with the derivatives
cl1.der.gam<-function(zlow,zupp,dzlow,dzupp,y,mu,gam,rho,
bivpairs,corstr,d,d2,u,link)
{ s<-rep(0,d)
  for(i in 1:d)
  { for(j in 1:d2)
  { k1<-bivpairs[j,][1]
    k2<-bivpairs[j,][2]
    if(corstr=="exch") rh<-rho else rh<-rho^(k2-k1)
    if(i==k1)
    { temp1<-bvn.deriv.margin(zlow[c(k1,k2)],zupp[c(k1,k2)],
      rh,1,1)*der.cdf.gam(y[k1],mu[k1],gam,u,link)/dzupp[k1]
      temp2<-bvn.deriv.margin(zlow[c(k1,k2)],zupp[c(k1,k2)],
      rh,1,-1)*der.cdf.gam(y[k1]-1,mu[k1],gam,u,link)/dzlow[k1]
      prob<-bvn(zlow[c(k1,k2)],zupp[c(k1,k2)],rh)
      derprob<-(temp1+temp2)/prob
    } else {
    if(i==k2)
    { temp1<-bvn.deriv.margin(zlow[c(k1,k2)],zupp[c(k1,k2)],
      rh,2,1)*der.cdf.gam(y[k2],mu[k2],gam,u,link)/dzupp[k2]
      temp2<-bvn.deriv.margin(zlow[c(k1,k2)],zupp[c(k1,k2)],
      rh,2,-1)*der.cdf.gam(y[k2]-1,mu[k2],gam,u,link)/dzlow[k2]
      prob<-bvn(zlow[c(k1,k2)],zupp[c(k1,k2)],rh)
      derprob<-(temp1+temp2)/prob
    } else {derprob<-0}}
    s[i]<-s[i]+derprob
  }}
  s
}

# derivative with respect to rho of the BCL
# y: the d-variate vector response vector
# ...: see weightedScores
# output:  the derivative
cl1.der.rho<-function(zlow,zupp,rho,bivpairs,corstr,d2)
{ s<-0
  for(j in 1:d2)
  { k1<-bivpairs[j,][1]
    k2<-bivpairs[j,][2]
    if(corstr=="exch")
    { der<-bvn.deriv.rho(zlow[bivpairs[j,]],zupp[bivpairs[j,]],rho)
      prob<-bvn(zlow[bivpairs[j,]],zupp[bivpairs[j,]],rho)
    } else {
    if(corstr=="ar")
    { t1<-k2-k1
      rhar<-rho^(k2-k1)
      der<-bvn.deriv.rho(zlow[bivpairs[j,]],zupp[bivpairs[j,]],
      rhar)*t1*rhar/rho
      prob<-bvn(zlow[bivpairs[j,]],zupp[bivpairs[j,]],rhar)
    }
    }
    s<-s+ der/prob
  }
  s
}


# d=dimension, K=#categories, ii = decimal
# return vector of size d, each element in 1..K
d2v=function(d, K, ii)
{ t=ii-1
  jj=rep(0,d)
  for(i in seq(d,1,-1))
  { jj[i]=t%%K; t=floor(t/K); }
  jj
}

# d-variate rectangle probability
# y: the d-variate vector response vector
# ...: see weightedScores
# output:  the bivariate rectangle probability
mrect.prob<-function(zlow,zupp,rhomat,bivpairs,corstr,mvncmp)
{ if(corstr=="exch" & min(rhomat)>0 )
  { prob<-exchmvn(zlow,zupp,rhomat[1,2])
  }
  else
  { d<-length(zlow)
    if(mvncmp==1)
    { prob<-mvnapp(zlow,zupp,rep(0,d),sigma=rhomat)$pr } else {
    set.seed(12345)
    prob<-pmvnorm(lower=zlow,upper=zupp,mean=rep(0,d),corr=rhomat)[1]
    }
  }
  prob
}




# derivative with respect to rho of the bivariate rectangle probability
rect.der.rho<-function(zlow,zupp,rho,corstr,j1,j2)
{ if(corstr=="ar")
  { tem1=j2-j1
    tem2=bvn.deriv.rho(zlow,zupp,rho)
    tem3=(rho)^(1/tem1)
    derprob=tem1*tem2*tem3^(tem1-1)
  } else { derprob=bvn.deriv.rho(zlow,zupp,rho) }
  derprob
 }


# derivative with respect to gam_j of the bivariate recatngle
rect.der.gam<-function(zlow,zupp,dzlow,dzupp,y,mu,gam,rho,j,link)
{ ub<-noCategories(gam)
  derprob<-rep(NA,ub-1)
  der1=bvn.deriv.margin(zlow,zupp,rho,j,1)
  der2=bvn.deriv.margin(zlow,zupp,rho,j,-1)
  for(u in 1:(ub-1))
  { temp1<-der1*der.cdf.gam(y[j],mu[j],gam,u,link)/dzupp[j]
    temp2<-der2*der.cdf.gam(y[j]-1,mu[j],gam,u,link)/dzlow[j]
    derprob[u]<-temp1+temp2
  }
  derprob
}




fisher.gam.rho<-function(mu,gam,rho,j,ub,corstr,j1,j2,link)
{ s<-rep(0,ub-1)
  nvect=ub^2
  for(ii in 1:nvect)
  { y=d2v(2,ub,ii)+1
    vlow<-pmargmodel.ord(y-1,mu,gam,link)
    tem<-dmargmodel.ord(y,mu,gam,link)
    vupp<-vlow+tem
    zlow=qnorm(vlow)
    zupp=qnorm(vupp)
    zlow[zlow < -10]<--10
    zupp[zupp > 10]<-10
    dzlow=dnorm(zlow)
    dzupp=dnorm(zupp)
    der1<-rect.der.gam(zlow,zupp,dzlow,dzupp,y,mu,gam,rho,j,link)
    der2<-rect.der.rho(zlow,zupp,rho,corstr,j1,j2)
    prob<-bvn(zlow,zupp,rho)
    s<-s+der1*der2/prob
  }
  s
}

# the inverse Godambe matrix for CL1 estimates
# input:
# ...: see weightedScores
clic1dePar<-function(nbcl,r,b,gam,xdat,id,tvec,corstr,WtScMat,link,mvncmp)
{ if(is.matrix(xdat))
  { dim<-dim(xdat)
    n<-1:dim[1]
    p<-dim[2]
  } else {n<-1:length(xdat); p<-1}
  
  uid<-unique(id)
  d<-id.size(id)
  maxd<-max(d)
  q<-length(gam)
    
  tvec<-id.time(tvec,d)
  bivpairs<-maxpairs(d)
  d2<-nrow(bivpairs)
  rmat<-cormat(maxd,r,bivpairs,corstr)
  t=p+q+1
  Dmat<-M<-matrix(0,t,t)
  m<-0
  for(i in uid)
  { m<-m+1
    cases<-id==i
    irow=n[cases]
    ti<-tvec[irow]
    if(length(ti)>1)
    {
    newmui<-rep(NA,maxd)
    if(is.matrix(xdat)) 
    { x<-xdat[cases,]
      newx<-matrix(NA,maxd,p) 
      newx[ti,]<-x}  else {
      x<-xdat[cases]
      newx<-rep(NA,maxd) 
      newx[ti]<-x 
    }
    
    
    mui<-ordreg.mu(x,b)
    newmui[ti]<-mui
    ub<-q+1
    nvect=ub^d[m] # change here nvect=ub^maxd
    nvect2=ub^2
    
    y=vlow=tem=matrix(NA,nvect,maxd)
    for(ii in 1:nvect)
    { y[ii,]=d2v(maxd,ub,ii)+1
      vlow[ii,]<-pmargmodel.ord(y[ii,]-1,newmui,gam,link)
      tem[ii,]<-dmargmodel.ord(y[ii,],newmui,gam,link)
    }
    vupp<-vlow+tem
    zlow=qnorm(vlow)
    zupp=qnorm(vupp)
    zlow[zlow < -10]<--10
    zupp[zupp > 10]<-10
    dzlow=dnorm(zlow)
    dzupp=dnorm(zupp) 
    D22i<-0
    Delta21i<-matrix(NA,d2,maxd*q)
    for(i1 in 1:d2)
    { s<-rep(0,q)
      k1<-bivpairs[i1,1]
      k2<-bivpairs[i1,2]
      if((sum(k1==ti)==1) & (sum(k2==ti)==1))
      {
      for(ii in 1:nvect2)
      { y2=d2v(2,ub,ii)+1
        vlow2<-pmargmodel.ord(y2-1,newmui[c(k1,k2)],gam,link)
        tem2<-dmargmodel.ord(y2,newmui[c(k1,k2)],gam,link)
        vupp2<-vlow2+tem2
        zlow2=qnorm(vlow2)
        zupp2=qnorm(vupp2)
        zlow2[zlow2 < -10]<--10
        zupp2[zupp2 > 10]<-10
        dzlow2=dnorm(zlow2)
        dzupp2=dnorm(zupp2) 
        der2<-rect.der.rho(zlow2,zupp2,rmat[k1,k2],corstr,k1,k2)
        prob<-bvn(zlow2,zupp2,rmat[k1,k2])
        D22i<-D22i+der2*der2/prob
        der1<-NULL
        for(i2 in 1:maxd)
        { if(i2==k1)
          { der1<-c(der1,rect.der.gam(zlow2,
            zupp2,dzlow2,dzupp2,
            y2,newmui[c(k1,k2)],gam,rmat[k1,k2],1,link))
          } else {
          if(i2==k2)
          { der1<-c(der1,rect.der.gam(zlow2,
            zupp2,dzlow2,dzupp2,
            y2,newmui[c(k1,k2)],gam,rmat[k1,k2],2,link)) 
          } else {
          der1<-c(der1,rep(0,q))}} # new change instead of c(der1,0)
        }
        s<-s+der1*der2/prob
        }
      }
        Delta21i[i1,]<-s
    }
    Delta21i<-apply(Delta21i,2,sum)
    
    if(d[m]<maxd)
    { # new code
    und<-d[m]
    und2<-choose(und,2)
    unpairs<-maxpairs(und)
    unrmat<-rmat[ti,ti]
    
    Omega22i<-0
    Omega12i<-rep(0,q*und)
    unnvect=ub^und
    for(ii in 1:unnvect)
    { y=d2v(und,ub,ii)+1
      vlow<-pmargmodel.ord(y-1,mui,gam,link)
      tem<-dmargmodel.ord(y,mui,gam,link)
      vupp<-vlow+tem
      zlow=as.vector(qnorm(vlow))
      zupp=as.vector(qnorm(vupp))
      zlow[zlow < -10]<--10
      zupp[zupp > 10]<-10
      dzlow=dnorm(zlow)
      dzupp=dnorm(zupp) 
      
      der2<-cl1.der.rho(zlow,zupp,r,unpairs,corstr,und2)
      prob<-mrect.prob(zlow,zupp,unrmat,unpairs,corstr,mvncmp)
      Omega22i<-Omega22i+der2*der2*prob
      mder<-rep(NA,und)
      dergam<-NULL
      for(jj in 1:und)
      { for(k in 1:q) #q=ub-1
        { dergam<-c(dergam,iderlik.gam.ord(mui[jj],gam,y[jj],k,link))}  
      }
      Omega12i<-Omega12i+dergam*der2*prob
    }
    } else {
      Omega22i<-0
      Omega12i<-rep(0,q*maxd)
      for(ii in 1:nvect)
      { der2<-cl1.der.rho(zlow[ii,],zupp[ii,],r,bivpairs,corstr,d2)
        prob<-mrect.prob(zlow[ii,],zupp[ii,],rmat,bivpairs,corstr,mvncmp)
        Omega22i<-Omega22i+der2*der2*prob
        mder<-rep(NA,maxd)
        dergam<-NULL
        for(jj in 1:maxd)
        { for(k in 1:q) #q=ub-1
        { dergam<-c(dergam,iderlik.gam.ord(newmui[jj],gam,y[ii,jj],k,link))}  
        }
        Omega12i<-Omega12i+dergam*der2*prob
      }
    }
      
      
    seli<-subselect.ord(ti,q)
    Xi<-WtScMat$X[,,m]
    Xi<-Xi[,seli]
    tXi<-t(Xi)
    Delta21i=Delta21i[seli]
    D21i<-Delta21i%*%tXi
    Delta11<-WtScMat$delta[,,m] # from the wtsc code: weightMat()
    Delta11<-Delta11[seli,seli]
    D11i<-Xi%*%Delta11%*%tXi
    d0<-dim(t(D21i))
    D12i<-matrix(0,d0[1],d0[2])
    Di<-rbind(cbind(D11i,D12i),cbind(D21i,D22i))
    Omega12i<-as.matrix(Omega12i)
    Omega11i<-WtScMat$omega[,,m] # from the wtsc code: weightMat() function
    Omega11i<-Omega11i[seli,seli]
    M11i<-Xi%*%Omega11i%*%tXi
    M12i<-Xi%*%Omega12i
    M21i<-t(M12i)
    M22i<-Omega22i
    Mi<-rbind(cbind(M11i,M12i),cbind(M21i,M22i))
    # summations
    Dmat<-Dmat+Di
    M<-M+Mi
 }}
 #inDmat<-solve(Dmat)
 #tr<-sum(diag(M%*%inDmat))
 tr<-sum(diag(solve(Dmat,M)))
 list(AIC=2*(nbcl+tr),BIC=2*nbcl+log(m)*tr)
}


cov3<-function(j1,j2,j3,mu,gam,rhomat,bivpairs,corstr,ub,link,mvncmp)
{ s<-0
  q=ub-1
  if(j1==j2)
  { nvect=ub^2
    for(ii in 1:nvect)
    { y=d2v(2,ub,ii)+1
      vlow<-pmargmodel.ord(y-1,mu[c(j2,j3)],gam,link)
      tem<-dmargmodel.ord(y,mu[c(j2,j3)],gam,link)
      vupp<-vlow+tem
      zlow=qnorm(vlow)
      zupp=qnorm(vupp)
      zlow[zlow < -10]<--10
      zupp[zupp > 10]<-10
      #scj1<-iderlik.nu(mu[j1],gam,invgam,y[1],margmodel,link)
      scj1<-NULL
      for(k in 1:q) #q=ub-1
      { scj1<-c(scj1,iderlik.gam.ord(mu[j1],gam,y[1],k,link)) }  
      scj2j3<-rect.der.rho(zlow,zupp,rhomat[j2,j3],corstr,j2,j3)
      #prob<-rect.prob(y,x[c(j2,j3),],b,gam,invgam,rhomat[j2,j3],marmodel,link)
      s<-s+scj1*scj2j3#/prob
    }
  } else {
    if(j1==j3)
    { nvect=ub^2
      for(ii in 1:nvect)
      { y=d2v(2,ub,ii)+1
        vlow<-pmargmodel.ord(y-1,mu[c(j3,j2)],gam,link)
        tem<-dmargmodel.ord(y,mu[c(j3,j2)],gam,link)
        vupp<-vlow+tem
        zlow=qnorm(vlow)
        zupp=qnorm(vupp)
        zlow[zlow < -10]<--10
        zupp[zupp > 10]<-10
        #scj1<-iderlik.nu(mu[j1],gam,invgam,y[1],margmodel,link)
        scj1<-NULL
        for(k in 1:q) #q=ub-1
        { scj1<-c(scj1,iderlik.gam.ord(mu[j1],gam,y[1],k,link)) }
        scj2j3<-rect.der.rho(zlow[c(2,1)],zupp[c(2,1)],rhomat[j3,j2],corstr,j3,j2)
        #prob<-rect.prob(y[c(2,1)],x[c(j2,j3),],b,gam,invgam,rhomat[j2,j3],marmodel,link)
        s<-s+scj1*scj2j3#/prob
      }
    }
    else
    { nvect=ub^3
      for(ii in 1:nvect)
      { y=d2v(3,ub,ii)+1
        vlow<-pmargmodel.ord(y-1,mu[c(j1,j2,j3)],gam,link)
        tem<-dmargmodel.ord(y,mu[c(j1,j2,j3)],gam,link)
        vupp<-vlow+tem
        zlow=qnorm(vlow)
        zupp=qnorm(vupp)
        zlow[zlow < -10]<--10
        zupp[zupp > 10]<-10
        #scj1<-iderlik.nu(mu[j1],gam,invgam,y[1],margmodel,link)  #not sure about y[1]
        scj1<-NULL
        for(k in 1:q) #q=ub-1
        { scj1<-c(scj1,iderlik.gam.ord(mu[j1],gam,y[1],k,link)) }
        scj2j3<-rect.der.rho(zlow[-1],zupp[-1],rhomat[j2,j3],corstr,j2,j3)
        prob2<-bvn(zlow[-1],zupp[-1],rhomat[j2,j3])
        prob3<-mrect.prob(zlow,zupp,rhomat[c(j1,j2,j3),c(j1,j2,j3)],bivpairs,corstr,mvncmp)
        s<-s+scj1*scj2j3/prob2*prob3
      }
    }}
  s
}

cov4<-function(j1,j2,j3,j4,mu,gam,rhomat,bivpairs,corstr,ub,link,mvncmp)
{ s<-0
  if(j1==j3 & j2==j4)
  { nvect=ub^2
    for(ii in 1:nvect)
    { y=d2v(2,ub,ii)+1
      vlow<-pmargmodel.ord(y-1,mu[c(j1,j2)],gam,link)
      tem<-dmargmodel.ord(y,mu[c(j1,j2)],gam,link)
      vupp<-vlow+tem
      zlow=qnorm(vlow)
      zupp=qnorm(vupp)
      zlow[zlow < -10]<--10
      zupp[zupp > 10]<-10
      scj1j2<-rect.der.rho(zlow,zupp,rhomat[j1,j2],corstr,j1,j2)
      prob<-bvn(zlow,zupp,rhomat[j1,j2])
      s<-s+scj1j2*scj1j2/prob
      #print(c(scj1j2,prob))
    }
  } else {
    if(j1==j3 & j2!=j4)
    { nvect=ub^3
      for(ii in 1:nvect)
      { y=d2v(3,ub,ii)+1
        vlow<-pmargmodel.ord(y-1,mu[c(j1,j2,j4)],gam,link)
        tem<-dmargmodel.ord(y,mu[c(j1,j2,j4)],gam,link)
        vupp<-vlow+tem
        zlow=qnorm(vlow)
        zupp=qnorm(vupp)
        zlow[zlow < -10]<--10
        zupp[zupp > 10]<-10
        scj1j2<-rect.der.rho(zlow[-3],zupp[-3],rhomat[j1,j2],corstr,j1,j2)
        prob21<-bvn(zlow[-3],zupp[-3],rhomat[j1,j2])
        scj3j4<-rect.der.rho(zlow[-2],zupp[-2],rhomat[j3,j4],corstr,j3,j4)
        prob22<-bvn(zlow[-2],zupp[-2],rhomat[j3,j4])
        prob3<-mrect.prob(zlow,zupp,rhomat[c(j1,j2,j4),c(j1,j2,j4)],
                          bivpairs,corstr,mvncmp)
        s<-s+scj1j2*scj3j4/prob21/prob22*prob3
      }
    } else {
      if(j1==j4 & j2!=j3)
      { nvect=ub^3
        for(ii in 1:nvect)
        { y=d2v(3,ub,ii)+1
          vlow<-pmargmodel.ord(y-1,mu[c(j1,j2,j3)],gam,link)
          tem<-dmargmodel.ord(y,mu[c(j1,j2,j3)],gam,link)
          vupp<-vlow+tem
          zlow=qnorm(vlow)
          zupp=qnorm(vupp)
          zlow[zlow < -10]<--10
          zupp[zupp > 10]<-10
          scj1j2<-rect.der.rho(zlow[-3],zupp[-3],rhomat[j1,j2],corstr,j1,j2)
          prob21<-bvn(zlow[-3],zupp[-3],rhomat[j1,j2])
          scj3j4<-rect.der.rho(zlow[c(3,1)],zupp[c(3,1)],rhomat[j3,j4],
                               corstr,j3,j4)
          prob22<-bvn(zlow[c(3,1)],zupp[c(3,1)],rhomat[j3,j4])
          prob3<-mrect.prob(zlow,zupp,rhomat[c(j1,j2,j3),c(j1,j2,j3)],
                            bivpairs,corstr,mvncmp)
          s<-s+scj1j2*scj3j4/prob21/prob22*prob3
        }
      } else {
        if(j2==j3 & j1!=j4)
        { nvect=ub^3
          for(ii in 1:nvect)
          { y=d2v(3,ub,ii)+1
            vlow<-pmargmodel.ord(y-1,mu[c(j1,j2,j4)],gam,link)
            tem<-dmargmodel.ord(y,mu[c(j1,j2,j4)],gam,link)
            vupp<-vlow+tem
            zlow=qnorm(vlow)
            zupp=qnorm(vupp)
            zlow[zlow < -10]<--10
            zupp[zupp > 10]<-10
            scj1j2<-rect.der.rho(zlow[-3],zupp[-3],rhomat[j1,j2],corstr,j1,j2)
            prob21<-bvn(zlow[-3],zupp[-3],rhomat[j1,j2])
            scj3j4<-rect.der.rho(zlow[-1],zupp[-1],rhomat[j3,j4],corstr,j3,j4)
            prob22<-bvn(zlow[-1],zupp[-1],rhomat[j3,j4])
            prob3<-mrect.prob(zlow,zupp,rhomat[c(j1,j2,j4),c(j1,j2,j4)],
                              bivpairs,corstr,mvncmp)
            s<-s+scj1j2*scj3j4/prob21/prob22*prob3
          }
        } else {
          if(j2==j4 & j1!=j3)
          { nvect=ub^3
            for(ii in 1:nvect)
            { y=d2v(3,ub,ii)+1
              vlow<-pmargmodel.ord(y-1,mu[c(j1,j2,j3)],gam,link)
              tem<-dmargmodel.ord(y,mu[c(j1,j2,j3)],gam,link)
              vupp<-vlow+tem
              zlow=qnorm(vlow)
              zupp=qnorm(vupp)
              zlow[zlow < -10]<--10
              zupp[zupp > 10]<-10
              scj1j2<-rect.der.rho(zlow[-3],zupp[-3],rhomat[j1,j2],corstr,j1,j2)
              prob21<-bvn(zlow[-3],zupp[-3],rhomat[j1,j2])
              scj3j4<-rect.der.rho(zlow[c(3,2)],zupp[c(3,2)],rhomat[j3,j4]
                                   ,corstr,j3,j4)
              prob22<-bvn(zlow[c(3,2)],zupp[c(3,2)],rhomat[j3,j4])
              prob3<-mrect.prob(zlow,zupp,rhomat[c(j1,j2,j3),c(j1,j2,j3)],
                                bivpairs,corstr,mvncmp)
              s<-s+scj1j2*scj3j4/prob21/prob22*prob3
            }
          } else
          { nvect=ub^4
            for(ii in 1:nvect)
            { y=d2v(4,ub,ii)+1
              vlow<-pmargmodel.ord(y-1,mu[c(j1,j2,j3,j4)],gam,link)
              tem<-dmargmodel.ord(y,mu[c(j1,j2,j3,j4)],gam,link)
              vupp<-vlow+tem
              zlow=qnorm(vlow)
              zupp=qnorm(vupp)
              zlow[zlow < -10]<--10
              zupp[zupp > 10]<-10
              scj1j2<-rect.der.rho(zlow[1:2],zupp[1:2],rhomat[j1,j2],corstr,j1,j2)
              prob21<-bvn(zlow[1:2],zupp[1:2],rhomat[j1,j2])
              scj3j4<-rect.der.rho(zlow[3:4],zupp[3:4],rhomat[j3,j4],corstr,j3,j4)
              prob22<-bvn(zlow[3:4],zupp[3:4],rhomat[j3,j4])
              prob4<-mrect.prob(zlow,zupp,rhomat[c(j1,j2,j3,j4),c(j1,j2,j3,j4)],
                                bivpairs,corstr,mvncmp)
              s<-s+scj1j2*scj3j4/prob21/prob22*prob4
            }
          }}}}}
  s
}

subselect2<-function(pairs,t)
{ 
  lr<-nrow(pairs)
  res=NULL
  for(j in 1:lr)
  { if(sum(pairs[j,1]==t)==1 & sum(pairs[j,2]==t)==1) res=c(res,j)
  }
  res
}




clic<-function(nbcl,r,b,gam,xdat,id,tvec,corstr,WtScMat,link,mvncmp)
{ if(is.matrix(xdat))
{ dim<-dim(xdat)
  n<-1:dim[1]
  p<-dim[2]
} else {n<-1:length(xdat); p<-1}

uid<-unique(id)
d<-id.size(id)
maxd<-max(d)
q<-length(gam)

tvec<-id.time(tvec,d)
bivpairs<-maxpairs(d)
d2<-nrow(bivpairs)
rmat<-cormat(maxd,r,bivpairs,corstr)
t=p+q+d2
Dmat<-M<-matrix(0,t,t)
m<-0
for(i in uid)
{ m<-m+1
  cases<-id==i
  irow=n[cases]
  ti<-tvec[irow]
  if(length(ti)>1)
  {
  newmui<-rep(NA,maxd)
  if(is.matrix(xdat)) 
  { x<-xdat[cases,]
    newx<-matrix(NA,maxd,p) 
    newx[ti,]<-x}  else {
      x<-xdat[cases]
      newx<-rep(NA,maxd) 
      newx[ti]<-x 
    }
  
  
  mui<-ordreg.mu(x,b)
  newmui[ti]<-mui
  ub<-q+1
  nvect=ub^maxd
  nvect2=ub^2
  
  y=vlow=tem=matrix(NA,nvect,maxd)
  for(ii in 1:nvect)
  { y[ii,]=d2v(maxd,ub,ii)+1
    vlow[ii,]<-pmargmodel.ord(y[ii,]-1,newmui,gam,link)
    tem[ii,]<-dmargmodel.ord(y[ii,],newmui,gam,link)
  }
  vupp<-vlow+tem
  zlow=qnorm(vlow)
  zupp=qnorm(vupp)
  zlow[zlow < -10]<--10
  zupp[zupp > 10]<-10
  dzlow=dnorm(zlow)
  dzupp=dnorm(zupp) 
  
  if(d[m]<maxd)
  { # new code
    und<-d[m]
    und2<-choose(und,2)
    unpairs<-maxpairs(und)
    unrmat<-rmat[ti,ti]
  
    
    D22i<-rep(0,und2)
    Delta21i<-matrix(NA,und2,und*q)
    for(i1 in 1:und2)
    { k1<-unpairs[i1,1]
      k2<-unpairs[i1,2]
      s<-rep(0,q)
      s2=0
        for(ii in 1:nvect2)
        { y2=d2v(2,ub,ii)+1
          vlow2<-pmargmodel.ord(y2-1,mui[c(k1,k2)],gam,link)
          tem2<-dmargmodel.ord(y2,mui[c(k1,k2)],gam,link)
          vupp2<-vlow2+tem2
          zlow2=qnorm(vlow2)
          zupp2=qnorm(vupp2)
          zlow2[zlow2 < -10]<--10
          zupp2[zupp2 > 10]<-10
          dzlow2=dnorm(zlow2)
          dzupp2=dnorm(zupp2) 
          der2<-rect.der.rho(zlow2,zupp2,unrmat[k1,k2],corstr,k1,k2)
          prob<-bvn(zlow2,zupp2,unrmat[k1,k2])
          s2<-s2+der2*der2/prob
          der1<-NULL
          for(i2 in 1:und)
          { if(i2==k1)
          { der1<-c(der1,rect.der.gam(zlow2,
                                      zupp2,dzlow2,dzupp2,
                                      y2,mui[c(k1,k2)],gam,unrmat[k1,k2],1,link))
          } else {
            if(i2==k2)
            { der1<-c(der1,rect.der.gam(zlow2,
                                        zupp2,dzlow2,dzupp2,
                                        y2,mui[c(k1,k2)],gam,unrmat[k1,k2],2,link)) 
            } else {
              der1<-c(der1,rep(0,q))}} # new change instead of c(der1,0)
          }
          s<-s+der1*der2/prob
        }
        D22i[i1]=s2
        Delta21i[i1,]<-s
      }
    
    
    
    Omega12i<-NULL #Omega12i<-matrix(NA,q*maxd,d2) 
    for(j1 in 1:und)
    { temp=NULL
      for(j2 in 1:und2)
      { temp=c(temp,cov3(j1,unpairs[j2,1],unpairs[j2,2],
                         mui,gam,unrmat,unpairs,corstr,ub,link,mvncmp))
      }
      subOmega12i=matrix(temp,ncol=und2)
      Omega12i<-rbind(Omega12i,subOmega12i) 
    } 
    
    Omega22i<-matrix(NA,und2,und2)
    if(und2==1){ Omega22i=D22i } else {
    for(j1 in 1:(und2-1))
    { for(j2 in (j1+1):und2)
    { temp<-cov4(unpairs[j1,1],unpairs[j1,2],unpairs[j2,1],unpairs[j2,2],
                 mui,gam,unrmat,unpairs,corstr,ub,link,mvncmp)
      Omega22i[j1,j2]<-temp
      Omega22i[j2,j1]<-temp
    }
    }}
    } else {
    
    D22i<-rep(0,d2)
    Delta21i<-matrix(NA,d2,maxd*q)
    for(i1 in 1:d2)
    { k1<-bivpairs[i1,1]
      k2<-bivpairs[i1,2]
      if((sum(k1==ti)==1) & (sum(k2==ti)==1))
      {
        s<-rep(0,q)
        s2=0
        for(ii in 1:nvect2)
        { y2=d2v(2,ub,ii)+1
          vlow2<-pmargmodel.ord(y2-1,newmui[c(k1,k2)],gam,link)
          tem2<-dmargmodel.ord(y2,newmui[c(k1,k2)],gam,link)
          vupp2<-vlow2+tem2
          zlow2=qnorm(vlow2)
          zupp2=qnorm(vupp2)
          zlow2[zlow2 < -10]<--10
          zupp2[zupp2 > 10]<-10
          dzlow2=dnorm(zlow2)
          dzupp2=dnorm(zupp2) 
          der2<-rect.der.rho(zlow2,zupp2,rmat[k1,k2],corstr,k1,k2)
          prob<-bvn(zlow2,zupp2,rmat[k1,k2])
          s2<-s2+der2*der2/prob
          der1<-NULL
          for(i2 in 1:maxd)
          { if(i2==k1)
          { der1<-c(der1,rect.der.gam(zlow2,
                                      zupp2,dzlow2,dzupp2,
                                      y2,newmui[c(k1,k2)],gam,rmat[k1,k2],1,link))
          } else {
            if(i2==k2)
            { der1<-c(der1,rect.der.gam(zlow2,
                                        zupp2,dzlow2,dzupp2,
                                        y2,newmui[c(k1,k2)],gam,rmat[k1,k2],2,link)) 
            } else {
              der1<-c(der1,rep(0,q))}} # new change instead of c(der1,0)
          }
          s<-s+der1*der2/prob
        }
        D22i[i1]=s2
        Delta21i[i1,]<-s
      }
    }
    
    
  Omega12i<-NULL #Omega12i<-matrix(NA,q*maxd,d2) 
  for(j1 in 1:maxd)
  { temp=NULL
    for(j2 in 1:d2)
    { temp=c(temp,cov3(j1,bivpairs[j2,1],bivpairs[j2,2],
                       newmui,gam,rmat,bivpairs,corstr,ub,link,mvncmp))
    }
    subOmega12i=matrix(temp,ncol=d2)
    Omega12i<-rbind(Omega12i,subOmega12i) 
  } 
  
  Omega22i<-matrix(NA,d2,d2)
  for(j1 in 1:(d2-1))
  { for(j2 in (j1+1):d2)
  { temp<-cov4(bivpairs[j1,1],bivpairs[j1,2],bivpairs[j2,1],bivpairs[j2,2],
               newmui,gam,rmat,bivpairs,corstr,ub,link,mvncmp)
    Omega22i[j1,j2]<-temp
    Omega22i[j2,j1]<-temp
  }
  }
  }
  
  seli<-subselect.ord(ti,q)
  if(d[m]>2) diag(Omega22i)=D22i
  #D22i<-diag(diag(Omega22i)) for j1=j2 
  #D22i=diag(D22i)
  Omega11i<-WtScMat$omega[,,m]
  Omega11i<-Omega11i[seli,seli]
  
  sel2i=subselect2(bivpairs,ti)
  
  Xi<-WtScMat$X[,,m]
  Xi<-Xi[,seli]
  tXi<-t(Xi)
  
  D21i<-Delta21i%*%tXi
  
  temp1=matrix(0,d2,p+q)
  temp1[sel2i,]=D21i
  D21i=temp1
  temp2=rep(0,d2)
  temp2[sel2i]=D22i
  D22i=diag(temp2)
  
  Delta11<-WtScMat$delta[,,m] # from the wtsc code: weightMat()
  Delta11<-Delta11[seli,seli]
  D11i<-Xi%*%Delta11%*%tXi
  
  d0<-dim(t(D21i))
  D12i<-matrix(0,d0[1],d0[2])
  Di<-rbind(cbind(D11i,D12i),cbind(D21i,D22i))
  
  
  M11i<-Xi%*%Omega11i%*%tXi
  M12i<-Xi%*%Omega12i
  temp1=matrix(0,p+q,d2)
  temp1[,sel2i]=M12i
  M12i=temp1
  M21i<-t(M12i)
  M22i<-Omega22i
  temp2=matrix(0,d2,d2)
  temp2[sel2i,sel2i]=M22i
  M22i=temp2
  Mi<-rbind(cbind(M11i,M12i),cbind(M21i,M22i))
  
  ##########################
  Dmat<-Dmat+Di
  M<-M+Mi
}}
#inDmat<-solve(Dmat)
#tr<-sum(diag(M%*%inDmat))
tr<-sum(diag(solve(Dmat,M)))
list(AIC=2*(nbcl+tr),BIC=2*nbcl+log(m)*tr)
}


