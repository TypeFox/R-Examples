

# Density  of the univariate marginal distribution
# input:
# y the vector of (non-negative integer) quantiles.
# mu the mean parameter of the univariate distribution.
# gam the parameter gamma of  the negative binomial distribution.
# output:
# the density of the univariate marginal  distribution
dmargmodel.ord<-function(y,mu,gam,link)
{ cuts<-c(-10,gam,10)
lb=cuts[y]+mu
ub=cuts[y+1]+mu
if(link=="probit") res<-pnorm(ub)-pnorm(lb) else res<-plogis(ub)-plogis(lb)
res[y<1]<-0
res
}

# CDF  of the univariate marginal distribution
# input:
# y the vector of (non-negative integer) quantiles.
# mu the mean parameter of the univariate distribution.
# gam the parameter gamma of  the negative binomial distribution.
# output:
# the cdf of the univariate marginal  distribution
pmargmodel.ord<-function(y,mu,gam,link)
{ cuts<-c(-10,gam,10)
  ub=cuts[y+1]+mu # for mprobit
  #ub=cuts[y+1]-mu # for polr
  if(link=="probit") res<-pnorm(ub) else res<-plogis(ub)
  res[y<1]<-0
  res
}



# Bivariate composite likelihood for multivariate normal copula with ordinal regression.
# input:
# r the vector of normal copula parameters
# b the regression coefficients
# gam the gamma parameter
# xdat the matrix of covariates (use the constant 1 for the first covariate)
# ydat the vector with the response
# id the the vector with the id
# tvec the time related vector
# link is the link function. Choices are ?log? for the log link function, 
# ?logit? for the logit link function, and ?probit? for the probit link function.
# corstr indicates the latent correlation structure of normal copula. 
# Choices are exch, ar, and unstr for exchangeable, ar(1) and 
# unstrucutred correlation structure, respectively. 
# output:
# negative bivariate composite likelihood for multivariate normal copula 
# with ordinal regression.
bcl.ord<-function(r,b,gam,xdat,ydat,id,tvec,corstr,link)
{ s<-0
  uid<-unique(id)
  d<-id.size(id)
  maxd<-max(d)
  tvec<-id.time(tvec,d)
  n<-1:length(ydat)    #### change here
  pairs<-maxpairs(d)
  rmat<-cormat(maxd,r,pairs,corstr)            
  c1<-(r[1] > 1 | r[1] < (-1/(maxd-1))) 
  c2<-(r[1] > 1 | r[1] < -1)            
  c3<-(det(rmat)<0 | min(rmat)< -1 | max(rmat)>1)
  if((corstr=="exch" & c1) | (corstr=="ar" & c2 ) | (corstr=="unstr" & c3)) {s<-1e10}
  else {  
  for(i in uid)
  { cases<-id==i
    irow=n[cases]
    yi<-ydat[irow]
    ti<-tvec[irow]
    newyi<-rep(NA,maxd)
    newmui<-rep(NA,maxd)
    newyi[ti]<-yi
    if(is.vector(xdat)) x<-xdat[cases]  else  x<-xdat[cases,]    ###change here
    mui<-ordreg.mu(x,b)
    newmui[ti]<-mui
    vlow<-pmargmodel.ord(newyi-1,newmui,gam,link)
    tem<-dmargmodel.ord(newyi,newmui,gam,link)
    vupp<-vlow+tem
    zlow=qnorm(vlow)
    zupp=qnorm(vupp)
    for(j in 1:dim(pairs)[1])
    { k1<-pairs[j,][1]
      k2<-pairs[j,][2]
      if((sum(k1==ti)==1) & (sum(k2==ti)==1))
      { if(corstr=="exch")
      { s<-s +bivlik(zlow[c(k1,k2)],zupp[c(k1,k2)],r) }
      else
      { if(corstr=="ar")
      { s<-s +bivlik(zlow[c(k1,k2)],zupp[c(k1,k2)],r^(k2-k1)) }
      else
      { # corstr=="unstr"
        s<-s +bivlik(zlow[c(k1,k2)],zupp[c(k1,k2)],r[j])
      }}}
  else { s<-s }
  }
  }
  -s
}
}

# optimization routine for composite likelihood for MVN copula
# b the regression coefficients
# gam the gamma parameter
# xdat the matrix of covariates (use the constant 1 for the first covariate)
# ydat the vector with the response
# id the the vector with the id
# tvec the time related vector
# link is the link function. Choices are  
# ?logit? for the logit link function, and ?probit? for the probit link function.
# corstr indicates the latent correlation structure of normal copula. 
# Choices are ?exch?, ?ar?, and ?unstr? for exchangeable, ar(1) and 
# unstrucutred correlation structure, respectively. 
# output: A list containing the following components:
# minimum the value of the estimated minimum of CL1 for MVN copula
# estimate the CL1 estimates
# gradient the gradient at the estimated minimum of CL1
# code an integer indicating why the optimization process terminated, see nlm.
cl1.ord<-function(b,gam,xdat,ydat,id,tvec,corstr,link)
{ d<-id.size(id)
  pairs<-maxpairs(d)
  if(corstr=="unstr")
  { nom<-dim(pairs)[1]
    nlm(bcl.ord,rep(0.1,nom),b,gam,xdat,ydat,id,tvec,corstr,link)
  }
  else
  { nlm(bcl.ord,0.1,b,gam,xdat,ydat,id,tvec,corstr,link) }
}


# derivative of the ordinal loglikelihood with respect to gamma
# input:
# gam the gamma parameter
# mu the mean parameter
# output:
# the vector with the derivatives of the NB loglikelihood with respect to gamma
derlik.gam.ord<-function(mu,gam,u,link)
{ K<-length(gam)+1
  k<-1:K
  cuts<-c(-10,gam,10)
  lb=cuts[k]+mu
  ub=cuts[k+1]+mu
  if(link=="probit") { dlatent=dnorm; platent=pnorm } else { dlatent=dlogis; platent=plogis }
  dlatentub<-dlatent(ub)
  dlatentlb<-dlatent(lb)
  den<-platent(ub)-platent(lb)
  res<-rep(NA,K)
  for(i in 1:K)
  { if(u==i)
  { res[i]=dlatentub[i]/den[i] }
  else
  { if(u==i-1)
  { res[i]=-dlatentlb[i]/den[i] }
    else {res[i]=0}}
}
res
}

# derivative of the NB loglikelihood with respect to gamma
# input:
# gam the gamma parameter
# mu the mean parameter
# y  the value of a non-negative integer quantile
# output:
# the derivative of the NB loglikelihood with respect to gamma
iderlik.gam.ord<-function(mu,gam,y,u,link)
{ cuts<-c(-10,gam,10)
  lb=cuts[y]+mu
  ub=cuts[y+1]+mu
  if(link=="probit") { dlatent=dnorm; platent=pnorm } else { dlatent=dlogis; platent=plogis }
  den<-platent(ub)-platent(lb)
  if(u==y) dlatent(ub)/den
  else if(u==y-1) -dlatent(lb)/den  else 0
}


der.dnorm<-function(x)
{ -x*dnorm(x) }

der.dlogis<-function(x)
{ expx=exp(x)
  expx*(1-expx)/(1+expx)^3 
}

# minus expectation of the second derivative of the marginal ordinal loglikelihood
# with resect to gamma
# input:
# mu the mean parameter
# gam the gamma parameter
# u the univariate cdfs
# output:
# the vector with the minus expectations of the margmodel loglikelihood
# with respect to gamma
fisher.gam.ord<-function(mu,gam,u,v,link)
{ cuts<-c(-10,gam,10)
  K<-length(gam)+1
  k<-1:K
  lb=cuts[k]+mu
  ub=cuts[k+1]+mu
  if(link=="probit") { dlatent=dnorm; platent=pnorm; der.dlatent=der.dnorm } else { 
    dlatent=dlogis; platent=plogis; der.dlatent=der.dlogis }
  den<-platent(ub)-platent(lb)
  dlatentub<-dlatent(ub)
  dlatentlb<-dlatent(lb)
  der.dlatentub<-der.dlatent(ub)
  der.dlatentlb<-der.dlatent(lb)
  h<-rep(NA,K)
  for(k in 1:K)
  { if(u==k & v==k)
  { num1<-der.dlatentub[k]
    num2<-dlatentub[k]
    tem<-num2/den[k]
    h[k]<--num1/den[k]+tem*tem  
  }
  else
  { if((u==k & v==k-1) | (u==k-1 & v==k))
  { h[k]<--dlatentub[k]*dlatentlb[k]/den[k]/den[k]  
  }
    else
    { if(u==k-1 & v==k-1)
    { num1<-der.dlatentlb[k]
    num2<-dlatentlb[k]
    tem<-num2/den[k]
    h[k]<-num1/den[k]+tem*tem  
    } else h[k]<-0}}
 }
 sum(h*den) 
}

# the mean values of the univariate marginal distribution 
# corresonding to the used link function 
# input:
# x the matix of the covariates 
# b the vector with the regression coefficients
# output:
# the mean values of the univariate marginal distribution 
ordreg.mu<-function(x,b)
{ if(length(b)!=1) mu<-x %*% b else mu<-x*b }


# select the present column and lines in Omega, Delta and X matrices
# for unbalanced data
# input:
# tvec a vector of the time for an individual
subselect.ord<-function(tvec,q)
{ sel<-NULL
for(i in 1:length(tvec))
{ tm<-tvec[i]
k<-((tm-1)*q+1):((tm-1)*q+q)
sel<-c(sel,k)
}
sel
}

# Calculating the number of categories
# input:
# gam the cutpoints
# output:
# the number of categories
noCategories<-function(gam)
{ length(gam)+1 }

# covariance matrix of the scores Omega_i
# input:
# scgam the array of the score functions with respect to gam
# index the bivariate pair
# pmf the matrix of rectangle probabilities
scoreCov.ord<-function(scgam,pmf,index)
{ j1<-index[1]
  j2<-index[2]
  q<-dim(scgam)[2]
  cov22<-matrix(NA,q,q)
  for(i1 in 1:q)
  { for(i2 in 1:q)
  { cov22[i1,i2]<-t(scgam[,i1,j1])%*%pmf%*%scgam[,i2,j2] }
  }
  cov22
}

# weight matrix fixed at values from the CL1 estimator
# input:
# b the regression coefficients
# gam the gamma parameter
# xdat the matrix of covariates (use the constant 1 for the first covariate)
# ydat the vector with the response
# id the the vector with the id
# tvec the time related vector
# rh the vector with CL1 estimates
# link is the link function. Choices are 
# ?logit? for the logit link function, and ?probit? for the probit link function.
# corstr indicates the latent correlation structure of normal copula. 
# Choices are ?exch?, ?ar?, and ?unstr? for exchangeable, ar(1) and 
# unstrucutred correlation structure, respectively. 
# output: A list containing the following components:
# omega the array with the Omega matrices
# delta the array with the Delta matrices
# X the array with the X matrices
weightMat.ord<-function(b,gam,rh,xdat,ydat,id,tvec,corstr,link)
{ uid<-unique(id)
  d<-id.size(id)
  maxd<-max(d)
  q<-length(gam)
  lid<-length(uid)
  if(is.matrix(xdat))
  { dim<-dim(xdat)
    n<-1:dim[1]
    p<-dim[2]
   } else {n<-1:length(xdat); p<-1}
   pairs<-maxpairs(d)
   tvec<-id.time(tvec,d)
   omega<-array(NA,c(maxd*q,maxd*q,lid))
   X<-array(NA,c(p+q,maxd*q,lid))
   delta<-array(NA,c(maxd*q,maxd*q,lid))
   dom<-maxd*q
   if(q>1) pos<-seq(1,dom-1,by=q) else pos<-seq(1,dom)    #for binary
   m<-0
   for(i in uid)
   { m<-m+1
     cases<-id==i
     irow=n[cases]
     newyi<-rep(NA,maxd)
     newmui<-rep(NA,maxd)
     ti<-tvec[irow]
     if(is.matrix(xdat)) 
     { x<-xdat[cases,]
       newx<-matrix(NA,maxd,p) 
       newx[ti,]<-x}  else {
       x<-xdat[cases]
       newx<-rep(NA,maxd) 
       newx[ti]<-x }
       mui<-ordreg.mu(x,b)
       newmui[ti]<-mui
       ub<-noCategories(gam)
       du<-matrix(NA,ub,maxd)
       scgam<-array(NA,c(ub,ub-1,maxd))
       for(j in ti)
       { du[,j]<-dmargmodel.ord(1:ub,newmui[j],gam,link)
         for(k in 1:(ub-1))
         { scgam[,k,j]<-derlik.gam.ord(newmui[j],gam,k,link) } 
        }
      u<-apply(du,2,cumsum)
      z<-qnorm(u)
      z[is.nan(z)]<-7
      z[z>10]<-7
      z<-rbind(-7,z)
      pz<-pnorm(z)
      dz<-dnorm(z)
      zs<-z*z
      zc<-z*zs
      zf<-zs*zs
      xi<-NULL
      diagonali<-array(NA,c(q,q,maxd))   
      for(j in 1:maxd)
      { if(is.matrix(newx)) 
      { temp1<-matrix(rep(newx[j,],each=q),q) }  else { 
        temp1<-matrix(rep(newx[j],each=q),q) }
      xi<-cbind(xi,t(cbind(temp1,diag(q))))
      fisher<-matrix(NA,ub-1,ub-1)
      for(k1 in 1:(ub-1))
      { for(k2 in 1:(ub-1))
      { fisher[k1,k2]<-fisher.gam.ord(newmui[j],gam,k1,k2,link)
      }
      }
      diagonali[,,j]<-fisher
      }
     deltai<-matrix(0,dom,dom)
     minus<-0
     for(j in pos)
     { deltai[j:(j+q-1),j:(j+q-1)]<-diagonali[,,j-minus]
       if(q>1) minus<-minus+q-1 else minus<-0
     }
     offi<-array(NA,c(q,q,dim(pairs)[1]))
     for(k in 1:dim(pairs)[1])
     { k1<-pairs[k,][1]
       k2<-pairs[k,][2]
       if((sum(k1==ti)==1) & (sum(k2==ti)==1))
       { x1=z[,k1]; x2=z[,k2]
         x1s=zs[,k1]; x2s=zs[,k2]
         x1c=zc[,k1]; x1f=zf[,k1]
         x2c=zc[,k2]; x2f=zf[,k2]
         t1=outer(pz[,k1],pz[,k2])
         t2=outer(dz[,k1],dz[,k2])
         if(corstr=="exch")
         { cdf<-approxbvncdf(rh,x1,x2,x1s,x2s,x1c,x2c,x1f,x2f,t1,t2) }   else { 
         if(corstr=="ar")
         { cdf<-approxbvncdf(rh^(k2-k1),x1,x2,x1s,x2s,x1c,x2c,x1f,x2f,t1,t2) } else { # corstr=="unstr"
           cdf<-approxbvncdf(rh[k],x1,x2,x1s,x2s,x1c,x2c,x1f,x2f,t1,t2)
         }
      }
      cdf1=apply(cdf,2,diff)
      pmf=apply(t(cdf1),2,diff)
      pmf=t(pmf)
      offi[,,k]<-scoreCov.ord(scgam,pmf,pairs[k,])
     }
     }
     omegai<-deltai
     ch1<-0
     ch2<-0
     ch3<-0
     if(d[m]>1)
     {
    for(j in 1:(d[m]-1))
    { for(r in pos[-(1:j)])
    { #print(c(j,r))
      omegai[(1+(j-1)*q):(j*q),r:(r+q-1)]<-offi[,,(j+ch2-ch1)]
      omegai[r:(r+q-1),(1+(j-1)*q):(j*q)]<-t(offi[,,(j+ch2-ch1)])
      ch2<-ch2+1
    }
      ch1<-ch1+1
    }
  }
  X[,,m]<-xi
  delta[,,m]<-deltai
  omega[,,m]<-omegai
  }
  list(omega=omega,X=X,delta=delta)
}


# the weigted scores equations
# input:
# param the vector of regression and not regression parameters
# xdat the matrix of covariates 
# ydat the vector with the response
# id the the vector with the id
# tvec the time related vector
# link is the link function. Choices are  
# ?logit? for the logit link function, and ?probit? for the probit link function.
# corstr indicates the latent correlation structure of normal copula. 
# Choices are ?exch?, ?ar?, and ?unstr? for exchangeable, ar(1) and 
# unstrucutred correlation structure, respectively. 
# WtScMat is a list containing the following components:
# omega the array with the Omega matrices
# delta the array with the Delta matrices
# X the array with the X matrices
# output
# the weigted scores equations
wtsc.ord<-function(param,WtScMat,xdat,ydat,id,tvec,link)
{ uid<-unique(id)
  d<-id.size(id)
  maxd<-max(d)
  if(is.matrix(xdat))
  { dim<-dim(xdat)
    n<-1:dim[1]
    p<-dim[2]
  } else {n<-1:length(xdat); p<-1}
  tvec<-id.time(tvec,d)
  b<-param[1:p]
  q<-length(unique(ydat))-1
  gam<-param[(p+1):(p+q)]
  g<-0
  m<-0
  for(i in uid)
  { #print(i)
    cases<-id==i
    irow=n[cases]
    m<-m+1
    if(is.matrix(xdat)) x<-xdat[cases,] else x<-xdat[cases]
    mu<-ordreg.mu(x,b)
    ub<-noCategories(gam)
    scgam<-array(NA,c(ub,ub-1,d[m]))
    for(j in 1:d[m])
    { for(k in 1:(ub-1))
    { scgam[,k,j]<-derlik.gam.ord(mu[j],gam,k,link) }
    }        
    y<-ydat[irow]
    sci<-NULL
    for(j in 1:d[m])
    { scgami<-NULL
      for(k in 1:(ub-1))
      { scgami<-c(scgami,scgam[y[j],k,j])}
      sci<-c(sci,scgami)
    }
    ti<-tvec[irow]
    seli<-subselect.ord(ti,q)
    Xi<-WtScMat$X[,,m]
    Xi<-Xi[,seli]
    deltai<-WtScMat$delta[,,m]
    deltai<-deltai[seli,seli]
    omegai<-WtScMat$omega[,,m]
    omegai<-omegai[seli,seli]
    gi<-Xi%*%t(deltai)%*%solve(omegai,sci)
    g<-g+gi
    }
  g
}


# solving the weigted scores equations
# input:
# start the starting values (IEE estimates) for the vector of
# regression and not regression parameters
# xdat the matrix of covariates (use the constant 1 for the first covariate)
# ydat the vector with the response
# id the the vector with the id
# tvec the time related vector
# link is the link function. Choices are 
# ?logit? for the logit link function, and ?probit? for the probit link function.
# corstr indicates the latent correlation structure of normal copula. 
# Choices are ?exch?, ?ar?, and ?unstr? for exchangeable, ar(1) and 
# unstrucutred correlation structure, respectively. 
# WtScMat is a list containing the following components:
# omega the array with the Omega matrices
# delta the array with the Delta matrices
# X the array with the X matrices
# output:
# the weighted scores estimates
solvewtsc.ord<-function(start,WtScMat,xdat,ydat,id,tvec,link)
{ multiroot(f=wtsc.ord,start,atol=1e-4,rtol=1e-4,ctol=1e-4,
            WtScMat=WtScMat,xdat=xdat,ydat=ydat,id=id,tvec=tvec,link=link) 
}


# inverse Godambe matrix with Delta and Omega evaluated at IEE estimator
# input:
# b the regression coefficients
# gam the gamma parameter
# xdat the matrix of covariates (use the constant 1 for the first covariate)
# ydat the vector with the response
# id the the vector with the id
# tvec the time related vector
# rh the vector with CL1 estimates
# WtScMat a list containing the following components:
# omega the array with the Omega matrices
# delta the array with the Delta matrices
# X the array with the X matrices
# link is the link function. Choices are
# ?logit? for the logit link function, and ?probit? for the probit link function.
# output:
# the inverse Godambe matrix
godambe.ord<-function(param,WtScMat,xdat,ydat,id,tvec,link)
{ uid<-unique(id)
  d<-id.size(id)
  maxd<-max(d)
  if(is.matrix(xdat))
  { dim<-dim(xdat)
    n<-1:dim[1]
    p<-dim[2]
  } else {n<-1:length(xdat); p<-1}
  tvec<-id.time(tvec,d)
  b<-param[1:p]
  q<-length(unique(ydat))-1
  gam<-param[(p+1):(p+q)]
  v<-v1<-v2<-v3<-0
  fv<-fv1<-fv2<-fv3<-0
  m<-0
  for(i in uid)
  { m<-m+1
    cases<-id==i
    irow=n[cases]
    ti<-tvec[irow]
    if(is.matrix(xdat))
    { x<-xdat[cases,]
      newx<-matrix(NA,maxd,p)
      newx[ti,]<-x}  else {
      x<-xdat[cases]
      newx<-rep(NA,maxd)
      newx[ti]<-x }
      newy<-rep(NA,maxd)
      newmui<-rep(NA,maxd)
      mui<-ordreg.mu(x,b)
      newmui[ti]<-mui
      ub<-noCategories(gam)
      scgam<-array(NA,c(ub,ub-1,maxd))
      for(j in ti)
      { for(k in 1:(ub-1))
      { scgam[,k,j]<-derlik.gam.ord(newmui[j],gam,k,link) }
      }
      y<-ydat[irow]
      newy[ti]<-y
      sci<-NULL
      for(j in 1:maxd) 
      { scgami<-NULL
        for(k in 1:(ub-1))
        { scgami<-c(scgami,scgam[newy[j],k,j])}
          sci<-c(sci,scgami)
        }
      sci<-sci[!is.na(sci)]
      seli<-subselect.ord(ti,q)
      Xi<-WtScMat$X[,,m]
      Xi<-Xi[,seli]
      fdeltai<-WtScMat$delta[,,m]
      fdeltai<-fdeltai[seli,seli]
      fomegai<-WtScMat$omega[,,m]
      fomegai<-fomegai[seli,seli]
      fci<-fdeltai%*%t(Xi)
      nrhs1=ncol(fdeltai)
      nrhs2=ncol(fci)
      tem.lin=solve(fomegai,cbind(fdeltai,fci))
      finvai=t(tem.lin[,1:nrhs1])
      xi.invai=Xi %*% finvai
      fai<-solve(finvai)
      tem=solve(t(fai),t(Xi))
      fv1i<-xi.invai %*% fci
      fv2i<-xi.invai %*% (sci%*% t(sci))  %*% tem
      fv3i<-t(fci) %*% tem
      fv1<-fv1+fv1i
      fv2<-fv2+fv2i
      fv3<-fv3+fv3i
  }
  solve(fv1,fv2)%*%solve(fv3)
}

# the weighted scores wrapper function: handles all the steps in the weighted scores
# xdat the matrix of covariates (use the constant 1 for the first covariate)
# ydat the vector with the response
# id the the vector with the id
# tvec the time related vector
# rh the vector with CL1 estimates
# WtScMat a list containing the following components:
# omega the array with the Omega matrices
# delta the array with the Delta matrices
# X the array with the X matrices
# corstr indicates the latent correlation structure of normal copula. 
# Choices are ?exch?, ?ar?, and ?unstr? for exchangeable, ar(1) and 
# unstrucutred correlation structure, respectively. 
# link is the link function. Choices are 
# ?logit? for the logit link function, and ?probit? for the probit link function.
# iprint indicates printing of some intermediate results, default FALSE
wtsc.ord.wrapper<-function(xdat,ydat,id,tvec,corstr,link,iprint=FALSE)
{ 
i.est<-iee.ord(xdat,ydat,link)
if(iprint)
{ cat("\niest: IEE estimates\n")
  print(c(i.est$reg,i.est$gam))
}
est.rho<-cl1.ord(b=i.est$reg,gam=i.est$gam,xdat,ydat,id,tvec,corstr,link)
if(iprint)
{ cat("\nest.rho: CL1 estimates\n")
  print(est.rho$e)
  cat("\nest.rho: CL1 likelihood\n")
  print(-est.rho$m)
}
WtScMat<-weightMat.ord(b=i.est$reg,gam=i.est$gam,rh=est.rho$e,
                   xdat,ydat,id,tvec,corstr,link)
ws<-solvewtsc.ord(start=c(i.est$reg,i.est$gam),WtScMat,xdat,ydat,id,
              tvec,link)
if(iprint)
{ cat("ws=parameter estimates\n")
  print(ws$r)
}
acov<-godambe.ord(ws$r,WtScMat,xdat,ydat,id,tvec,link)
se<-sqrt(diag(acov))
if(iprint)
{ cat("\nacov: inverse Godambe matrix with W based on first-stage wt matrices\n")
  print(acov)
  cat("\nse: robust standard errors\n")
  print(se)
  res<-round(cbind(ws$r,se),3)
  cat("\nres: Weighted scores estimates and standard errors\n")
  print(res)
}
list(IEEest=c(i.est$reg,i.est$gam),CL1est=est.rho$e,CL1lik=-est.rho$m,WSest=ws$r, asympcov=acov)
}


