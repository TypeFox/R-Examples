# ML for ordinal probit model, using modified Newton-Raphson



iee.ord<- function(x,y,link,iprint=0,maxiter=20,toler=1.e-6)
{
  if(!is.vector(x))
  { if(nrow(x)!=length(y)) stop("x, y not same length") }     
  else if(length(x)!=length(y)) { stop("x, y not same length") }

  if(is.vector(x)) x=as.matrix(x)
  n=length(y)
  # assume y in 1,...,norc
  norc=length(unique(y))
  npred=ncol(x)
  np=norc-1+npred
  # centering of x so that can use start of 0 for beta
  xmn=apply(x,2,"mean")
  xc=scale(x,center=xmn,scale=F)

  # starting point for NR
  cum=(1:(norc-1))
  cutp=rep(0,norc-1)
  for(k in cum)
  { pr=sum(y<=k)
    if (pr==0) { pr=1 }
    cutp[k]=qlogis(pr/n)
  }
  b=rep(0,npred)
  if(link=="probit") { dlatent=dnorm; platent=pnorm; der.dlatent=der.dnorm } else { 
    dlatent=dlogis; platent=plogis; der.dlatent=der.dlogis }
  # loop
  mxdif=1
  iter=0
  while(iter<maxiter & mxdif>toler)
  { tem=xc%*%b
    cutb=c(-10,cutp,10)
    ub=tem+cutb[y+1]
    lb=tem+cutb[y]
    ucdf=platent(ub)
    lcdf=platent(lb)
    updf=dlatent(ub)
    lpdf=dlatent(lb)
    # score vector
    dbeta=rep(0,npred)
    dcut=rep(0,norc+1)
    # Hessian matrix
    d2beta=matrix(0,npred,npred)
    d2bcut=matrix(0,npred,norc+1)
    d2cut=matrix(0,norc+1,norc+1)
    for(i in 1:n)
    { uderi=der.dlatent(ub[i])
      lderi=der.dlatent(lb[i])
      xx=xc[i,]
      pri=ucdf[i]-lcdf[i]
      prderi=updf[i]-lpdf[i]
      dbeta=dbeta+xx*prderi/pri 
      k=y[i]
      dcut[k+1]=dcut[k+1]+updf[i]/pri
      dcut[k]=dcut[k]-lpdf[i]/pri
      pr2=pri^2
      d2beta=d2beta+ outer(xx,xx)*((uderi-lderi)*pri-prderi^2)/pr2 
      d2bcut[,k+1]=d2bcut[,k+1]+xx*(uderi*pri-updf[i]*prderi)/pr2
      d2bcut[,k]=d2bcut[,k]+xx*(-lderi*pri+lpdf[i]*prderi)/pr2
      d2cut[k+1,k+1]=d2cut[k+1,k+1]+ (uderi*pri-updf[i]^2)/pr2
      d2cut[k,k]=d2cut[k,k]+ (-lderi*pri-lpdf[i]^2)/pr2
      tem2=updf[i]*lpdf[i]/pr2
      d2cut[k,k+1]=d2cut[k,k+1]+tem2
      d2cut[k+1,k]=d2cut[k+1,k]+tem2
    }  

    sc=c(dcut[2:norc],dbeta)
    if(npred==1) d2bcut=matrix(c(d2bcut[,2:norc]),npred,norc-1)
    else d2bcut=d2bcut[,2:norc]
    d2cut=d2cut[2:norc,2:norc]
    h=cbind(d2cut,t(d2bcut))
    h=rbind(h,cbind(d2bcut,d2beta))

    dif=solve(h,sc)
    mxdif=max(abs(dif))
    cutp=cutp-dif[1:(norc-1)]
    b=b-dif[norc:np]
    # modification for cutp out of order
    chk=cutp[-1]-cutp[1:(norc-2)]
    ibad=sum(chk<=0)
    while(ibad>0)
    { dif=dif/2
      mxdif=mxdif/2
      cutp=cutp+dif[1:(norc-1)]
      b=b+dif[norc:np]
      chk=cutp[-1]-cutp[1:(norc-2)]
      ibad=sum(chk<=0)
    }
    iter=iter+1
    if(iprint==1)
    { cat("iter=",iter,", (with centered x's) cutp=", cutp, ", b=",b,"\n")
      cat("         scorevec=", sc,"\n\n")
    }
  }

  if(iter>=maxiter) cat("*** did not converge, check with iprint=1\n")
  # cutpoints with original x
  for(j in 1:npred)
  { cutp=cutp-b[j]*xmn[j] }
  if(iprint==1) cat("(with original x's) cutp=", cutp,"\n")

  # Hessian with original x's, repeat of previous code with x instead of xc
  tem=x%*%b
  cutb=c(-10,cutp,10)
  ub=tem+cutb[y+1]
  lb=tem+cutb[y]
  ucdf=platent(ub)
  lcdf=platent(lb)
  updf=dlatent(ub)
  lpdf=dlatent(lb)
  nllk=0
  dbeta=rep(0,npred)
  dcut=rep(0,norc+1)
  d2beta=matrix(0,npred,npred)
  d2bcut=matrix(0,npred,norc+1)
  d2cut=matrix(0,norc+1,norc+1)
  for(i in 1:n)
  { uderi=-updf[i]*ub[i]
    lderi=-lpdf[i]*lb[i]
    xx=x[i,]
    pri=ucdf[i]-lcdf[i]
    nllk=nllk-log(pri)
    prderi=updf[i]-lpdf[i]
    dbeta=dbeta+xx*prderi/pri 
    k=y[i]
    dcut[k+1]=dcut[k+1]+updf[i]/pri
    dcut[k]=dcut[k]-lpdf[i]/pri
    pr2=pri^2
    d2beta=d2beta+ outer(xx,xx)*((uderi-lderi)*pri-prderi^2)/pr2 
    d2bcut[,k+1]=d2bcut[,k+1]+xx*(uderi*pri-updf[i]*prderi)/pr2
    d2bcut[,k]=d2bcut[,k]+xx*(-lderi*pri+lpdf[i]*prderi)/pr2
    d2cut[k+1,k+1]=d2cut[k+1,k+1]+ (uderi*pri-updf[i]^2)/pr2
    d2cut[k,k]=d2cut[k,k]+ (-lderi*pri-lpdf[i]^2)/pr2
    tem2=updf[i]*lpdf[i]/pr2
    d2cut[k,k+1]=d2cut[k,k+1]+tem2
    d2cut[k+1,k]=d2cut[k+1,k]+tem2
  }  

  sc=c(dcut[2:norc],dbeta)
  # print(sc)

  if(npred==1) d2bcut=matrix(c(d2bcut[,2:norc]),npred,norc-1)
  if(npred>1) d2bcut=d2bcut[,2:norc]
  d2cut=d2cut[2:norc,2:norc]
  h=cbind(d2cut,t(d2bcut))
  h=rbind(h,cbind(d2bcut,d2beta))
  h=-h
  #print(h)
  covm=solve(h)
  #print(covm)

  if(iprint==1)
  { cat("nllk= ", nllk,"\n")
    cat("cutpts= ", cutp,"\n")
    cat("beta= ", b,"\n")
    cat("SEs : ",sqrt(diag(covm)),"\n\n")
  }
  list(negloglik=nllk, gam=cutp, reg=b, cov=covm)
}
