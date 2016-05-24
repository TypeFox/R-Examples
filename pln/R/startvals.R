startalphas<- function (x, ncat, nitem=NULL) {

  ## input checking & data prep
  myInput<-check.input(x, ncat, nitem)
  
  out <- startpln.func(myInput$nitem, myInput$ncat, myInput$nrec, myInput$myX)

  return(out$alphas)

}

startpln.func<-function(nitem, ncat, nrec, myX){

  ## prep return variables
  alphas=rep(0,nitem*(ncat-1))

  ## TO DO: add PACKAGE argument here
  out <- .C("Rstartpln",
            as.integer(nitem), as.integer(ncat), as.integer(nrec), as.double(myX),
            alphas=as.double(alphas))
  return(out)
}

startbetas<-function(x, ncat, nitem=NULL){
  myInput<-check.input(x, ncat, nitem)
  x<-myInput$myX
  betas<-startbetas.func(x)
  return(betas)
}

startbetas.func<-function(x){
  nn<-ncol(x)
  nitem<-nn-1
  y<-x[,-nn]
  fr<-x[,nn]
  tot<-sum(fr)
  mnvec<-apply(y*fr,2,sum)/tot ## means
  vvec<-apply(y*y*fr,2,sum)/tot ## variances
  vvec<-vvec-mnvec^2
  cc<-matrix(0,nitem,nitem)
  for(j in 1:(nitem-1))
  { for(k in (j+1):nitem)
    { ss<-sum(y[,j]*y[,k]*fr)/tot
      den<-vvec[j]*vvec[k]
      cc[j,k]<-(ss-mnvec[j]*mnvec[k])/sqrt(den)
   	  cc[k,j]<-cc[j,k]
    }
  }
  avcc<-apply(cc,2,sum)/(nitem-1)
  #print(avcc)
  bvec<-avcc*10; bvec[bvec>=1]=1; bvec[bvec<= -.2]= -.2
  #print(bvec)
  bvec
}