nrbcpln<- function (x, ncat, nitem=NULL, alphas=NULL, betas=NULL, abound=c(-10,10),
    bbound=c(-1,10), nq=48, mxiter=200, se=TRUE, iprint=FALSE) {

  myInput<-check.input(x, ncat, nitem, nq, mxiter, iprint)

  ## get starting values if not present already
  if(!check.alphas(alphas, myInput$nitem, myInput$ncat)){
    alphas<-startpln.func(myInput$nitem, myInput$ncat, myInput$nrec, myInput$myX)$alphas
  }

  ## prep betas
  if(!check.betas(betas, myInput$nitem)){
    betas<-startbetas.func(myInput$myX)
  }

  ## check bounds
  abound<-check.bounds(alphas, abound)
  bbound<-check.bounds(betas, bbound)
  
  nrbcplnout <- nrbcpln.func(myInput$ncat, myInput$nitem, myInput$nrec, myInput$myX, alphas,
        betas, abound, bbound, myInput$nq, myInput$mxiter, myInput$iprint)

  alphas<-nrbcplnout$bcpln[1:((myInput$ncat-1)*myInput$nitem)]
  betas<-nrbcplnout$bcpln[((myInput$ncat-1)*myInput$nitem+1):(myInput$ncat*myInput$nitem)]

  out<-list(alphas=alphas,betas=betas,nllk=nrbcplnout$nrbcpln,conv=nrbcplnout$iconv)
    
  if(se){
    nrbcplncov<-nrbcplncov.func(myInput$nitem, myInput$ncat, myInput$nrec, alphas, betas,
      myInput$N, myInput$nq, myInput$iprint)
      
    V<-matrix(nrbcplncov$V, nrow=myInput$nitem*myInput$ncat,ncol=myInput$nitem*myInput$ncat)
    seVec<-sqrt(diag(V))
      
    sealphas<-seVec[1:((myInput$ncat-1)*myInput$nitem)]
    sebetas<-seVec[((myInput$ncat-1)*myInput$nitem+1):(myInput$ncat*myInput$nitem)]
      
    out<-append(out,list(sealphas=sealphas, sebetas=sebetas, vcov=V)) 
  }
    
  return(out)
 
}

nrbcpln.func<-function(ncat, nitem, nrec, myX, alphas, betas, abound, bbound, nq,
    mxiter, iprint){
  nrbcplnout<-0
  iconv<-0
  np<-ncat*nitem
  bcplnout<-rep(0,np)
  out <- .C("Rnrbcpln",
    as.integer(nitem), as.integer(ncat), as.integer(nrec), as.double(myX),
    as.double(alphas), as.double(betas), as.double(abound), as.double(bbound),
    nbcplnout=as.double(nrbcplnout), bcplnout=as.double(bcplnout), as.integer(nq),
    as.integer(mxiter), iconv=as.integer(iconv), as.integer(iprint))
  list(nrbcpln=out$nbcplnout,bcpln=out$bcplnout, iconv=out$iconv)
}

## asymptotic covariance matrix
nrbcplncov.func<-function(nitem, ncat, nrec, alphas, betas, N, nq, iprint) {
  V<-rep(0,nitem*ncat*nitem*ncat)
  params<-c(alphas,betas)
  out<-.C("Rbclcov",
          as.integer(nitem), as.integer(ncat), as.integer(N),
          as.double(params), V=as.double(V), as.integer(nq), as.integer(iprint))
  return(out)
}
