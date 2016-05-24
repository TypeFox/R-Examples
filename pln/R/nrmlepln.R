nrmlepln <- function (x, ncat, nitem=NULL, alphas=NULL, betas=NULL, abound=c(-10,10),
    bbound=c(-1,10), nq=48, mxiter=200, m2=TRUE, iprint=FALSE) {

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

  nrmleplnout <- nrmlepln.func(myInput$nitem, myInput$ncat, myInput$nrec, myInput$myX, alphas,
        betas, abound, bbound, myInput$nq, myInput$mxiter, myInput$iprint)
  alphas<-nrmleplnout$mlePlnOut[1:((myInput$ncat-1)*myInput$nitem)]
  betas<-nrmleplnout$mlePlnOut[((myInput$ncat-1)*myInput$nitem+1):(myInput$ncat*myInput$nitem)]
  
  V<-matrix(nrmleplnout$invHesOut, nrow=myInput$nitem*myInput$ncat,ncol=myInput$nitem*myInput$ncat)
  seVec<-nrmleplnout$seVecOut
    
  sealphas<-seVec[1:((myInput$ncat-1)*myInput$nitem)]
  sebetas<-seVec[((myInput$ncat-1)*myInput$nitem+1):(myInput$ncat*myInput$nitem)]
    
    
  out<-list(alphas=alphas,betas=betas,nllk=nrmleplnout$nllkOut,sealphas=sealphas,
    sebetas=sebetas,invhes=V)

  if(m2){
    m2out<-m2pln.func(myInput$nitem, myInput$ncat, myInput$nrec, myInput$myX, alphas, betas,
      myInput$nq, myInput$iprint)
    pval<-pchisq(m2out$m2stat, m2out$df, lower.tail=FALSE)
      ##out<-append(out, list(samplemout=m2out$samplem, m2stat=m2out$m2stat, df=m2out$df, pval=pval))      
    out<-append(out, list(teststat=m2out$m2stat, df=m2out$df, pval=pval))
  }

  return(out)
}

nrmlepln.func <- function(nitem, ncat, nrec, myX, alphas, betas, abound, bbound, nq, mxiter,
    iprint){
  
  ## prep return variables
  nllkOut<-0
  iconv<-0
  np<-ncat*nitem
  mlePlnOut<-rep(0,np)
  seVecOut<-rep(0,np)
  invHesOut<-matrix(0,nrow=np,ncol=np)

  ## TO DO: add PACKAGE argument here
  out <- .C("Rnrmlepln",
            as.integer(nitem), as.integer(ncat), as.integer(nrec), as.double(myX),
            as.double(alphas), as.double(betas), as.double(abound), as.double(bbound),
            nllkOut=as.double(nllkOut), mlePlnOut=as.double(mlePlnOut),
            seVecOut=as.double(seVecOut), invHesOut=as.double(invHesOut),
            as.integer(nq), as.integer(mxiter), iconv=as.integer(iconv),
            as.integer(iprint))
  return(out)
  
}
