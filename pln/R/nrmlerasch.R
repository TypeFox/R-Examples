nrmlerasch <- function (x, ncat, nitem=NULL, alphas=NULL, abound=c(-10,10),
    bbound=c(-1,10), nq=48, mxiter=200,  m2=TRUE, iprint=FALSE) {

  myInput<-check.input(x, ncat, nitem, nq, mxiter, iprint)

  ## get starting values if not present already
  if(!check.alphas(alphas, myInput$nitem, myInput$ncat)){
    alphas<-startpln.func(myInput$nitem, myInput$ncat, myInput$nrec, myInput$myX)$alphas
  }

  ## prep betas
#  if(is.null(betas) || betas%%1!=0 || length(betas)>1){
    betas<-1
#  }

  ## check bounds
  abound<-check.bounds(alphas, abound)
  bbound<-check.bounds(betas, bbound)      
        
  nrmleraschout <- nrmlerasch.func(myInput$nitem, myInput$ncat,
        myInput$nrec, myInput$myX, alphas, betas, abound, bbound,
        myInput$nq, myInput$mxiter, myInput$iprint)

  alphas<-nrmleraschout$mlePlnOut[1:((myInput$ncat-1)*myInput$nitem)]
  betas<-nrmleraschout$mlePlnOut[((myInput$ncat-1)*myInput$nitem)+1]
  
  V<-matrix(nrmleraschout$invHesOut, nrow=((myInput$ncat-1)*myInput$nitem)+1,
        ncol=((myInput$ncat-1)*myInput$nitem)+1)
  seVec<-nrmleraschout$seVecOut
    
  sealphas<-seVec[1:((myInput$ncat-1)*myInput$nitem)]
  sebetas<-seVec[((myInput$ncat-1)*myInput$nitem)+1]
    
    
  out<-list(alphas=alphas,betas=betas,nllk=nrmleraschout$nllkOut,sealphas=sealphas,
    sebetas=sebetas,invhes=V)

  if(m2){
    m2out<-m2rasch.func(myInput$nitem, myInput$ncat, myInput$nrec, myInput$myX,
      alphas, betas, myInput$nq, myInput$iprint)
    pval<-pchisq(m2out$m2stat, m2out$df, lower.tail=FALSE)
      ##out<-append(out, list(samplemout=m2out$samplem, m2stat=m2out$m2stat, df=m2out$df,pval=pval))
    out<-append(out, list(teststat=m2out$m2stat, df=m2out$df,pval=pval))
  }

  return(out)
  
}

nrmlerasch.func <- function(nitem, ncat, nrec, myX, alphas, betas, abound, bbound,
    nq, mxiter,iprint){
  
  ## prep return variables
  nllkOut<-0
  iconv<-0
  np<-ncat*nitem-nitem+1
  mlePlnOut<-rep(0,np)
  seVecOut<-rep(0,np)
  invHesOut<-matrix(0,nrow=np,ncol=np)

  ## TO DO: add PACKAGE argument here
  ## TO DO: coerce invHes into a matrix?
  out <- .C("Rnrmlerasch",
            as.integer(nitem), as.integer(ncat), as.integer(nrec), as.double(myX),
            as.double(alphas), as.double(betas), as.double(abound),
            as.double(bbound),nllkOut=as.double(nllkOut),
            mlePlnOut=as.double(mlePlnOut), seVecOut=as.double(seVecOut),
            invHesOut=as.double(invHesOut), as.integer(nq), as.integer(mxiter),
            iconv=as.integer(iconv), as.integer(iprint))
  return(out)
  
}
