"Gibbs.sample.coeff.gp" <-
function(object,z,sig2e,meanVal=0,sdVal=1,returnHastings=FALSE,...){
  # takes a Gibbs sample, following the approach of Wikle (2002)
  
  m1=object$gridsize[1]
  m2=object$gridsize[2]

  sig2e.precMatrix=matrix(2*sdVal*sdVal/sig2e,nrow=m1,ncol=m2) 
  sig2e.precMatrix[1,1]=(1/2)*sig2e.precMatrix[1,1]  
  sig2e.precMatrix[(m1/2+1),1]=(1/2)*sig2e.precMatrix[(m1/2+1),1]
  if(object$d==2){
    sig2e.precMatrix[(m1/2+1),(m2/2+1)]=(1/2)*sig2e.precMatrix[(m1/2+1),(m2/2+1)]
    sig2e.precMatrix[1,(m2/2+1)]=(1/2)*sig2e.precMatrix[1,(m2/2+1)]
  }
  coeff.var=1/(1/object$variances+sig2e.precMatrix)
  coeff.mean=coeff.var*sig2e.precMatrix*fft(matrix((z-meanVal)/sdVal,nrow=m1,ncol=m2,byrow=FALSE), inverse = FALSE)/(sqrt(m1*m2)) # division by sqrt(m1*m2) ensures proper scaling
  object$coeff=matrix(rnorm(m1*m2,Re(coeff.mean),sqrt(c(coeff.var))),nrow=m1,ncol=m2)+(0+1i)*matrix(rnorm(m1*m2,Im(coeff.mean),sqrt(c(coeff.var))),nrow=m1,ncol=m2)
  if(object$const.fixed){
    object$coeff[1,1]=0
  } else{
    object$coeff[1,1]=Re(object$coeff[1,1])
  }
  object$coeff[(m1/2+1),1]=Re(object$coeff[(m1/2+1),1])
  object$coeff[m1:(m1/2+2),1]=Conj(object$coeff[2:(m1/2),1])
  if(object$d==2){
    object$coeff[1,(m2/2+1)]=Re(object$coeff[1,(m2/2+1)])
    object$coeff[(m1/2+1),(m2/2+1)]=Re(object$coeff[(m1/2+1),(m2/2+1)])
    object$coeff[1,m2:(m2/2+2)]=Conj(object$coeff[1,2:(m2/2)])
    object$coeff[m1:(m1/2+2),m2:(m2/2+1)]=Conj(object$coeff[2:(m1/2),2:(m2/2+1)])
    object$coeff[(m1/2+1):2,m2:(m2/2+2)]=Conj(object$coeff[(m1/2+1):m1,2:(m2/2)])
  }
  updateprocess(object)
  if(!returnHastings){
    return(NULL)
  } else{
    # this block of code determines which coefficients are actually proposed and not just determined as complex conjugates
    screenr=matrix(1,nrow=m1,ncol=m2)
    screenr[m1:(m1/2+2),1]=0
    if(object$d==2){
      screenr[(m1/2+2):m1,(m2/2+1):m2]=0
      screenr[1:(m1/2+1),(m2/2+2):m2]=0
    }
    if(object$const.fixed){
      screenr[1,1]=0
    }
    screeni=screenr
    screeni[1,1]=0
    screeni[(m1/2+1),1]=0
    if(object$d==2){
      screeni[1,(m2/2+1)]=0
      screeni[(m1/2+1),(m2/2+1)]=0
    }
    tmpr=c(dnorm(Re(object$coeff),Re(coeff.mean),sqrt(coeff.var),log=T))
    tmpi=c(dnorm(Im(object$coeff),Im(coeff.mean),sqrt(coeff.var),log=T))
    return(sum(tmpr[c(screenr==1)])+sum(tmpi[c(screeni==1)]))  # calculate logdensity of proposal
  }
}

