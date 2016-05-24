"propose.coeff.gp" <-
function(object,block=0,proposal.sd=1,...){
  # proposes coefficients as a normal proposal with the proposal std dev the product of proposal.sd and the prior std dev
  # if block=0, all coefficients are proposed
  m1=object$gridsize[1]
  m2=object$gridsize[2]
  if(block){  # block is specified
    if(is.null(object$blocks)){
      stop("No blocks have been specified in the gp object")
    }
    if(block>object$num.blocks || block<1){
      stop("The block number supplied does not exist.")
    }
    elements=object$blocks[[block]]  # subset of coefficients being proposed
    if(object$d==2){
      n.elements=dim(elements)[1]
    } else{
      n.elements=length(elements)
    }
    proposal.sd=proposal.sd*sqrt(object$variances[elements])  # proposal std devs
    object$coeff[elements]=object$coeff[elements]+rnorm(n.elements,0,proposal.sd)+1i*rnorm(n.elements,0,proposal.sd)  # propose the coefficients
  } else{  # all coefficients are proposed
    proposal.sd=proposal.sd*sqrt(object$variances)
    object$coeff=object$coeff+rnorm(m1*m2,0,proposal.sd)+1i*rnorm(m1*m2,0,proposal.sd)
  }
  # deterministically calculates some coefficients so processes are real-valued
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
  updateprocess(object) # recalculate process values based on new coefficients
  return(NULL)
}
