"add.blocks.gp" <-
function(object,breaks=NULL,...){
  # adds block structure for sampling coefficients by block in an MCMC framework
  # if breaks is supplied, it should be a vector of the largest allowed frequency for the blocks; largest frequency should be one half the largest of the grid dimensions to ensure that all coefficients are sampled
  # return value is the gp with a new component, a list giving the indices of elements of each block
  object$blocks=list()
  d=object$d
  if(d==1){
    indices=matrix(1:object$gridsize[1],ncol=1)
    maxvals=abs(object$omega)
  } else{ # d=2
    indices=as.matrix(expand.grid(1:object$gridsize[1],1:object$gridsize[2]))
    maxvals=apply(cbind(abs(object$omega[,1]),abs(object$omega[,2])),1,max)
  }
  if(is.null(breaks)){
    object$num.blocks=log(max(object$gridsize),2) 
    breaks=2^(seq(0,object$num.blocks-1,by=1))  # default breaks, another possibility is 2^(0,1,2,4,...)-1
  }
  object$num.blocks=length(breaks)
  if(d==1){
    object$blocks[[1]]=indices[maxvals<=breaks[1] & object$omega>=0]
    for(b in 2:object$num.blocks){
      object$blocks[[b]]=indices[maxvals<=breaks[b] & maxvals>breaks[b-1] & object$omega>=0]
    }
  } else{
    object$blocks[[1]]=indices[maxvals<=breaks[1] & object$omega[,2]>=0,]
    if(!is.matrix(object$blocks[[1]])){
      object$blocks[[1]]=matrix(object$blocks[[1]],ncol=2)
    }
    for(b in 2:object$num.blocks){
      object$blocks[[b]]=indices[maxvals<=breaks[b] & maxvals>breaks[b-1] & object$omega[,2]>=0,]
      if(!is.matrix(object$blocks[[b]])){
        object$blocks[[b]]=matrix(object$blocks[[b]],ncol=2)
      }
    }
    object$blockSizes=NULL
    for(block in 1:object$num.blocks){
      if(d==2){
        object$blockSizes[block]=length(object$blocks[[block]][,1])
      } else{
        object$blockSizes[block]=length(object$blocks[[block]])
      }
    }
  }
  return(NULL)
}
