NullModelTest <-
function(null.out,x,xte,type=c("mle","deseq","quantile")){
  type <- match.arg(type)
  if(type=="mle"){
    sizeste <- rowSums(xte)/sum(x)
    nste <- outer(sizeste, colSums(x), "*")
  } else if (type=="quantile"){
    sizeste <- pmax(1, apply(xte,1,quantile,.75)) # don't want to scale by anything less than 1...
    sizeste <- sizeste/sum(apply(x, 1, quantile, .75))
    nste <- outer(sizeste, colSums(x), "*")
  } else if (type=="deseq"){
    countste <- t(xte)
    sizeste <- apply(countste,2,function(cnts) median((cnts/null.out$geomeans)[null.out$geomeans>0], na.rm=TRUE))/sum(null.out$rawsizestr)
    nste <- outer(sizeste, colSums(x), "*")
  }
  return(list(nste=nste,sizeste=sizeste))
}

