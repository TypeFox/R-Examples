NullModel <-
function(x,type=c("mle","deseq","quantile")){
  type <- match.arg(type)
  # x MUST be a n times p matrix - i.e. observations on the rows and features on the columns
  rowsum <- rowSums(x)
  colsum <- colSums(x)
  mle <- outer(rowsum, colsum, "*")/sum(rowsum)
  if(type=="mle"){
    return(list(n=mle, sizes=rowSums(x)/sum(x)))
  } else if (type=="quantile"){
    #This is quantile normalization idea of Bullard et al 2010 -- quantile-normalize using 75th quantile of observed counts for each sample, excluding zero-counts
    sample.qts <- apply(x,1,quantile,.75)
    sample.qts <- pmax(sample.qts, 1) # Don't wait to standardize by 0... min allowed is 1
    sample.qts <- sample.qts/sum(sample.qts)
    fit <- outer(sample.qts,colsum,"*")
    return(list(n=fit, sizes=sample.qts))
  } else if (type=="deseq"){ #Trying to implement idea from Anders and Huber Genome Biology 2010 paper.
    # I stole the next 3 lines from the DESeq bioconductor package.... it was in the function estimateSizeFactorsForMatrix
    counts <- t(x)
    geomeans <- exp(rowMeans(log(counts)))
    sizes <- apply(counts, 2, function(cnts) median((cnts/geomeans)[geomeans > 0]))
    rawsizestr <- sizes
    sizes <- sizes/sum(sizes)
    fit <- outer(sizes, colsum, "*")
    return(list(n=fit,sizes=sizes,geomeans=geomeans,rawsizestr=rawsizestr))
  } 
}

