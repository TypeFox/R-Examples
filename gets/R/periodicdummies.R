periodicdummies <-
function(x, values=1)
{
  ##prepare:
  if(!is.regular(x, strict=TRUE)) stop("Vector/matrix not strictly regular")
  iFreq <- frequency(x)
  if(iFreq==1) stop("Frequency must be greater than 1")
  if(!is.zoo(x)){ x <- as.zooreg(x) }
  vCycle <- as.numeric(cycle(x))

  ##values argument:
  if(length(values)==1){ values <- rep(1,iFreq) }
  if(length(values)!=iFreq) stop("length(values) must be 1 or equal to frequency")

  ##make dummies:
  mDums <- matrix(0,NROW(x),iFreq)
  colnames(mDums) <- paste("dum", 1:iFreq, sep="")
  for(i in 1:NCOL(mDums)){
    whereIs <- which(vCycle==i)
    mDums[whereIs,i] <- values[i]
  }

  ##out:
  mDums <- zoo(mDums, order.by=index(x), frequency=iFreq)
  return(mDums)
}
