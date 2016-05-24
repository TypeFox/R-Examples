`mu.rank.nna` <-
function(x) 
{
  if (sum(dim(x)>1)< 2) {
    if ((n<-length(x))>0) {
      sx <- x[sl <- sort.list(x)]                          # sorted !NA
      ntie <- diff(c(0, (1:n)[c(sx[-1] != sx[-n], TRUE)])) # fast.rle
      x[sl] <- rep(cumsum(ntie) - (ntie-1)/2, ntie)# ranks
    }
    return(x)
  } else {
    lenx <- length(x)
    uxx <- mu.Sums(mu.GE(x))
    rxx <-(uxx[[1]]+(lenx+1))/2
    return(rxx)
  }
}
