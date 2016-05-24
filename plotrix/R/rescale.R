# linearly transforms a vector or matrix of numbers to a new range

rescale<-function(x,newrange) {
 if(missing(x) | missing(newrange)) {
  usage.string<-paste("Usage: rescale(x,newrange)\n",
   "\twhere x is a numeric object and newrange is the new min and max\n",
   sep="",collapse="")
  stop(usage.string)
 }
 if(is.numeric(x) && is.numeric(newrange)) {
  xna<-is.na(x)
  if(all(xna)) return(x)
  if(any(xna)) xrange<-range(x[!xna])
  else xrange<-range(x)
  # if x is constant, just return it
  if(xrange[1] == xrange[2]) return(x)
  mfac<-(newrange[2]-newrange[1])/(xrange[2]-xrange[1])
  return(newrange[1]+(x-xrange[1])*mfac)
 }
 else {
  warning("Only numeric objects can be rescaled")
  return(x)
 }
}
