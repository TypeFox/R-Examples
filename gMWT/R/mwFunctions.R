# Version: 30-11-2012, Daniel Fischer

mwNullDist <- function(x,y,nper){
  permTestValues <- c(rep(0,nper))
  
  for(i in 1:nper){
	permvalues <- permObs2.C(x,y)
	permTestValues[i] <- getP.Cnaive(permvalues$x,permvalues$y)
  }
    return(permTestValues)

}