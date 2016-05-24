# Version: 30-11-2012, Daniel Fischer

mwSingleNullDist <- function(X,g,t,goi,nper){
  permTestValues <- c(rep(0,nper))
  
  for(i in 1:nper){
	permgroups <- sample(g)
	permTestValues[i] <- oneVsOther(X,permgroups,t,goi)
  }
    return(permTestValues)

}