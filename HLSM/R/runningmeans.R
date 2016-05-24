###Running Means####
runningmeans <- function(x, start){
  n.iter=length(x)
  val=sapply(start:n.iter, function(ii)  mean(x[start:ii]))
  return(val)
}
