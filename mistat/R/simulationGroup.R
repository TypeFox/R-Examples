simulationGroup <- function(x, n){
  if(inherits(x, what="mistatSimulation", which=FALSE)){
    res <- .simulationGroup(x, n)
    return(res)
  } else
    stop("x must inherits from class mistatSimulation")
}

.simulationGroup <- function(x, n){
  group <- rep(
    ceiling(nrow(x)/n):1, 
    each=n)[nrow(x):1]
  res <- data.frame(x, group=factor(group))
  n2 <- aggregate(res$group, by=res["group"], FUN=length)
  if(n2$x[1] < n){
    res <- res[-(1:n2$x[1]),]
    warning(paste("discarded first", n2$x[1],
                  "observations because nrow(x) is not a multiple of n"))
  }
  return(res)
}