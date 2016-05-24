multinomial.test <- 
function(x,p=rep(1/length(x),length(x))) {
  call <- match.call()
  data.name <- if(length(call$x)==1) {call$x} else {paste(call$x[1],"(",paste(call$x[-1],collapse=","),")",sep="")}
  if (!is.vector(x)) {stop("'x' must be a vector")}
  if (sum(p)!=1) {stop("sum of probabilities must be 1")}
  if (length(x)!=length(p)) {stop("'x' and 'p' lengths differ")}
  size <- sum(x)
  groups <- length(x)
  numEvents <- choose(size+groups-1,groups-1)
  pObs <- dmultinom(x,size,p)
  findVectors <- function(groups,size) {
    if (groups==1) {
	mat <- size
    } else {
	mat <- matrix(rep(0,groups-1),nrow=1)
	for (i in 1:size) {
	  mat <- rbind(mat,findVectors(groups-1,i))
	}
	mat <- cbind(mat,size-rowSums(mat))
    }
    mat
  }
  eventMat <- findVectors(groups,size)
  eventProb <- apply(eventMat,1,function(x) dmultinom(x,size,p))
  p.val <- sum(eventProb[eventProb<=pObs])
  method <- "Exact multinomial test"
  result <- list(method=method,p.value=p.val,data.name=data.name,observed=x,expected=p*size)
  class(result) <- "htest"
  return(result)
}
