rating.prob <- function(x,g,type=c("prob","cumprob","class")) {
  if (!is.ordered(x)) {stop(paste0(deparse(substitute(x))," is not an ordered factor"))}
  if (!is.factor(g)) {g <- factor(g)}
  x <- x[complete.cases(x,g)]
  g <- g[complete.cases(x,g)]
  type <- match.arg(type)
  if (type=="prob") {
    reps <- tapply(x,g,table)
    result <- matrix(0,nrow=nlevels(g),ncol=nlevels(x),dimnames=list(levels(g),levels(x)))
    for (i in 1:nrow(result)) {
	result[i,] <- reps[[i]]/sum(reps[[i]])
    }
  } else if (type=="cumprob") {
    reps <- tapply(x,g,table)
    result <- matrix(0,nrow=nlevels(g),ncol=nlevels(x),dimnames=list(levels(g),levels(x)))
    for (i in 1:nrow(result)) {
	result[i,] <- cumsum(reps[[i]])/sum(reps[[i]])
    }
  } else if (type=="class") {
    result <- tapply(x,g,mod)
  }
  return(result)
}
