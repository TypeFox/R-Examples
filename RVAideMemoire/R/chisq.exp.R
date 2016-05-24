chisq.exp <- function(data,p,graph=FALSE){
  if (length(p)!=nrow(data)){stop("number of expected proportions and populations differ")}
  n <- integer(nrow(data))
  for (i in 1:nrow(data)) {n[i] <- sum(data[i,])}
  n.theo1 <- n*p
  n.theo2 <- n*(1-p)
  n.theo.mat <- matrix(c(n.theo1,n.theo2),nrow=nrow(data),dimnames=list(rownames(data),colnames(data)))
  cochran.max <- ceiling(0.8*length(data))
  cochran.min <- length(data)-cochran.max
  result <- list(p.theo=p,mat=n.theo.mat,cochran=cochran.min)
  class(result) <- "chisq.exp"
  if (graph) {mosaicplot(t(n.theo.mat),main="Expected distribution",col=TRUE)}
  return(result)
}

print.chisq.exp <-
function(x,...) {
  cat("\nExpected counts\n")
  print(x$mat)
  cat(paste("\nCochran's rule: maximum",x$cochran,"count(s) can be < 5\n\n"))
}
