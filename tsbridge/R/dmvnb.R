dmvnb <-
function(bug, sims, MU, COV, P){ 
  if(is.null(colnames(sims)))
    stop("columns of sims can not be NULL. names should correspond to parameters")
  if(class(bug)!="tsbugs")
   stop("bug must be a object of class tsbugs")
  nn<-paste0("delta.", 1:bug$info$n)
  nn<-nn[nn %in% colnames(sims)]
  nn<-colnames(sims)[!(colnames(sims) %in% nn)]
  mvn <- dmvnorm(sims[,nn], mean=MU[nn], sigma=COV[nn,nn], log=TRUE)
  delta<-NULL
  delta<-theta.it(bug, sims)$delta
  mvb <- dbinom(delta, 1, prob=P, log=TRUE)
  mvnb <- mvn + rowSums(mvb)
  return(mvnb)
}
