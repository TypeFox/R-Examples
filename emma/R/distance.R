distance  <- function(xpop,xspace,yspace,weight,opt)
{
#  library(clusterSim)
  yspace.norm <- as.matrix(data.Normalization(yspace,type="n4"))

  obj <- ifelse(opt=="mn",0,1)

  yspace.weight <- matrix(0,nrow(yspace),ncol(yspace))
  rownames(yspace.weight) <- rownames(yspace)
  for(i in 1:ncol(yspace.norm)){
    yspace.weight[,i] <- weight[i]*(yspace.norm[,i]-rep(obj[i],nrow(yspace.norm)))
  }

  fit <- apply(yspace.weight^2,1,sum)^(1/2)

  ## non normalized target response values
  obj.nn <- rep(NA,ncol(yspace))
  for(i in 1:ncol(yspace)){
    obj.nn[i] <- yspace[which.min(as.matrix(yspace[,i])),i]
    names(obj.nn)[i] <- which.min(as.matrix(yspace[,i]))
  }

  d <- list(fit=fit,obj.nn=obj.nn)
  return(d)
}

