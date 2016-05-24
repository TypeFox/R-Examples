generate  <- function(in.name,nlev,lower,upper)
{
  X <- list()

  for(i in 1:length(in.name)){
    X[[i]] <- seq(lower[i],upper[i],length=nlev[i])
  }
  names(X) <- in.name
  expe <- expand.grid(X)
#  expe2 <- expe[round(expe[,2]+expe[,3]+expe[,6],3)==0.005 & expe[,2]>=0.00225,]

  return(data.frame(expe))
}

