weights_gausscopula <-
function(data,prior_type = "uniform",a=1,b=1,nbcores = 1){
  p <- ncol(data)
  n <- nrow(data)
  data.cdf <- apply(data,2,function(x){ecdf(x)(x) - 1/(2*n)})
  if (prior_type == "beta"){
    prior.cor <- function(rho){0.5*dbeta(0.5*(rho+1),a,b)}      
  } else if (prior_type == "uniform"){
    prior.cor <- function(rho){0.5}    
  } else {
    message("Unknown prior. Uniform used instead.",appendLF=FALSE)
    prior.cor <- function(rho){0.5}
  }
  uptri <- upper.tri(matrix(0,p,p))
  
  if (requireNamespace("parallel",quietly = TRUE) &&
        (nbcores > 1)){
    W.res <- matrix(parallel::mcmapply(function(y, i, j) if (y){integrate(function(x) gauss.copula(x,data.cdf[,i],data.cdf[,j])*prior.cor(x),-1,1)$value} else 0,
                             uptri, 
                             row(uptri), 
                             col(uptri),
                             mc.cores = nbcores), 
                    nrow = nrow(uptri)) 
  } else {
    W.res <- matrix(mapply(function(y, i, j) if (y){integrate(function(x) gauss.copula(x,data.cdf[,i],data.cdf[,j])*prior.cor(x),-1,1)$value} else 0,
                           uptri, 
                           row(uptri), 
                           col(uptri)), 
                    nrow = nrow(uptri))
  }
  
  return(W.res + t(W.res))
}
