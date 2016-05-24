# Function to sample a VAR sequence

simulateVAR <- function(
  pars, # a single VAR matrix or a list of VAR matrices for each lag
  lags = 1, # A sequence of lags that corresponds to the lags of each matrix in the 'pars' argument
  Nt, # Number of time points to sample
  init, # Initial state, defaults to zeros
  errorSD = 0.1
)
{
  if (is.matrix(pars)) pars <- list(pars)
  if (any(sapply(pars,function(x) length(unique(dim(x))) > 1 ))) stop ("non-square graph detected.")
  
  Ni <- ncol(pars[[1]])
  maxLag <- max(lags)
  
  if (missing(init)) 
  {
    init <- matrix(0,maxLag,Ni)
  }
  

    Res <- matrix(,Nt,Ni)
    Res[1:maxLag,] <- init
    for (t in (maxLag+1):(Nt))
    {
      Res[t,] <- rowSums(do.call(cbind,lapply(seq_along(lags),function(i)pars[[i]] %*% Res[t-lags[i],]))) + rnorm(Ni, 0, errorSD)
    }    

  return(as.data.frame(Res))
}