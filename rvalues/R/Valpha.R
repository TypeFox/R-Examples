Valpha <- function(V, alpha.grid, smooth = "none")  {
### allow for alpha.upper < 1
### To do this, easier just to write another C version of VVcut

  ########################################################
  ## Last edited: 4/6/14
  ## Input
  ##   V             nunits x ngrid matrix containing posterior tail probs
  ##   alpha.grid    ngrid vector
  ##   smooth        smoothing parameter
  ####
  ## Output
  ##   a vector of r-values
  ########################################################
  if(ncol(V)!=length(alpha.grid)) {
    stop("The number of columns must equal the length of the alpha grid")
  }
  
  aa <- (alpha.grid > 0) & (alpha.grid < 1)
  alpha.grid <- alpha.grid[aa]
  V <- V[,aa]
  ngrid <- length(alpha.grid)

  lam <- numeric(ngrid)
  for( j in 1:ngrid )  {
      lam[j] <- quantile(V[,j], probs = 1 - alpha.grid[j], names = FALSE, type = 1)
  }
  ## Smooth \lambda_{\alpha}
  if(smooth=="none") {
      cc2 <- approxfun(alpha.grid, lam, yleft = 1, yright = 0)
      lam.smooth.eval <- cc2(alpha.grid)
      lam.smooth <- approxfun( c(0,alpha.grid,1), c(1,lam.smooth.eval,0))
  }
  else {
      cc2 <- supsmu( alpha.grid, lam, bass= smooth )
      lam.smooth.eval <- cc2$y
      lam.smooth <- approxfun( c(0,cc2$x,1), c(1,lam.smooth.eval,0))
  }
  ### For each row of Valpha, determine the index at which
  ### Valpha[i,] intersects lambda_{\alpha}
  rvalues <- VVcut(V, lam.smooth.eval, nrow(V), ngrid, alpha.grid)
 
  ans <- list()
  ans$rvalues <- rvalues
  ans$Vmarginals <- lam
  ans$Vmarginals.smooth <- lam.smooth
  return(ans)
}
