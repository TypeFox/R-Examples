safe.band <-
function(x, k) 
{
  ### k = number of subdiagonals to keep
  ### x should be centered
  n = dim(x)[1]
  p = dim(x)[2]
  ## protect k to avoid singularities
  k=min(c(k, p-1));
  cholesky = diag(p)
  resid = x[,1,drop=FALSE]
  sigma2 = mean(resid^2)
  for (j in 2:p) 
  {
    if( k >= 1)
    {
      keep = (max(j-k, 1):(j-1))
      newy <- x[, j]
      newx <- resid[,keep]
      
      if( length(keep) < (n-1))
      {	
        xtx <- crossprod(newx)
        xty <- crossprod(newx,newy)			
        phi.new = qr.solve(xtx, xty, tol=1e-200)
      }else
      {
        o=svd(newx, nv=(n-1), nu=(n-1))
		d=o$d[-n]
		w=1/d
        phi.new=o$v%*%(w*crossprod(o$u, newy))
	  }
	  thisresid=newy - newx %*% phi.new
	  sigma.new=0
	  if( length(keep) < (n-1) )
        sigma.new <- mean((thisresid)^2)
      cholesky[j, keep] <- phi.new
      resid = cbind(resid, thisresid)
      sigma2 <- c(sigma2, sigma.new)
    }else
    {
      sigma2 <- c(sigma2, mean(x[,j]^2) )
    }
  }
  sigma <- tcrossprod(cholesky %*% diag(sigma2),cholesky)
  return(sigma)
}

