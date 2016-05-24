
'lambda.stem.ml' <- 
function(n, tb, eps=0)
{
  betaF <- function(r, t1)
  {
    xf <- (exp(r*t1)-1)/(exp(r*t1)-eps)
    xf;
  }
  LF <- function(r)
  {
   -(log(1 - betaF(r, tb)) + (n-1)*log(betaF(r, tb)))
  }
  res <- suppressWarnings(nlm(function(p) LF(p[1]), .05))  
  res$lambda <- res$estimate/(1-eps)
  res$r <- res$estimate
  res
}


