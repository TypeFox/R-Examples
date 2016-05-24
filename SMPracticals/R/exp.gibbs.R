"exp.gibbs" <-
function(u1=NULL, u2=NULL, B, I=100, S=100 )
{
  H <- function(u0, B, S=ncol(u0) )
  {
    rexp.trunc <- function( S, lambda, B) 
            -log( 1-runif( n=S, max=1-exp(-lambda*B)) )/lambda
    u1 <- rexp.trunc( S, u0[2,], B)
    u2 <- rexp.trunc( S, u1, B)
    rbind( u1, u2 )
  }
  out <- array(NA, dim = c(2, S, I))
  if (is.null(u1)) u1 <- runif(n=S, max=B)
  if (is.null(u2)) u2 <- runif(n=S, max=B)
  out[,,1] <- rbind( u1, u2 )
  for (i in 2:I)  out[,,i] <- H( out[,,i-1], B )
  out
}

