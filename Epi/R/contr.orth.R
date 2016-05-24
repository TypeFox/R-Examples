contr.orth <-
  function( n )
  {
  if( is.numeric( n ) && length( n )==1 )
  levs <- 1:n
  else {
  levs <- n
  n <- length( n )
  }
  Z <- contr.sum( n )
  L <- 1:n - mean(1:n)
  contr <- Z - L%*%( ( t(L) %*% L )^(-1) ) %*% ( t(L) %*% Z )
  contr[,1:(n-2)]
  }
