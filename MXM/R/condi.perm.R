condi.perm <- function(ind1, ind2, cs, dat, type = "pearson", rob = FALSE, R = 999) {

  x1 = dat[, ind1]
  x2 = dat[, ind2 ]
  n <- length(x1) 
  d <- sum(cs == 0)
  
  if ( d > 0 ) {  ## There are no conditioning variables

    if (rob == TRUE) { ## robust correlation
      #cont = robust::lmRob.control(mxr = 2000, mxf = 2000, mxs = 2000 )  ## only used in robust linear regression
      #mod1 <- robust::lmRob( x1 ~ x2, control = cont )
      #mod2 <- robust::lmRob( x2 ~ x1, control = cont ) 
      mod1 <- MASS::rlm( x1 ~ x2, maxit = 2000 )
      mod2 <- MASS::rlm( x2 ~ x1, maxit = 2000 ) 
      
      b1 <- coef( mod1 )[2]
      b2 <- coef( mod2 )[2]
      r <- sqrt( abs(b1 * b2) )
      stat <- abs( 0.5 * log( (1 + r) / (1 - r) ) )  ## absolute of the test statistic
      
      e1 <- resid(mod1)
      e2 <- resid(mod2)
      res <- permcor( cbind(e1, e2), R )
      stat <- res[1]
      pvalue <- res[2]

    }else{
      res <- permcor( cbind(x1, x2), R )
      stat <- res[1]
      pvalue <- res[2]
    }

  }else{  ## there are conditioning variables

   if (rob == TRUE) { ## robust correlation
      # e1 <- resid( robust::lmRob( x1 ~ ., data = as.data.frame(dat[, cs]) ), control = cont )
      # e2 <- resid( robust::lmRob( x2 ~.,  as.data.frame(dat[, cs]) ), control = cont )
      e1 <- resid( MASS::rlm( x1 ~ ., data = as.data.frame(dat[, cs]) ), maxit = 2000 )
      e2 <- resid( MASS::rlm( x2 ~.,  data = as.data.frame(dat[, cs]) ), maxit = 2000 )
      res <- permcor( cbind(e1, e2), R)
      stat <- res[1]
      pvalue <- res[2]

   }else{
      er <- resid( lm( cbind( x1, x2 ) ~ ., data = as.data.frame(dat[, cs]) ) )
      res <- permcor( er, R )
      stat <- res[1]
      pvalue <- res[2]
   }
  }
  #lets calculate the stat and p-value which are to be returned

  dof <- n - d - 3; #degrees of freedom
  stat <- stat / dof
  pvalue <- log( pvalue )
  
  result <- c(stat, pvalue, dof)
  names(result) <- c('test', 'logged.p-value', 'df')
  result
}
