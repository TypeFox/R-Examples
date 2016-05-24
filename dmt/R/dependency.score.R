# (C) 2008-2012 Leo Lahti and Olli-Pekka Huovilainen          
# All rights reserved. 
# FreeBSD License (keep this notice)     


# "The important thing in science is not so much to obtain new facts as
#  to discover new ways of thinking about them."
#  - Sir William Bragg 



dependency.score <- function ( model ) {

  # (C) 2008-2012 Leo Lahti and Olli-Pekka Huovilainen          
  # All rights reserved. 
  # FreeBSD License (keep this notice)     

  W   <- model$W
  phi <- model$phi

  if (!is.null(W$X)){

    # this equals to the trace of the full Phi
    noise <- sum(diag(as.matrix(phi$X))) + sum(diag(as.matrix(phi$Y))) 
    W$total <- rbind(W$X, W$Y)

  } else {

    # For single data case, check proportion between 
    # latent covariance and expected diagonal noise
    noise <- sum(diag(as.matrix(phi$total)))

  }

  signal <- sum( diag(W$total %*% t(W$total)) ) # trace of full WWt covariance 

  return( signal/noise )

}

