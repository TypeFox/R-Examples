
# (C) 2008-2011 Leo Lahti and Olli-Pekka Huovilainen          
# All rights reserved. 
# FreeBSD License (keep this notice)     

# "I am among those who think that science has great beauty. A scientist
# in his laboratory is not only a technician: he is also a child placed
# before natural phenomena which impress him like a fairy tale."
# - Marie Curie


z.expectation <- function (model, X, Y = NULL) {

  # (C) 2008-2011 Leo Lahti and Olli-Pekka Huovilainen          
  # All rights reserved. 
  # FreeBSD License (keep this notice)     

    W <- getW(model)
  phi <- getPhi(model)

  # Center data
  X <- t(centerData(t(X), rm.na = TRUE))
  
  if (!is.null(Y)){

    Y <- t(centerData(t(Y), rm.na = TRUE))

    phi.inv <- list(X = solve(phi$X), Y = solve(phi$Y))
    S <- set.M.full2(W, phi.inv)
    #S <- solve(t(W$X)%*%solve(phi$X)%*%W$X + t(W$Y)%*%solve(phi$Y)%*%W$Y + diag(ncol(W$X)))

    return( S%*%(t(W$X)%*%phi.inv$X%*%X + t(W$Y)%*%phi.inv$Y%*%Y) )
  } else {
    phi.inv <- solve(phi$total)
    S <- set.M.full(W$total, phi.inv)
    return(S %*% (t(W$total) %*% phi.inv %*% X))
  }
}