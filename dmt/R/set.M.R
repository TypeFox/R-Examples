# (C) 2008-2011 Leo Lahti and Olli-Pekka Huovilainen          
# All rights reserved. 
# FreeBSD License (keep this notice)     



# "Do not worry about your difficulties in Mathematics. I can assure you
# mine are still greater."
# - Albert Einstein


set.M.isotropic <- function (wtw, sigma.sq, dx) {

  # (C) 2008-2011 Leo Lahti and Olli-Pekka Huovilainen          
  # All rights reserved. 
  # FreeBSD License (keep this notice)     

  solve(wtw + diag(sigma.sq, dx, dx)) 
  
}

set.M <- function (W, phi) {

  # (C) 2008-2011 Leo Lahti and Olli-Pekka Huovilainen          
  # All rights reserved. 
  # FreeBSD License (keep this notice)     

  solve(t(W)%*%W/phi + diag(ncol(W)))

}


set.M.full <- function (W, phi.inv) {#

  # (C) 2008-2011 Leo Lahti and Olli-Pekka Huovilainen          
#  # All rights reserved. 
  # FreeBSD License (keep this notice)     #

  # for full marginal covariance
  solve(t(W)%*%phi.inv%*%W + diag(ncol(W)))
}


set.M.full2 <- function (W, phi.inv) {

  # (C) 2008-2011 Leo Lahti and Olli-Pekka Huovilainen          
  # All rights reserved. 
  # FreeBSD License (keep this notice)     

  # modified from G in Bishop's book
  # when phi$total is block-diagonal, we can use sums of the two blocks

  # This corresponds to 
  # M <- set.M.full(W$total, phi.inv$total) # for non-matched case
  # but should be faster in general

  # for full marginal covariance
  solve(t(W$X)%*%phi.inv$X%*%W$X + t(W$Y)%*%phi.inv$Y%*%W$Y + diag(ncol(W$X)))

}
 
