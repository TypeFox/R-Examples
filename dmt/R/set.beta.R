
# (C) 2008-2011 Leo Lahti and Olli-Pekka Huovilainen          
# All rights reserved. 
# FreeBSD License (keep this notice)     


set.beta.isotropic <- function (M, W, phi) {

  # (C) 2008-2011 Leo Lahti and Olli-Pekka Huovilainen          
  # All rights reserved. 
  # FreeBSD License (keep this notice)     

  # assuming isotropic marginal covariance
  M%*%t(W)/phi

}


set.beta.fullcov <- function (M, W, phi.inv) {

  # (C) 2008-2011 Leo Lahti and Olli-Pekka Huovilainen          
  # All rights reserved. 
  # FreeBSD License (keep this notice)     

  # assuming full marginal covariance
  # as in section 4.1 EM algorithm from BachJordan probCCA paper
  M%*%t(W)%*%phi.inv

}


