`propdiff.mblwoc` <-
function(len,c1,d1,c2,d2,level=0.95)
{
  #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  #**************************************************************************
  # Define initial step
  # (as a function of any frequentist sample size estimate)

  step <- ceiling(log(propdiff.freq0(len,c1,d1,c2,d2,level))/log(2))

  # Also define the threshold to cross for the quantity under study (the
  # length or the coverage probability of an HPD region)

  threshold <- level

  # and define a factor, which is +/- 1, depending on if the quantity under
  # study is (almost) surely too large or too small when making no
  # observations [-1 if the quantity to measure is DEcreasing with n
  #               +1 if the quantity to measure is INcreasing with n]
  #
  # [ -1 if threshold_len, +1 if thresold_level ]

  factor <- +1

  #**************************************************************************
  #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  quantity.to.measure <- ifelse(factor == +1,0,2*threshold)

  step <- 2^step

  n1 <- 0

  found.upper.bound <- FALSE
  direction <- +1

  while(step>=1)
  {
    while(sign(factor*(threshold-quantity.to.measure)) == direction && step >= 1)
    {
      step[found.upper.bound] <- max(1,step/2)

      n1 <- n1+direction*step

      #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      #************************************************************************
      # Define n2 from n1

      n2 <- n1

      #*********************************

      # Here compute the worst coverage

      # worst x's are defined by (17) in
      # Joseph, du Berger, Belisle, for general c's and d's, here set to 1.

      x1 <- ceiling(n1/2)
      x2 <- ceiling(n2/2)

      posterior.mean <- (1+x2)/(2+n2)-(1+x1)/(2+n1)
      posterior.var <- (x1+1)*(n1-x1+1)/(n1+2)^2/(n1+3)+
                    (x2+1)*(n2-x2+1)/(n2+2)^2/(n2+3)

      # We make the approximation of the difference between two betas by
      # another beta, with parameters a and b given by

      a <- -1/2 * (posterior.mean + 1) * (posterior.mean^2 - 1 + posterior.var)/
         posterior.var
      b <- 1/2*(1-posterior.mean)*(-1*posterior.mean^2+1-posterior.var)/
         posterior.var

      # and approximate the coverage probability of the HPD region by the
      # coverage obtained by a symetric region around the posterior mean

      pi.min <- posterior.mean-len/2
      pi.max <- posterior.mean+len/2
      pi.max[pi.min< -1] <- -1+len
      pi.min[pi.min< -1] <- -1
      pi.min[pi.max>1] <- 1-len
      pi.max[pi.max>1] <- 1
      coverage.around.mean <- pbeta((pi.max+1)/2,a,b)-pbeta((pi.min+1)/2,a,b)

      quantity.to.measure <- mean(coverage.around.mean)

      #************************************************************************
      #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    }

    found.upper.bound <- TRUE

    direction <- -direction
    step[step==1] <- 0
  }

  direction[n1==0] <- 0

  n1[direction==+1] <- n1+1

  # Return

  c(n1,n1)
}

