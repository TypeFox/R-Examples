`propdiff.mblmodwoc` <-
function(len,c1,d1,c2,d2,level=0.95,worst.level=.95)
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
  quant <- qchisq(worst.level,2)

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

      #*****************************************
      # Compute the worst (worst-level) coverage

      # worst x's are defined by (17) in
      # Joseph, du Berger, Belisle

      x1 <- ifelse(n1 >= abs(c1-d1),ceiling((n1-c1+d1)/2),n1)
      x2 <- ifelse(n2 >= abs(c2-d2),ceiling((n2-c2+d2)/2),n2)

      # Preposterior moments of x1 and x2

      m1 <- n1*c1/(c1+d1)
      v1 <- n1*c1*d1*(n1+c1+d1)/(c1+d1)^2/(c1+d1+1)
      m2 <- n2*c2/(c2+d2)
      v2 <- n2*c2*d2*(n2+c2+d2)/(c2+d2)^2/(c2+d2+1)

      # (x1,x2) is the worst outcome possible;
      # see if it is in the HPD region of level
      # worst-level (approximately ellipsoidal)

      distance <- (x1-m1)^2/v1+(x2-m2)^2/v2
      if (distance > quant)
      {
        # If the worst possible outcome is not in HPD region of indicated
        # level, then the worst outcome of that level is on the border
        # on the HPD region (on the border of the ellipse)

        # Normal approximation for distributions of X1 and X2

        x11 <- seq(floor(m1-sqrt(v1*quant)),ceiling(m1+sqrt(v1*quant)))
        # Eliminate impossible values and doubles
        x11[x11<0] <- 0
        x11[x11>n1] <- n1
        x11[x11[-length(x11)]!=x11[-1]]
        # Prepare the lower part of the ellipse
        x12 <- x11[length(x11):1]

        # Higher part of the ellipse
        s1 <- v2*(quant-(x11-m1)^2/v1)
        s1[s1<0] <- 0
        x21 <- ceiling(m2+sqrt(s1))
        # Lower part of the ellipse
        s2 <- v2*(quant-(x12-m1)^2/v1)
        s2[s2<0] <- 0
        x22 <- floor(m2-sqrt(s2))

        x1 <- c(x11,x12,x11[1])
        x2 <- c(x21,x22,x21[1])
        x2[x2<0] <- 0
        x2[x2>n2] <- n2
      }

      # Compute the posterior mean & variance
      # (for points on the border of the ellipse if worst possible outcome
      #  is not in the HPD region)
      # setting c's and d's to 1's.

      posterior.mean <- (1+x2)/(2+n2)-(1+x1)/(2+n1)
      posterior.var <- (x1+1)*(n1-x1+1)/(n1+2)^2/(n1+3)+
                    (x2+1)*(n2-x2+1)/(n2+2)^2/(n2+3)

      # Compute the maximal posterior mean & variance
      # (necessary if considering points on the border of
      #  the ellipse)

      posterior.mean <- posterior.mean[posterior.var==max(posterior.var)][1]
      posterior.var <- max(posterior.var)[1]

      pi.min <- posterior.mean-len/2
      pi.max <- posterior.mean+len/2
      pi.max[pi.min< -1] <- -1+len
      pi.min[pi.min< -1] <- -1
      pi.min[pi.max>1] <- 1-len
      pi.max[pi.max>1] <- 1

      # We make the approximation of the difference between two betas by
      # another beta, with parameters a and b given by

      a <- -1/2 * (posterior.mean + 1) * (posterior.mean^2 - 1 + posterior.var)/
         posterior.var
      b <- 1/2*(1-posterior.mean)*(-1*posterior.mean^2+1-posterior.var)/
         posterior.var

      quantity.to.measure <- pbeta((1+pi.max)/2,a,b)-pbeta((1+pi.min)/2,a,b)

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

