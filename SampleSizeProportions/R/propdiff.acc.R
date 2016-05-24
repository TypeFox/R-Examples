`propdiff.acc` <-
function(len,c1,d1,c2,d2,level=0.95,equal=TRUE,m=10000,mcs=3)
{
  min.for.possible.return <- 2^ceiling(1.5*mcs)

  # If we always allow a return, there is a risk of making bad steps
  # when we are close to the answer.
  # Thus, we should not allow any return once some arbitrary 'step' (which is
  # 'min.for.possible.return') is reached.

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

  history.ns <- 0
  history.steps <- 0

  n1 <- 0

  max.cons.steps.same.dir <- mcs
  found.upper.bound <- FALSE
  possible.to.move.back <- TRUE
  cons.steps.same.dir <- 0
  direction <- +1

  while(step>=1)
  {
    while(sign(factor*(threshold-quantity.to.measure)) == direction && step >= 1)
    {
      step[found.upper.bound] <- max(1,step/2)
      possible.to.move.back[step < min.for.possible.return &&
                            found.upper.bound] <- FALSE



      n1 <- n1+direction*step

      cons.steps.same.dir <- cons.steps.same.dir+1

      history.ns <- c(n1,history.ns)
      history.steps <- c(step*direction,history.steps)


      #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      #************************************************************************
      # Define n2 from n1

      n2 <- ifelse(equal,n1,ceiling(sqrt(c2*d2/c1/d1*(c1+d1)*(c1+d1+1)/(c2+d2)/
                               (c2+d2+1))*(c1+d1+n1)-c2-d2))

      n2[n2<0] <- 0

      #**********************************
      # Here compute the average coverage

      pi1 <- rbeta(m,c1,d1)
      pi2 <- rbeta(m,c2,d2)
      x1 <- rbinom(m,n1,pi1)
      x2 <- rbinom(m,n2,pi2)

      # Precaution: if n1 or n2 is 0, then the correponding x given
      # by rbinom is NA. We correct the situation by setting it to 0.

      x1[is.na(x1)] <- 0
      x2[is.na(x2)] <- 0

      # Posterior variance of the difference between the two proportions

      posterior.mean <- (c2+x2)/(c2+d2+n2)-(c1+x1)/(c1+d1+n1)
      posterior.var <- (x1+c1)*(n1-x1+d1)/(n1+c1+d1)^2/(n1+c1+d1+1)+
        (x2+c2)*(n2-x2+d2)/(n2+c2+d2)^2/(n2+c2+d2+1)

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

      if(found.upper.bound &&
         cons.steps.same.dir == max.cons.steps.same.dir+1 &&
         possible.to.move.back)
      {

        if(sign(factor*(threshold-quantity.to.measure)) == direction)
        {
          # There was (most likely) a mistake, look for n's in the other direction
          at.n <- seq(along=history.ns)[history.ns == n1 ]
          hs.an <- history.steps[at.n]

          step <- abs(hs.an[sign(hs.an) != direction][1])

          cons.steps.same.dir <- 0

          # and if there has never been a step coming from the other direction...
          step[is.na(step)] <- max(abs(hs.an))
        }
        else
        {
          # There was (most likely) no mistake; keep looking around the same n's

          direction <- -direction
          cons.steps.same.dir <- 0
        }
      }

      if(found.upper.bound &&
         cons.steps.same.dir==max.cons.steps.same.dir &&
         sign(factor*(threshold-quantity.to.measure))==direction &&
         possible.to.move.back)
      {
        step <- 2*step
      }

    }

    found.upper.bound <- TRUE

    direction <- -direction
    cons.steps.same.dir <- 0
    step[step==1] <- 0


  }

  direction[n1==0] <- 0

  n1[direction==+1] <- n1+1

  #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  #****************************************************************************
  # Once again, define n2

  n2 <- ifelse(equal,n1,ceiling(sqrt(c2*d2/c1/d1*(c1+d1)*(c1+d1+1)/(c2+d2)/
                           (c2+d2+1))*(c1+d1+n1)-c2-d2))

  n2[n2<0] <- 0

  #****************************************************************************
  #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  # Return

  c(n1,n2)
}

