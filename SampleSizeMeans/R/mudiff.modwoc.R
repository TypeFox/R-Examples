`mudiff.modwoc` <-
function(len,alpha1,beta1,alpha2,beta2,n01,n02,level=0.95,worst.level=0.95,equal=TRUE,m=50000,mcs=3)
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

  step <- ceiling(log(mudiff.freq(len,alpha1/beta1,alpha2/beta2,level,equal)[1])/log(2))

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
  # history.cons.steps_0

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
      if(n1 <= 0) {
        found.lower.bound <- TRUE
        n1 <- 0
      }

      cons.steps.same.dir <- cons.steps.same.dir+1

      history.ns <- c(n1,history.ns)
      history.steps <- c(step*direction,history.steps)

      #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      #************************************************************************
      # Define n2 from n1

      n2 <- mudiff.samplesize1(alpha1,beta1,alpha2,beta2,n01,n02,equal,n1)

      #*********************************

      # The worst.level% 'worst' case is when the posterior variance takes
      # its worst.level-quantile value, that is, when bn takes
      # its worst.level-quantile
      # Estimate it (worstbn)

      x1 <- rt(m,2*alpha1)/sqrt(n01*n1/(n01+n1)*alpha1/beta1)
      x2 <- rt(m,2*alpha2)/sqrt(n02*n2/(n02+n2)*alpha2/beta2)

      b1 <- n1*n01/(n1+n01)*x1^2+2*beta1
      ns21 <- rvgg(alpha1+1/2,b1,(n1-1)/2)
      b2 <- n2*n02/(n2+n02)*x2^2+2*beta2
      ns22 <- rvgg(alpha2+1/2,b2,(n2-1)/2)

      bn1 <- beta1+ns21/2+1/2/(n01+n1)*n1*n01*x1^2
      bn2 <- beta2+ns22/2+1/2/(n02+n2)*n2*n02*x2^2

      cn1 <- (n1+n01)*(alpha1+n1/2)/bn1
      cn2 <- (n2+n02)*(alpha2+n2/2)/bn2

      dn <- (2*alpha1+n1)/(2*alpha1+n1-2)/cn1+(2*alpha2+n2)/(2*alpha2+n2-2)/cn2

      worstdn <- sort(dn)[worst.level*m]

      quantity.to.measure <- 2*pnorm(len/2/sqrt(worstdn))-1

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

          # and if there has never been a step coming from the other direction

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

  n2 <- mudiff.samplesize1(alpha1,beta1,alpha2,beta2,n01,n02,equal,n1)

  #****************************************************************************
  #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


  # Return

  c(n1,n2)

}

