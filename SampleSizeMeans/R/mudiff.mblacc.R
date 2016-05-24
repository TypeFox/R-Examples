`mudiff.mblacc` <-
function(len,alpha1,beta1,alpha2,beta2,level=0.95,m=10000,mcs=3)
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

  step <- ceiling(log(mudiff.freq(len,alpha1/beta1,alpha2/beta2,level)[1])/log(2))

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
      if(n1 <= 2) {
        found.lower.bound <- TRUE
        n1 <- 2
      }

      cons.steps.same.dir <- cons.steps.same.dir+1

      history.ns <- c(n1,history.ns)
      history.steps <- c(step*direction,history.steps)

      #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
      #************************************************************************
      # Define n2 from n1

      n2 <- n1

      #*********************************

      # Let total.var=n1*s21/n1/(n1-1)+n2*s22/n2/(n2-1)

      total.var <- rgg(m,alpha1,2*beta1/n1/(n1-1),(n1-1)/2)+
                rgg(m,alpha2,2*beta2/n2/(n2-1),(n2-1)/2)

      # Average coverage is

      quantity.to.measure <- 2*mean(pnorm(len/2/sqrt(total.var)))-1

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

  # Return

  c(n1,n1)
}

