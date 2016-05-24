## MCMC updates via Hamiltonian (Hybrid) Monte Carlo (HMC).
## f  = lik
## q  = x.init
## fq = y.init
## w, lower, upper, control: as usual.
##
## f(x) is supposed to compute y with an attribute 'gr' that is the
## gradient.
sampler.hmc.bounded <- function(f, q, fq, w, lower, upper,
                                control) {
  ## The issue here is that I don't get a chance to do any error
  ## checking on these parameters without directly altering mcmc(),
  ## which seems silly.  I might go with a similar sort of thing to
  ## find.mle() where rather than passing in samplers directly, we can
  ## just pass in strings, and work off those?  Might be easier,
  ## especially to do the error checking.
  epsilon <- control$epsilon
  L <- control$L
  m <- control$m
  
  if ( m != 1 )
    stop("Varying m not allowed yet.")
  if ( any(upper < Inf) )
    stop("Probably should not try upper yet.")

  q.in <- q
  fq.in <- fq
  gq <- attr(fq, "gr")

  p.in <- p <- rnorm(length(q), 0, 1)

  ## Make a half step for momentum at the beginning
  p <- p + epsilon * gq / 2

  ## See p. 36-37 and figure 8 of Neal 2011 to see how to handle
  ## constraints.  These are added here, I think.

  ## Alternate full steps for position and momentum
  for (i in seq_len(L)) {
    q <- q + epsilon * p # position

    ## Bounds check:
    nok <- q < lower
    if ( any(nok) ) {
      q[nok] <- 2 * lower[nok] - q[nok]
      p[nok] <- -p[nok]
    }
    
    if ( i != L )
      p <- p + epsilon * attr(f(q), "gr") # momentum
  }

  ## Make a half step for momentum at the end.
  fq.out <- f(q)
  p <- p + epsilon * attr(fq.out, "gr") / 2

  k.in <- sum(p.in^2) / 2
  k.out <- sum(p^2) / 2

  alpha <- exp(-fq.in - -fq.out + k.in - k.out)  

  if ( runif(1) < alpha )
    list(q, fq.out) # accept
  else
    list(q.in, fq.in) # reject
}
