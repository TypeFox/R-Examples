################################################################################
# 
# rlaptrans taken from Martin Ridout: 
# Ridout, M.S. (2009) Generating random numbers from a distribution specified 
# by its Laplace transform. Statistics and Computing, 19, 439-450.
# http://www.kent.ac.uk/smsas/personal/msr/rlaptrans.html
# 
#======================================================================================
rlaptrans <- function(n, ltpdf, ..., tol=1e-7, x0=1, xinc=2, m=11, L=1, A=19, nburn=38)
#======================================================================================
{
  
  #---------------------------------------------------------
  # Function for generating a random sample of size n from a 
  # distribution, given the Laplace transform of its p.d.f.
  #---------------------------------------------------------
  
  maxiter = 500
  
  # -----------------------------------------------------
  # Derived quantities that need only be calculated once,
  # including the binomial coefficients
  # -----------------------------------------------------
  nterms = nburn + m*L
  seqbtL = seq(nburn,nterms,L)
  y = pi * (1i) * seq(1:nterms) / L
  expy = exp(y)
  A2L = 0.5 * A / L
  expxt = exp(A2L) / L
  coef = choose(m,c(0:m)) / 2^m
  
  
  # --------------------------------------------------
  # Generate sorted uniform random numbers. xrand will
  # store the corresponding x values
  # --------------------------------------------------
  u = sort(runif(n), method="qu")
  xrand = u
  
  #------------------------------------------------------------
  # Begin by finding an x-value that can act as an upper bound
  # throughout. This will be stored in upplim. Its value is
  # based on the maximum value in u. We also use the first
  # value calculated (along with its pdf and cdf) as a starting
  # value for finding the solution to F(x) = u_min. (This is
  # used only once, so doesn't need to be a good starting value
  #------------------------------------------------------------
  t = x0/xinc
  cdf = 0   
  kount0 = 0
  set1st = FALSE
  while (kount0 < maxiter & cdf < u[n]) {
    t = xinc * t
    kount0 = kount0 + 1
    x = A2L / t
    z = x + y/t
    ltx = ltpdf(x, ...)
    ltzexpy = ltpdf(z, ...) * expy
    par.sum = 0.5*Re(ltx) + cumsum( Re(ltzexpy) )
    par.sum2 = 0.5*Re(ltx/x) + cumsum( Re(ltzexpy/z) )
    pdf = expxt * sum(coef * par.sum[seqbtL]) / t
    cdf = expxt * sum(coef * par.sum2[seqbtL]) / t
    if (!set1st & cdf > u[1]) {
      cdf1 = cdf
      pdf1 = pdf
      t1 = t
      set1st = TRUE
    }
  }
  if (kount0 >= maxiter) {
    stop('Cannot locate upper quantile')
  }
  upplim = t
  
  #--------------------------------
  # Now use modified Newton-Raphson
  #--------------------------------
  
  lower = 0
  t = t1
  cdf = cdf1
  pdf = pdf1
  kount = numeric(n)
  
  maxiter = 1000
  
  for (j in 1:n) {
    
    #-------------------------------
    # Initial bracketing of solution
    #-------------------------------
    upper = upplim
    
    kount[j] = 0
    while (kount[j] < maxiter & abs(u[j]-cdf) > tol) {
      kount[j] = kount[j] + 1
      
      #-----------------------------------------------
      # Update t. Try Newton-Raphson approach. If this 
      # goes outside the bounds, use midpoint instead
      #-----------------------------------------------
      t = t - (cdf-u[j])/pdf 
      if (t < lower | t > upper) {
        t = 0.5 * (lower + upper)
      }
      
      #----------------------------------------------------
      # Calculate the cdf and pdf at the updated value of t
      #----------------------------------------------------
      x = A2L / t
      z = x + y/t
      ltx = ltpdf(x, ...)
      ltzexpy = ltpdf(z, ...) * expy
      par.sum = 0.5*Re(ltx) + cumsum( Re(ltzexpy) )
      par.sum2 = 0.5*Re(ltx/x) + cumsum( Re(ltzexpy/z) )
      pdf = expxt * sum(coef * par.sum[seqbtL]) / t
      cdf = expxt * sum(coef * par.sum2[seqbtL]) / t
      
      #------------------
      # Update the bounds 
      #------------------
      if (cdf <= u[j]) {
        lower = t}
      else {
        upper = t}
    }
    if (kount[j] >= maxiter) {
      warning('Desired accuracy not achieved for F(x)=u')
    }
    xrand[j] = t
    lower = t
  }
  
  if (n > 1) {
    rsample <- sample(xrand) }
  else {
    rsample <- xrand} 
  rsample
}
