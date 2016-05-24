getARL = function(distr=NULL, K=NULL, H=NULL,
    Mean=NULL, std=NULL, prob=NULL, Var=NULL, mu=NULL, lambda=NULL, 
    samp.size=NULL, is.upward=NULL, winsrl=NULL, winsru=NULL) {
  if (is.null(distr)) {
    cat('1 = Normal location, \n 
         2 = Normal variance, \n
         3 = Poisson, \n 
         4 = Binomial, \n 
         5 = Negative binomial, \n
         6 = Inv Gaussian mean. \n')
    stop("Sepcify a distribution.")
  }
  distrs = c('Normal location','Normal variance',
    'Poisson', 'Binomial', 'Negative binomial',
    'Inv Gaussian mean')
  nonmen = length(distrs)
  if (!(distr %in% seq(nonmen))) stop ("'distr' is not a valid choice.")

  # parameter settings
  eps = 0.001
  denrat = 0.0
  isint = FALSE
  plus = 1.0
  arg2 = 0.0

  ## Winsorizing constants
  if (is.null(winsrl)) winsrl = -1.0E8 else winsrl = as.double(winsrl)
  if (is.null(winsru)) winsru = 1.0E8 else winsru = as.double(winsru)

  ## read parameters for different distributions
  if (distr == 1) {
    # Normal location
    if (is.null(Mean)) stop ("'Mean' is missing.")
    if (is.null(std)) stop ("Standard deviation 'std' is missing.")
    if (is.null(K)) stop ("Reference value 'K' is missing.")
    if (is.null(H)) stop ("Decision interval 'H' is missing.")
    amu = as.double(Mean)
    asig = as.double(std)
    ref = as.double(K)
    di = as.double(H)
    ref = (ref - amu) / asig
    di  = di / asig
    if (is.upward == TRUE || is.null(is.upward))
      cat('The cusum is assumed upward. \n') else
      cat('Always assume is.upward = TRUE for normal mean. \n')
    ndis = as.integer(1)
  }

  if (distr == 2) {
    ## Normal variance
    if (is.null(samp.size)) stop ("Rational group size 'samp.size' 
      is missing.") else nsamp = as.integer(samp.size)
    if (nsamp < 2) {
      cat("Charting squared deviations from true mean. \n")
      nsamp = as.integer(2)
    }
    integ = nsamp - as.integer(1)
    # change the common variable 'ndf' in varup & vardn subroutines
    .Fortran("cmpar2", integ)
    if (is.null(std)) stop ("Standard deviation 'std' is missing.")
    if (is.null(K)) stop ("Reference value 'K' is missing.")
    if (is.null(H)) stop ("Decision interval 'H' is missing.")
    sigma = as.double(std)
    ref = as.double(K)
    di = as.double(H)
    ref = ref / (sigma * sigma)
    di = di / (sigma * sigma)
    if (is.null(is.upward)) 
      stop ("Upward or Downward? 'is.upward' is missing.")
    ndis = as.integer(2)
    if (is.upward == FALSE) {
      plus = -1
      ndis = as.integer(3)
    }
  }

  if (distr == 3) {
    ## Poisson
    if (is.null(Mean)) stop ("Poisson mean 'Mean' is missing.")
    if (Mean <= 0) stop ("Poisson mean must be positive: Mean > 0.")
    arg2 = as.double(Mean)
    if (is.null(K)) stop ("Reference value 'K' is missing.")
    if (is.null(H)) stop ("Decision interval 'H' is missing.")
    ref = as.double(K)
    di = as.double(H)
    ndis = as.integer(4)
    if (is.null(is.upward)) 
      stop ("Upward or Downward? 'is.upward' is missing.")
    if (is.upward == FALSE) {
      plus = -1
      ndis = as.integer(5)
    }
    isint = TRUE
  }

  if (distr == 4) {
    ## Binomial
    if (is.null(samp.size)) stop ("Sample size 'samp.size' is missing.")
    if (is.null(prob)) stop ("Success probability 'prob' is missing.")
    integ = as.integer(samp.size)
    # change the common variable 'nbig' in binup & bindn subroutines
    .Fortran("cmpar2", integ)
    arg2 = as.double(prob)
    if (is.null(K)) stop ("Reference value 'K' is missing.")
    if (is.null(H)) stop ("Decision interval 'H' is missing.")
    ref = as.double(K)
    di = as.double(H)
    ndis = as.integer(6)
    if (is.null(is.upward)) 
      stop ("Upward or Downward? 'is.upward' is missing.")
    if (is.upward == FALSE) {
      plus = -1
      ndis = as.integer(7)
    } 
    isint = TRUE
  }

  if (distr == 5) {
    ## Negative binomial
    if (is.null(Mean)) stop ("'Mean' is missing.")
    if (is.null(Var)) stop ("Variance 'Var' is missing.")
    aver = as.double(Mean)
    varian = as.double(Var)
    arg2 = aver / (varian - aver)
    realno = aver * arg2
    if(max(realno, arg2) <= 0) 
      stop("Mean must be > 0 and variance > mean.")
    # change the common variable 'r' in nbinup & nbindn subroutines
    .Fortran("cmpar", realno)
    if (is.null(K)) stop ("Reference value 'K' is missing.")
    if (is.null(H)) stop ("Decision interval 'H' is missing.")
    ref = as.double(K)
    di = as.double(H)
    ndis = as.integer(8)
    if (is.null(is.upward)) 
      stop ("Upward or Downward? 'is.upward' is missing.")
    if (is.upward == FALSE) {
      plus = -1
      ndis = as.integer(9)
    }
    isint = TRUE
  }

  if (distr == 6) {
  ## Inverse Gaussian mean
    if (is.null(mu)) stop ("Inverse Gaussian mean 'mu' is missing.")
    if (is.null(lambda)) stop ("Shape parameter 'lambda' is missing.")
    if(min(mu, lambda) <= 0) 
      stop("All parameters must be strictly positive.")
    arg2 = as.double(mu)
    realno = as.double(lambda)
    # change the common variable 'alam' in gauiup & gauidn subroutines
    .Fortran("cmpar", realno)
    if (is.null(K)) stop ("Reference value 'K' is missing.")
    if (is.null(H)) stop ("Decision interval 'H' is missing.")
    ref = as.double(K)
    di = as.double(H)
    ndis = as.integer(10)
    if (is.null(is.upward)) 
      stop ("Upward or Downward? 'is.upward' is missing.")
    if (is.upward == FALSE){
      plus = -1
      ndis = as.integer(11)
    }
  }

  if (isint) {
    den = .Fortran('getden', ref=as.double(ref), di=as.double(di), 
      denrat=as.double(denrat), as.logical(TRUE), as.logical(FALSE))
    ref = den$ref
    di = den$di
    denrat = den$denrat
  }
  # Integer - find rational denominator
  cat(sprintf("( k = %9.4f,  h = %9.4f). \n", ref, di))
  winlo = winsrl
  winhi = winsru
  if (plus < 0) {
    winlo = -winsru
    winhi = -winsrl
  }
  cal = .Fortran('calcus', as.double(di), as.double(plus*ref), 
    as.integer(ndis), as.double(denrat), as.double(arg2), 
    as.double(winlo), as.double(winhi), regarl=as.double(0), 
    firarl=as.double(0), ssarl=as.double(0), eps=eps, 
    esterr=as.double(0), ifault=as.integer(0))
  regarl = cal$regarl
  firarl = cal$firarl
  ssarl = cal$ssarl
  esterr = cal$esterr
  ifault = cal$ifault
  if (ifault > 0) {
    cat(sprintf("Error code %5i was returned. \n", ifault))
  } 
  cat(sprintf('zero start, FIR, steady state ARLs 
    %10.2f %10.2f %10.2f \n', regarl, firarl, ssarl))
  if (ifault != 0) {
    good = -log10(abs(esterr))
    cat(sprintf("I think this answer has only %7.1f accurate digits. \n", 
      good))
  }
  ## Output the result
  list(ARL_Z = regarl, ARL_F = firarl, ARL_S = ssarl)
} 

