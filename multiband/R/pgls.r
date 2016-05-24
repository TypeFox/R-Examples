#' Penalized Generalized Lomb-Scargle
#' 
#' \code{pgls} estimates periods for a collection of lightcurves sampled over multiple bands.
#' It borrows strength across multiple bands via shrinkage penalties on amplitudes an phases.
#' 
#' @param lclist list of lightcurve data frames
#' @param period_min minimum period
#' @param period_max maximum period
#' @param periods grid of periods
#' @param gamma1 vector of Amplitude regularization parameter
#' @param gamma2 vector of Phase regularization parameter
#' @param at amplitude prior parameter
#' @param LS_flag boolean whether to run Lomb-Scargle algorithm
#' @param sol_ls Lomb-Scargle solution, used if LS_flag=FALSE
#' @param BCD_flag boolean whether to run bcd algorithm
#' @param fast_BCD_flag boolean whether to run BCD on relevant subset of periods
#' @param max_iter maximum number of outer iterations - passed to bcd_inexact
#' @param tol tolerance on relative change in loss - passed to bcd_inexact
#' @param mm_iter number of MM iterations for rho update - passed to bcd_inexact
#' @param verbose boolean whether to print progress
#' @export
#' @examples
#' period_min <- 3.5
#' period_max <- 3.65
#' out <- pgls(cepii,period_min=period_min,period_max=period_max)
#' out$best_fitBCD
pgls <- function(lclist,period_min=NULL,period_max=NULL,periods=NULL,gamma1=0,gamma2=20,at=rep(1,length(lclist)),LS_flag=TRUE,sol_ls=NULL,BCD_flag=TRUE,fast_BCD_flag=TRUE,max_iter=1e2,tol=1e-4,mm_iter=5,verbose=FALSE) {

  ## Check consistency of inputs
  ## 1. atilde has length(lclist)
  ## 2. atilde > 0
  ## 3. gamma1, gamma2 nonnegative
  ## 4. period_min, period_max > 0
  ## 5. period sequence - ascending

    ## Normalize atilde
    at <- at/sqrt(sum(at**2))
  if (is.null(periods) & (is.null(period_min) | is.null(period_max))){
      stop("either period range (period_min and period_max) or a vector of periods (periods) must be non NULL")
  }

  if (is.null(periods)) {
    # get_freqs and output 'omega_seq'
    freq_del <- get_freq_del(lclist)
    omega_seq <- get_freqs(period_min,period_max,freq_del)
  } else {
    omega_seq <- sort(2*pi/periods)
    # convert periods to frequencies 'omega_seq'
  }
  
  nOmega <- length(omega_seq)
  B <- length(lclist)
  
  # Step 1: Run standard Lomb-Scargle
  rssLS <- NULL
  best_fitLS <- NULL
  if (LS_flag) {
    if(verbose) print('Running Lomb-Scargle')
    sol_ls <- single_band_lomb_scargle(lclist,omega_seq)  
}  else {
    if (is.null(sol_ls)){
        stop("when LS_flag is FALSE, sol_ls must be non-null")
    }
    ## sol_ls from user will be in period ascending, reverse this
    ix_rev_ls <- seq(nOmega,1)
    sol_ls=sol_ls[,ix_rev_ls,,drop=FALSE]
}
  ## Compute RSS for Lomb-Scargle
    rssLS <- rep(0,length=nOmega)
    rssLS <- colSums(sol_ls[,,4,drop=FALSE])
  ## Extract best Lomb-Scargle
    if(verbose) print('Extracting Best Lomb-Scargle')
    param.names <- c("beta0","A","rho")

    ix_min <- which.min(rssLS)
    best_fitLS <- matrix(0,nrow=B,ncol=3)
    colnames(best_fitLS) <- param.names
    rownames(best_fitLS) <- names(lclist)
    best_fitLS <- sol_ls[,ix_min,1:3]


  
  
  # Step 2: Run BCD
  rss_bcd <- NULL
  sol_bcd <- NULL
  best_fitBCD <- NULL
  if (BCD_flag) {
    if (fast_BCD_flag) {
      if(verbose) print('Running fast BCD')
      sol_bcd <- bcd_express(lclist,sol_ls,omega_seq,gamma1,gamma2,at,max_iter)
    } else {
      if(verbose) print('Running regular BCD')
      sol_bcd <- bcd_all(lclist,sol_ls,omega_seq,gamma1,gamma2,at,max_iter)
    }
    rss_bcd <- sol_bcd[[3]]
    
    ## Extract best BCD model
    ix_min <- which.min(rss_bcd)
    best_fitBCD <- matrix(0,nrow=B,ncol=3)
    colnames(best_fitBCD) <- param.names
    rownames(best_fitBCD) <- names(lclist)
#    best_fitBCD[,1] <- sol_bcd[[2]][[ix_min]]$beta0
#    best_fitBCD[,2] <- sol_bcd[[2]][[ix_min]]$A
#    best_fitBCD[,3] <- sol_bcd[[2]][[ix_min]]$rho
#     if (length(rss_bcd) > 1) {
#       best_fitBCD <- sol_bcd[[2]][,ix_min,]
#     } else {
#       best_fitBCD <- sol_bcd[[2]]
#     }
    best_fitBCD <- sol_bcd[[2]][,ix_min,]
    ## Note: sol_bcd[[2]] is a 3-way array
  }
  
  # Output:
  # 1. RSS for LS and BCD
  # 2. Parameters for each fit?
  # 3. Model for best RSS period
  ix_rev_ls <- seq(nOmega,1)
  ix_rev_bcd <- seq(length(sol_bcd[[1]]),1)
  return(list(rss_ls=rev(rssLS), rss_bcd=rev(rss_bcd), 
              best_fitLS=best_fitLS, best_fitBCD=best_fitBCD,
              period_seq_all= rev(2*pi/omega_seq), period_seq_bcd = rev(2*pi/sol_bcd[[1]]), 
              sol_ls=sol_ls[,ix_rev_ls,,drop=FALSE], sol_bcd=sol_bcd[[2]][,ix_rev_bcd,,drop=FALSE]))
}
