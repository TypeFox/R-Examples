#' Run BCD on a subset of frequencies
#' 
#' \code{bcd_express} runs bcd_inexact on a subset of frequencies determined by Lomb-Scargle fits that is guaranteed to include the minimum.
#' 
#' @param tms list of matrices whose rows are the triple (t,m,sigma) for each band.
#' @param sol_ls list of lists where sol_ls is output from lomb_scargle
#' @param omega_seq vector of freqencies used for producing sol
#' @param gamma1 nonnegative regularization parameter for shrinking amplitudes
#' @param gamma2 nonnegative regularization parameter for shrinking phases
#' @param at prior vector
#' @param max_iter maximum number of outer iterations for bcd
#' @param verbose boolean flag; if TRUE diagnostic info is printed
#' @export
bcd_express <- function(tms,sol_ls,omega_seq,gamma1,gamma2,at,
                        max_iter=100,verbose=FALSE){
  nOmega <- length(omega_seq)
  B <- length(tms)
  
  ## for each omega, sum rss across bands
#  rssLS <- rep(0,length=nOmega)
#  for(sol_ls in sol){
#    rssLS <- rssLS + vapply(sol_ls,function(x){x$RSS},c(0))
#  }
  rssLS <- rep(0,length=nOmega)
  rssLS <- colSums(sol_ls[,,4,drop=FALSE])
  
  ## order least square solutions by
  ## sum of rsses across bands
  ord <- order(rssLS)
  sol_ls <-sol_ls[,ord,,drop=FALSE]
#  for(ii in 1:length(sol)){
#    sol[[ii]] <- sol[[ii]][ord]
#  }
  rssLS <- rssLS[ord]
  omega_seq <- omega_seq[ord]
  
  ## run bcd on subset of omega

  sol.names <- list(band=names(tms),
                  NULL,
                  param=c("beta0","amp","rho"))
  sol_bcd <- array(NA,dim=c(B,nOmega,3),dimnames=sol.names)

#  sol_bcd <- vector(mode="list",length=nOmega)
  rss_bcd <- double(nOmega)
  ii <- 1
  max_freqs <- nOmega # maximum number of frequencies to try
  while (ii <= max_freqs){
    if(verbose){
      print(paste("bcd express iteration / max: ",ii,"/",max_freqs,sep=""))
      print(paste("frequency: ",omega_seq[ii]))
    }
    ## get inital params for bcd
#    beta0 <- vapply(sol,function(x){x[[ii]]$beta0},c(0))
#    A <- vapply(sol,function(x){x[[ii]]$A},c(0))
#    rho <- vapply(sol,function(x){x[[ii]]$rho},c(0))   
    beta0 <- sol_ls[,ii,1]
    A <- sol_ls[,ii,2]
    rho <- sol_ls[,ii,3]
    ## compute bcd solution and pnll
    tmp <- bcd_inexact(tms,
    #sol_bcd[[ii]] <- bcd_inexact(tms,
                                 beta0,
                                 A,
                                 at,
                                 rho,
                                 omega_seq[ii],
                                 gamma1,
                                 gamma2,
                                 max_iter=max_iter)
    sol_bcd[,ii,1] <- tmp$beta0
    sol_bcd[,ii,2] <- tmp$A
    sol_bcd[,ii,3] <- tmp$rho
    rss_bcd[ii] <- tmp$loss[tmp$iter]
#    rss_bcd[ii] <- sol_bcd[[ii]]$loss[sol_bcd[[ii]]$iter]
    
    ## only evaluate bcd at omega where
    ## unpenalized likelihood (rssLS)
    ## is less than rss_bcd[ii]
    max_freqs <- sum(rss_bcd[ii] > rssLS[1:max_freqs])
    ii <- ii + 1
  }
  ## return omegas where bcd was evaluated
  omega_bcd <- omega_seq[1:max_freqs]
#  sol_bcd <- sol_bcd[1:max_freqs] ## what happens if max_freqs==1, could this return error
  sol_bcd <- sol_bcd[,1:max_freqs,,drop=FALSE]
  rss_bcd <- rss_bcd[1:max_freqs]

  omega_ord <- order(omega_bcd,decreasing = FALSE)
  omega_bcd <- omega_bcd[omega_ord]
  sol_bcd <- sol_bcd[,omega_ord,,drop=FALSE]
  rss_bcd <- rss_bcd[omega_ord]
  return(list(omega_bcd,sol_bcd,rss_bcd))
}
