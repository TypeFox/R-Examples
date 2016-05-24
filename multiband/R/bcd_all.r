#' Run BCD on all frequencies
#' 
#' \code{bcd_all} runs bcd_inexact on all input frequencies.
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
bcd_all <- function(tms,sol_ls,omega_seq,gamma1,gamma2,at,
                        max_iter=100,verbose=FALSE){
  nOmega <- length(omega_seq)
  B <- length(tms)  
    
  ## run bcd on all omega
  sol.names <- list(band=names(tms),
                    NULL,
                    param=c("beta0","amp","rho"))
  sol_bcd <- array(NA,dim=c(B,nOmega,3),dimnames=sol.names)
#  sol_bcd <- vector(mode="list",length=nOmega)
  rss_bcd <- double(nOmega)
  for (ii in 1:nOmega) {
    if(verbose){
      print(paste("bcd all iteration / max: ",ii,"/",nOmega,sep=""))
      print(paste("frequency: ",omega_seq[ii]))
    }
    ## get inital params for bcd
#    beta0 <- vapply(sol_ls,function(x){x[[ii]]$beta0},c(0))
#    A <- vapply(sol_ls,function(x){x[[ii]]$A},c(0))
#    rho <- vapply(sol_ls,function(x){x[[ii]]$rho},c(0))
    beta0 <- sol_ls[,ii,1]
    A <- sol_ls[,ii,2]
    rho <- sol_ls[,ii,3]
    ## compute bcd solution and pnll
    tmp <- bcd_inexact(tms,
#    sol_bcd[[ii]] <- bcd_inexact(tms,
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
  }
  ## return omegas where bcd was evaluated
  return(list(omega_seq,sol_bcd,rss_bcd))
}
