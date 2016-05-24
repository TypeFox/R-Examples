#' Single band Lomb-Scargle
#' 
#' \code{single_band_lomb_scargle} runs Lomb-Scargle on all bands in tms at all frequencies in omega_seq
#' 
#' @param tms List of lightcurves
#' @param omega_seq Sequence of frequencies
#' @export
single_band_lomb_scargle <- function(tms, omega_seq) {
  nOmega <- length(omega_seq)
  B <- length(tms)
  sol.names <- list(band=names(tms),
                    NULL,
                    param=c("beta0","amp","rho","rss"))
  sol <- array(NA,dim=c(B,nOmega,4),dimnames=sol.names)
#  sol <- vector(mode="list",length=length(tms))
#  names(sol) <- names(tms)
  last_success <- c(0,1,0)
  for(ii in 1:B) {
#    sol[[ii]] <- vector(mode="list",length=nOmega)
    for(jj in 1:nOmega) {
#      sol[[ii]][[jj]] <- NA
      tmp <- try(lomb_scargle(tms[[ii]][,1],
                              tms[[ii]][,2],
                              tms[[ii]][,3],
                              omega_seq[jj]),silent=TRUE)
      if (!inherits(tmp, "try-error")) {
#        sol[[ii]][[jj]] <- tmp
        sol[ii,jj,] <- unlist(tmp)
        ## put on same scale as bcd        
        sol[ii,jj,4] <- sol[ii,jj,4] / 2
#        ## put on same scale as bcd
#        sol[[ii]][[jj]]$RSS <- sol[[ii]][[jj]]$RSS / 2
        last_success <- sol[ii,jj,1:3]
      } else {
        sol[ii,jj,1:3] <- last_success        
        sol[ii,jj,4] <- 0
      }
#      } else {
#        print(paste(ii,jj))
#      }
    }
  }
  return(sol)
}