#' @title Rippl analysis
#' @description Computes the Rippl no-failure storage for given time series of inflows and releases using the sequent peak algorithm. 
#' @param Q             vector or time series object. Net inflow totals to the reservoir.
#' @param target        a target release constant in same volumteric units as Q. Can be omitted if R is given.
#' @param R             a time series or vector of target releases (volumetric). Must be the same length as Q.
#' @param double_cycle  logical. If TRUE the Q and R time series will be replicated and placed end-to-end to double the simulation. Recommended if the critical period occurs at the end of the sequence.
#' @param plot          logical. If TRUE (the default) the storage behavior diagram is plotted.
#' @return Returns the no-fail storage capacity and corresponding storage behaviour time series.
#' @examples # define a release vector for a constant release equal to 90 % of the mean inflow
#' no_fail_storage <- Rippl(resX$Q_Mm3, target = 0.9 * mean(resX$Q_Mm3))$No_fail_storage
#' @import stats
#' @export
Rippl <- function(Q, target, R, double_cycle = FALSE, plot = TRUE) {
  
    if (is.ts(Q) == FALSE && is.vector(Q) == FALSE)
        stop("Q must be time series or vector object")

    if (!missing(target) && !missing(R))
        stop("Supply target or R, not both!")
    
    if (!missing(target))
        R <- rep(target, length(Q))
    
    if (length(Q) != length(R))
      stop("Q and R must be of equal length")
    
    if (is.ts(R) == FALSE && is.vector(R) == FALSE)
      stop("R must be time series or vector object")
    
    if (double_cycle) {
      Q <- ts(c(Q, Q), start = start(Q), frequency = frequency(Q))
      R <- c(R, R)
    }
    
    
    K <- vector("numeric", length = length(Q) + 1)
    Spill <- vector("numeric", length = length(Q))
    for (t in 1:length(Q)) {
        if (R[t] - Q[t] + K[t] > 0) {
            K[t + 1] <- R[t] - Q[t] + K[t]
        } else {
            K[t + 1] <- 0
            Spill[t] <- abs(R[t] - Q[t] + K[t])
        }
    }
    K <- ts(K[1:(length(K) - 1)], start = start(Q), frequency = frequency(Q))
    Spill <- ts(Spill, start = start(Q), frequency = frequency(Q))
    if (plot) {
        plot(max(K) - K, ylab = "Storage", main = paste0("No Fail Storage = ", round(max(K),3)))
        plot(Q, ylab = "Inflow totals")
        plot(Spill, ylab = "Uncontrolled spill")
    }
    results <- list(max(K),(max(K) - K), Spill)
    names(results) <- c("No_fail_storage","Storage_behavior", "Uncontrolled_spill")
    return(results)
}
