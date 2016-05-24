#' @title Storage-Reliability-Yield (SRY) relationships: Yield computation
#' @description Returns the yield for given inflow time series, reservoir capacity, and required time-based reliability. Assumes standard operating policy. Yield is computed iteratively using the bi-section method.
#' @param Q               vector or time series object. Net inflow totals to the reservoir. 
#' @param capacity        numerical. The reservoir storage capacity. Must be in the same volumetric units as Q.
#' @param reliability     numerical. The required time-based reliability.
#' @param demand_profile  a vector of factors with length = frequency(Q). Represents within-year demand profile. Defaults to constant release if left blank.
#' @param plot            logical. If TRUE (the default) the storage behavior diagram and release time series are plotted.
#' @param S_initial       numeric. The initial storage as a ratio of capacity (0 <= S_initial <= 1). The default value is 1.
#' @param max_iterations  Maximum number of iterations for yield computation.
#' @param double_cycle    logical. If TRUE the input series will be replicated and placed end-to-end to double the simulation. (Recommended if the critical period occurs at the end of the recorded inflow time series)
#' @return Returns yield of a reservoir with specified storage capacity and time-based reliability. 
#' @examples # Compute yield for 0.95 reliability
#' layout(1:3)
#' yield_ResX <- yield(resX$Q_Mm3, capacity = 500, reliability = 0.95)
#' # Compute yield for quarterly time series with seasonal demand profile
#' 
#' quart_ts <- aggregate(resX$Q_Mm3, nfrequency = 4)
#' yld <- yield(quart_ts,
#' capacity = 500, reliability = 0.9, demand_profile = c(0.8, 1.2, 1.2, 0.8))
#' @import stats
#' @export
yield <- function(Q, capacity, reliability, demand_profile,
                  plot = TRUE, S_initial = 1, max_iterations = 50, double_cycle = FALSE) {
    
    if (missing(demand_profile))
      demand_profile <- rep(1, frequency(Q))
    if (length(demand_profile) != frequency(Q)) 
        stop("demand_profile must have length equal to the time series frequency")
    if (reliability < 0 || reliability > 1)
        stop("Reliability must be between 0 and 1")
  
    if (double_cycle) {
        Q <- ts(c(Q, Q), start = start(Q), frequency = frequency(Q))
    }
  
    S <- vector("numeric", length = length(Q) + 1); S[1] <- capacity * S_initial
    R <- vector("numeric", length = length(Q))
    min_yield <- 0
    max_yield <- 10 * mean(Q)
    i <- 0
    repeat {
        Spill <- vector("numeric", length(Q))
        i <- i + 1
        mid_yield <- (min_yield + max_yield) / 2
        R_target <- rep(mid_yield, length(Q)) * rep(demand_profile, length(Q) / frequency(Q))
        for (t in 1:length(Q)) {
            x <- Q[t] - R_target[t] + S[t]
            if (x < 0) {
                S[t + 1] <- 0
                R[t] <- Q[t] + S[t]
            } else {
                if (x > capacity) {
                  S[t + 1] <- capacity
                  R[t] <- R_target[t]
                  Spill[t] <- Q[t] - R_target[t] + S[t] - capacity
                } else {
                  S[t + 1] <- x
                  R[t] <- R_target[t]
                }
            }
        }
        reliability_0 <- 1 - (sum( (R / R_target) < 1) / length(Q))
        if (reliability_0 >= reliability) {
            min_yield <- mid_yield
        }
        if (reliability_0 < reliability) {
            max_yield <- mid_yield
        }
        if (max_yield - min_yield < 0.01) {
            break
        }
        if (i >= max_iterations) {
            break
        }
    }
    
    S <- ts(S[1:(length(S) - 1)], start = start(Q), frequency = frequency(Q))
    R <- ts(R, start = start(Q), frequency = frequency(Q))
    Spill <- ts(Spill, start = start(Q), frequency = frequency(Q))
    
    
    if (plot) {
        plot(S, ylab = "Storage", ylim = c(0,max(S)), main = (paste0("Yield = ",round(mid_yield, 2), " at ", reliability * 100, " % reliability")))
        plot(R, ylim = c(0, max(R)), ylab = "Water supplied")
        plot(Spill, ylab = "Spill")
    }
    message(paste0("Converged after ", i,
                 " iterations and time-based reliability equal to ", reliability_0))
    
    results <- list(mid_yield, S, R, Spill)
    names(results) <- c("Yield", "Storage", "Water_supplied", "Spill")
    return(results)
}
