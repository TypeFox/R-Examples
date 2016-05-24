#' @title Storage-Reliability-Yield (SRY) relationships: Storage computation
#' @description Returns the required storage for given inflow time series, yield, and target time-based reliability. Assumes standard operating policy. Storage is computed iteratively using the bi-section method.
#' @param Q               vector or time series object. Net inflow totals to the reservoir. Recommended units: Mm^3 (Million cubic meters).
#' @param yield           the required yield. Must be same volumetric units as Q.
#' @param reliability     numerical. The required time-based reliability.
#' @param demand_profile  a vector of factors with length = frequency(Q). Represents within-year demand profile. Defaults to constant release if left blank.
#' @param plot            logical. If TRUE (the default) the storage behavior diagram and release time series are plotted.
#' @param S_initial       numeric. The initial storage as a ratio of capacity (0 <= S_initial <= 1). The default value is 1.
#' @param max_iterations  Maximum number of iterations for yield computation.
#' @param double_cycle    logical. If TRUE the input series will be replicated and placed end-to-end to double the simulation. (Recommended if the critical period occurs at the end of the recorded inflow time series)
#' @return Returns the required storage capacity necessary to supply specified yield with specified reliability.
#' @examples # Determine the required storage for 95 % reliability and yield equal to 80 % of the mean inflow.
#' layout(1:3)
#' storage(resX$Q_Mm3 * 20, yield = 0.9 * mean(resX$Q_Mm3), reliability = 0.95)
#' @import stats
#' @export
storage <- function(Q, yield, reliability, demand_profile,
                    plot = TRUE, S_initial = 1, max_iterations = 50, double_cycle = FALSE) {
    
    if (missing(demand_profile))
        demand_profile <- rep(1, frequency(Q))
    if (length(demand_profile) != frequency(Q))
        stop("demand_profile must have length equal to the time series frequency")
    if (reliability < 0 || reliability > 1)
        stop("Reliability must be between 0 and 1")
    if (length(yield) > 1)
        stop("yield must be a scalar, not a vector")
    if (double_cycle) {
        Q <- ts(c(Q, Q), start = start(Q), frequency = frequency(Q))
    }

    S <- vector("numeric", length = length(Q) + 1)
    R <- vector("numeric", length = length(Q))
    min_storage <- 0
    max_storage <- 100 * mean(Q) * frequency(Q)  #Provides upper bound of 100 years' storage
    i <- 0
    repeat {
        Spill <- vector("numeric", length(Q))
        i <- i + 1
        mid_storage <- (min_storage + max_storage) / 2
        S[1] <- mid_storage * S_initial
        R_target <- rep(yield, length(Q)) * rep(demand_profile, length(Q) / frequency(Q))
        for (t in 1:length(Q)) {
            x <- Q[t] - R_target[t] + S[t]
            if (x < 0) {
                S[t + 1] <- 0
                R[t] <- Q[t] + S[t]
            } else {
                if (x > mid_storage) {
                  S[t + 1] <- mid_storage
                  R[t] <- R_target[t]
                  Spill[t] <- Q[t] - R_target[t] + S[t] - mid_storage
                } else {
                  S[t + 1] <- x
                  R[t] <- R_target[t]
                }
            }
        }
        
        reliability_0 <- 1 - (sum((R / R_target) < 1) / length(Q))
        if (reliability_0 >= reliability) {
            max_storage <- mid_storage
        }
        if (reliability_0 < reliability) {
            min_storage <- mid_storage
        }
        
        if (max_storage - min_storage < 0.01) {
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
        plot(S, ylab = "Storage",main = (paste0("Storage = ", 
            round(mid_storage, 3), " at ", round(reliability_0, 3) * 100, " % reliability")))
        plot(R, ylim = c(0, max(R)), ylab = "Water supplied")
        plot(Spill, ylab = "Spill")
    }
    message(paste0("Converged after ", i,
                 " iterations and time-based reliability equal to ", reliability))
    
    results <- list(mid_storage, S, R, Spill)
    names(results) <- c("Required_storage", "Storage", "Water_supplied", "Spill")
    return(results)
    }
