#' @title Dynamic Programming (Deprecated function; use 'dp_supply' instead)
#' @description Determines the optimal sequence of releases from the reservoir to minimise a penalty cost function based on water supply defict.
#' @param Q             vector or time series object. Net inflows to the reservoir.
#' @param capacity      numerical. The reservoir storage capacity (must be the same volumetric unit as Q and the target release).
#' @param target        numerical. The target release constant.
#' @param S_disc        integer. Storage discretization--the number of equally-sized storage states. Default = 1000.
#' @param R_disc        integer. Release discretization. Default = 10 divisions.
#' @param loss_exp      numeric. The exponent of the penalty cost function--i.e., Cost[t] <- ((target - release[t]) / target) ^ **loss_exp**). Default value is 2.
#' @param S_initial     numeric. The initial storage as a ratio of capacity (0 <= S_initial <= 1). The default value is 1. 
#' @param plot          logical. If TRUE (the default) the storage behavior diagram and release time series are plotted.
#' @param rep_rrv       logical. If TRUE then reliability, resilience and vulnerability metrics are computed and returned.
#' @return Returns the time series of optimal releases and, if requested, the reliability, resilience and vulnerability of the system.
#' @references Loucks, D.P., van Beek, E., Stedinger, J.R., Dijkman, J.P.M. and Villars, M.T. (2005) Water resources systems planning and management: An introduction to methods, models and applications. Unesco publishing, Paris, France.
#' @seealso \code{\link{sdp}} for Stochastic Dynamic Programming
#' @import stats 
#' @export
dp <- function(Q, capacity, target, S_disc = 1000,
               R_disc = 10, loss_exp = 2, S_initial = 1,
               plot = TRUE, rep_rrv = FALSE) {
  
  .Deprecated("dp_supply")
  
  if (is.ts(Q) == FALSE && is.vector(Q) == FALSE) {
    stop("Q must be time series or vector object")
  }
  
  S_states <- seq(from = 0, to = capacity, by = capacity / S_disc)
  R_disc_x <- seq(from = 0, to = target, by = target / R_disc)
  R_costs <- ( ( (target - R_disc_x) / target) ^ loss_exp)
  State_mat <- matrix(0, nrow = length(S_states), ncol = length(R_disc_x))
  State_mat <- apply(State_mat, 2, "+", S_states)
  State_mat <- t(apply(State_mat, 1, "-", R_disc_x))
  Cost_to_go <- vector("numeric", length = length(S_states))
  Bellman <- matrix(0, nrow = length(S_states), ncol = length(Q))
  R_policy <- matrix(0, ncol = length(Q), nrow = length(S_states))
  
  # POLICY OPTIMIZATION----------------------------------------------------------------
  
  for (t in length(Q):1) {
    Cost_mat <- matrix(rep(R_costs, length(S_states)),
                       ncol = length(R_costs), byrow = TRUE)
    Balance_mat <- State_mat + Q[t]
    Cost_mat[which(Balance_mat < 0)] <- NaN
    Balance_mat[which(Balance_mat < 0)] <- NaN
    Balance_mat[which(Balance_mat > capacity)] <- capacity
    Implied_S_state <- round(1 + ((Balance_mat / capacity) *
                                    (length(S_states) - 1)))
    Cost_mat2 <- Cost_mat + matrix(Cost_to_go[Implied_S_state],
                                   nrow = length(S_states))
    Cost_to_go <- apply(Cost_mat2, 1, min, na.rm = TRUE)
    Bellman[, t] <- Cost_to_go
    R_policy[, t] <- apply(Cost_mat2, 1, which.min)
  }
  # ===================================================================================
  
  # POLICY SIMULATION------------------------------------------------------------------
  
  S <- vector("numeric", length(Q) + 1)
  S[1] <- S_initial * capacity
  R <- vector("numeric", length(Q))
  for (t in 1:length(Q)) {
    S_state <- round(1 + ( (S[t] / capacity) *
                             (length(S_states) - 1)))
    R[t] <- R_disc_x[R_policy[S_state, t]]
    if ( (S[t] - R[t] + Q[t]) > capacity) {
      S[t + 1] <- capacity
    } else {
      S[t + 1] <- S[t] - R[t] + Q[t]
    }
  }
  S <- ts(S[2:length(S)], start = start(Q), frequency = frequency(Q))
  R <- ts(R, start = start(Q), frequency = frequency(Q))
  # ===================================================================================
  
  
  # COMPUTE RRV METRICS FROM SIMULATION RESULTS---------------------------------------
  
  if (rep_rrv == TRUE){
    
    deficit <- ts(round(1 - (R / target),5), start = start(Q), frequency = frequency(Q))
    rel_ann <- sum(aggregate(deficit, FUN = mean) == 0) /
      length(aggregate(deficit, FUN = mean))
    rel_time <- sum(deficit == 0) / length(deficit)
    rel_vol <- sum(R) / (target * length(deficit))
    fail.periods <- which(deficit > 0)
    if (length(fail.periods) == 0) {
      resilience <- NA
      vulnerability <- NA
    } else {
      if (length(fail.periods) == 1) {
        resilience <- 1
        vulnerability <- max(deficit)
      } else {
        resilience <- (sum(diff(which(deficit > 0)) > 1) + 1) / (length(which(deficit > 0)))
        fail.refs <- vector("numeric", length = length(fail.periods))
        fail.refs[1] <- 1
        for (j in 2:length(fail.periods)) {
          if (fail.periods[j] > (fail.periods[j - 1] + 1)) {
            fail.refs[j] <- fail.refs[j - 1] + 1
          } else {
            fail.refs[j] <- fail.refs[j - 1]
          }
        }
        n.events <- max(fail.refs)
        event.starts <- by(fail.periods, fail.refs, FUN = min)
        event.ends <- by(fail.periods, fail.refs, FUN = max)
        max.deficits <- vector("numeric", length = n.events)
        for (k in 1:n.events) {
          max.deficits[k] <- max(deficit[event.starts[k]:event.ends[k]])
        }
        
        vulnerability <- mean(max.deficits)
      }
    }
    
    results <- list(S, R ,rel_ann, rel_time, rel_vol, resilience, vulnerability)
    names(results) <- c("storage", "releases", "annual_reliability",
                        "time_based_reliability", "volumetric_reliability",
                        "resilience", "vulnerability")
    
    
    
  } else {
    results <- list(S, R)
    names(results) <- c("storage", "releases")
  }
  
  #===============================================================================
  
  
  if (plot) {
    plot(S, ylab = "storage", ylim = c(0, capacity))
    plot(R, ylab = "release", ylim = c(0, target))
  }
  return(results)
} 