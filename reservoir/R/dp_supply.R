#' @title Dynamic Programming for water supply reservoirs
#' @description Determines the optimal sequence of releases from the reservoir to minimise a penalty cost function based on water supply defict.
#' @param Q             vector or time series object. Net inflow totals to the reservoir. Recommended units: Mm^3 (Million cubic meters).
#' @param capacity      numerical. The reservoir storage capacity. Recommended units: Mm^3 (Million cubic meters).
#' @param target        numerical. The target release constant. Recommended units: Mm^3 (Million cubic meters).
#' @param surface_area  numerical. The reservoir water surface area at maximum capacity. Recommended units: km^2 (square kilometers).
#' @param max_depth     numerical. The maximum water depth of the reservoir at maximum capacity. If omitted, the depth-storage-area relationship will be estimated from surface area and capacity only. Recommended units: meters.
#' @param evap          vector or time series object of length Q, or a numerical constant.  Evaporation from losses from reservoir surface. Varies with level if depth and surface_area parameters are specified. Recommended units: meters, or kg/m2 * 10 ^ -3.
#' @param S_disc        integer. Storage discretization--the number of equally-sized storage states. Default = 1000.
#' @param R_disc        integer. Release discretization. Default = 10 divisions.
#' @param loss_exp      numeric. The exponent of the penalty cost function--i.e., Cost[t] <- ((target - release[t]) / target) ^ **loss_exp**). Default value is 2.
#' @param S_initial     numeric. The initial storage as a ratio of capacity (0 <= S_initial <= 1). The default value is 1. 
#' @param plot          logical. If TRUE (the default) the storage behavior diagram and release time series are plotted.
#' @param rep_rrv       logical. If TRUE then reliability, resilience and vulnerability metrics are computed and returned.
#' @return Returns the reservoir simulation output (storage, release, spill), total penalty cost associated with the objective function, and, if requested, the reliability, resilience and vulnerability of the system.
#' @examples layout(1:3)
#' dp_supply(resX$Q_Mm3, capacity = resX$cap_Mm3, target = 0.3 * mean(resX$Q_Mm3), S_disc = 100)
#' @seealso \code{\link{sdp_supply}} for Stochastic Dynamic Programming for water supply reservoirs
#' @import stats 
#' @export
dp_supply <- function(Q, capacity, target, surface_area, max_depth, evap,
                      S_disc = 1000, R_disc = 10, loss_exp = 2, S_initial = 1,
                      plot = TRUE, rep_rrv = FALSE) {
  
  if (is.ts(Q) == FALSE && is.vector(Q) == FALSE) {
    stop("Q must be time series or vector object")
  }

  if (missing(evap)) {
    evap <- rep(0, length(Q))
  }
  if(length(evap) == 1) {
    evap <- rep(evap, length(Q))
  }
  if (length(evap) != length(Q)){
    stop("Evaporation must be either a vector (or time series) length Q, or a single numeric constant")
  }

  if (missing(surface_area)) {
    surface_area <- 0
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
  
  
  if (missing(max_depth)){
    c <- sqrt(2) / 3 * (surface_area * 10 ^ 6) ^ (3/2) / (capacity * 10 ^ 6)
    GetLevel <- function(c, V){
      y <- (6 * V / (c ^ 2)) ^ (1 / 3)
      return(y)
    }
    GetArea <- function(c, V){
      Ay <- (((3 * c * V) / (sqrt(2))) ^ (2 / 3))
      return(Ay)
    }
  } else {
    c <- 2 * capacity / (max_depth * surface_area)
    GetLevel <- function(c, V){
      y <- max_depth * (V / (capacity * 10 ^ 6)) ^ (c / 2)
      return(y)
    }
    GetArea <- function(c, V){
      Ay <- ((2 * (capacity * 10 ^ 6)) / (c * max_depth * (V / (capacity * 10 ^ 6)) ^ (c / 2))) * ((V / (capacity * 10 ^ 6)) ^ (c / 2)) ^ (2 / c)
      Ay[which(is.nan(Ay) == TRUE)] <- 0
      return(Ay)
    }
  }
  
  GetEvap <- function(s, q, r, ev){
    e <- GetArea(c, V = s * 10 ^ 6) * ev / 10 ^ 6
    n <- 0
    repeat{
      n <- n + 1
      s_plus_1 <- max(min(s + q - r - e, capacity), 0)
      e_x <- GetArea(c, V = ((s + s_plus_1) / 2) * 10 ^ 6) * ev / 10 ^ 6
      if (abs(e_x - e) < 0.001 || n > 20){
        break
      } else {
        e <- e_x
      }
    }
    return(e)
  }

  S_area_rel <- GetArea(c, V = S_states * 10 ^ 6)

  # POLICY OPTIMIZATION----------------------------------------------------------------
  
  for (t in length(Q):1) {
    Cost_mat <- matrix(rep(R_costs, length(S_states)),
                       ncol = length(R_costs), byrow = TRUE)
    Balance_mat <- State_mat + Q[t] - (evap[t] * S_area_rel / 10 ^ 6)
    Cost_mat[which(Balance_mat < 0)] <- NaN
    Balance_mat[which(Balance_mat < 0)] <- NaN
    Cost_mat[which(is.nan(Balance_mat[,1]))] <- 0           ## Correction to allow zero release under high evap
    Balance_mat[which(is.nan(Balance_mat[,1]))] <- 0        ## Correction to allow zero release under high evap
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
  E <- vector("numeric", length(Q))
  y <- vector("numeric", length(Q))
  Spill <- vector("numeric", length(Q))
  for (t in 1:length(Q)) {
    S_state <- round(1 + ( (S[t] / capacity) *
                             (length(S_states) - 1)))
    R[t] <- R_disc_x[R_policy[S_state, t]]
    E[t] <- GetEvap(s = S[t], q = Q[t], r = R[t], ev = evap[t])
    
    
    
    y[t] <- GetLevel(c, S[t] * 10 ^ 6)

    if ( (S[t] - R[t] + Q[t] - E[t]) > capacity) {
      S[t + 1] <- capacity
      Spill[t] <- S[t] - R[t] + Q[t] - capacity - E[t]
    } else {
      S[t + 1] <- max(0, S[t] - R[t] + Q[t] - E[t])
    }
  }
  S <- ts(S[1:(length(S) - 1)], start = start(Q), frequency = frequency(Q))
  R <- ts(R, start = start(Q), frequency = frequency(Q))
  E <- ts(E, start = start(Q), frequency = frequency(Q))
  y <- ts(y, start = start(Q), frequency = frequency(Q))
  Spill <- ts(Spill, start = start(Q), frequency = frequency(Q))
  total_penalty <- sum( ( (target - R) / target) ^ loss_exp)
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
    
    results <- list(S, R, E, y, Spill, rel_ann, rel_time, rel_vol, resilience, vulnerability, total_penalty)
    names(results) <- c("storage", "releases", "evap_loss", "water_level",
                        "spill", "annual_reliability",
                        "time_based_reliability", "volumetric_reliability",
                        "resilience", "vulnerability", "total_penalty_cost")
    
    
    
  } else {
    results <- list(S, R, E, y, Spill, total_penalty)
    names(results) <- c("storage", "releases", "evap_loss",
                        "water_level", "spill", "total_penalty_cost")
  }
  
  #===============================================================================
  
  
  if (plot) {
    plot(R, ylab = "Controlled release", ylim = c(0, target))
    plot(S, ylab = "Storage", ylim = c(0, capacity))
    plot(Spill, ylab = "Uncontrolled spill")
  }
  return(results)
}
