#' @title Dynamic Programming for hydropower reservoirs
#' @description Determines the optimal sequence of turbined releases to maximise the total energy produced by the reservoir.
#' @param Q             time series object. Net inflows to the reservoir. Must be in volumetric units of Mm^3.
#' @param capacity      numerical. The total reservoir storage capacity (including unusable "dead" storage). Must be in Mm^3.
#' @param capacity_live numerical. The volume of usable water in the reservoir ("live capacity" or "active storage"). capacity_live <= capacity. Default capacity_live = capacity. Must be in Mm^3.
#' @param surface_area  numerical. The reservoir surface area at full capacity. Must be in square kilometers (km^2), or Mm^2.
#' @param max_depth     numerical. The maximum water depth of the reservoir at maximum capacity. If omitted, the depth-storage-area relationship will be estimated from surface area and capacity only. Recommended units: meters.
#' @param evap          vector or time series object of length Q, or a numerical constant, representing evaporation loss potential from reservoir surface. Varies with level if depth and surface_area parameters are specified. Must be in meters, or kg/m2 * 10 ^ -3.
#' @param installed_cap numerical. The hydropower plant electric capacity (MW).
#' @param efficiency    numerical. The hydropower plant efficiency. Default is 0.9, but, unless user specifies an efficiency, it will be automatically re-estimated if head and qmax are supplied.
#' @param head          numerical. The maximum hydraulic head of the hydropower plant (m). Can be omitted and estimated if qmax is supplied.
#' @param qmax          numerical. The maximum flow into the hydropower plant. Can be omitted and estimated if head is supplied. Must be in volumetric units of Mm^3.
#' @param S_disc        integer. Storage discretization--the number of equally-sized storage states. Default = 1000.
#' @param R_disc        integer. Release discretization. Default = 10 divisions.
#' @param S_initial     numeric. The initial storage as a ratio of capacity (0 <= S_initial <= 1). The default value is 1. 
#' @param plot          logical. If TRUE (the default) the storage behavior diagram and release time series are plotted.
#' @return Returns the time series of optimal releases and simulated storage, evaporation, depth, uncontrolled spill, and power generated. Total energy generated is also returned.
#' @examples layout(1:4)
#' dp_hydro(resX$Q_Mm3, resX$cap_Mm3, surface_area = resX$A_km2,
#' installed_cap = resX$Inst_cap_MW, qmax = mean(resX$Q_Mm3), S_disc = 100)
#' @seealso \code{\link{sdp_hydro}} for Stochastic Dynamic Programming for hydropower reservoirs.
#' @import stats 
#' @export
dp_hydro <- function(Q, capacity, capacity_live = capacity, surface_area, evap,
                     installed_cap, head, qmax, max_depth, efficiency = 0.9,
                     S_disc = 1000, R_disc = 10, S_initial = 1, plot = TRUE) {
  
  if (is.ts(Q) == FALSE) {
    stop("Q must be time series object")
  }
  if ((missing(head) || is.na(head)) && (missing(qmax) || is.na(qmax))) {
    stop("You must enter a value for either head or qmax")
  }
  if (!missing(head) && !missing(qmax) && missing(efficiency) && !is.na(head) && !is.na(qmax)) {
    efficiency <- installed_cap / (9.81 * 1000 * head * (qmax / ((365.25/frequency(Q)) * 24 * 60 * 60)))
    if (efficiency > 1) {
      warning("Check head, qmax and installed_cap: calculated efficiency exceeds 100 %")
    }
  }
  if (missing(head) || is.na(head)) {
    head <- installed_cap / (efficiency * 9.81 * 1000 * (qmax / ((365.25/frequency(Q)) * 24 * 60 * 60)))
  }
  if (missing(qmax) || is.na(qmax)) {
    qmax <- (installed_cap / (efficiency * 9.81 * 1000 * head)) * ((365.25/frequency(Q)) * 24 * 60 * 60)
  }
  
  if (missing(evap)) {
    evap <- rep(0, length(Q))
  }
  if(length(evap) == 1) {
    evap <- ts(rep(evap, length(Q)), start = start(Q), frequency = frequency(Q))
  }
  if (length(evap) != length(Q)){
    stop("Evaporation must be either a vector (or time series) length Q, or a single numeric constant")
  }
  
  S_states <- seq(from = 0, to = capacity, by = capacity / S_disc)
  R_disc_x <- seq(from = 0, to = qmax, by = qmax / R_disc)
  State_mat <- matrix(0, nrow = length(S_states), ncol = length(R_disc_x))
  State_mat <- apply(State_mat, 2, "+", S_states)
  State_mat <- t(apply(State_mat, 1, "-", R_disc_x))
  Rev_to_go <- vector("numeric", length = length(S_states))
  Bellman <- matrix(0, nrow = length(S_states), ncol = length(Q))
  R_policy <- matrix(0, ncol = length(Q), nrow = length(S_states))
  
  if (missing(max_depth) || is.na(max_depth)){
    c <- sqrt(2) / 3 * (surface_area * 10 ^ 6) ^ (3/2) / (capacity * 10 ^ 6)
    GetLevel <- function(c, V){
      y <- (6 * V / (c ^ 2)) ^ (1 / 3)
      return(y)
    }
    GetArea <- function(c, V){
      Ay <- (((3 * c * V) / (sqrt(2))) ^ (2 / 3))
      return(Ay)
    }
    yconst <- head - GetLevel(c, capacity * 10 ^ 6)
  } else {
    c <- 2 * (capacity * 10 ^ 6) / (max_depth * surface_area * 10 ^ 6)
    GetLevel <- function(c, V){
      y <- max_depth * (V / (capacity * 10 ^ 6)) ^ (c / 2)
      return(y)
    }
    GetArea <- function(c, V){
      Ay <- ((2 * (capacity * 10 ^ 6)) / (c * max_depth * (V / (capacity * 10 ^ 6)) ^ (c / 2))) * ((V / (capacity * 10 ^ 6)) ^ (c / 2)) ^ (2 / c)
      Ay[which(is.nan(Ay) == TRUE)] <- 0
      return(Ay)
    }
    yconst <- head - max_depth
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
    Release_mat <- matrix(rep(R_disc_x, length(S_states)),
                          ncol = length(R_disc_x), byrow = TRUE)
    Balance_mat <- State_mat + Q[t] - (evap[t] * S_area_rel / 10 ^ 6)
    Release_mat[, 2:ncol(Balance_mat)][which(Balance_mat[, 2:ncol(Balance_mat)] < (capacity - capacity_live))] <- NaN
    Balance_mat[, 2:ncol(Balance_mat)][which(Balance_mat[, 2:ncol(Balance_mat)] < (capacity - capacity_live))] <- NaN
    Balance_mat[which(Balance_mat < 0)] <- 0
    Balance_mat[which(Balance_mat > capacity)] <- capacity
    Implied_S_state <- round(1 + ((Balance_mat / capacity) * (length(S_states) - 1)))
    H_mat <- GetLevel(c, ((Balance_mat + S_states) * (10 ^ 6))  / 2) + yconst
    Rev_mat <- Release_mat * H_mat
    Rev_mat2 <- Rev_mat + matrix(Rev_to_go[Implied_S_state],
                                 nrow = length(S_states))
    Rev_to_go <- apply(Rev_mat2, 1, max, na.rm = TRUE)
    Bellman[, t] <- Rev_to_go
    R_policy[, t] <- apply(Rev_mat2, 1, which.max)
  }
  # ===================================================================================
  
  # POLICY SIMULATION------------------------------------------------------------------
  
  S <- vector("numeric", length(Q) + 1)
  S[1] <- S_initial * capacity
  R <- vector("numeric", length(Q))
  E <- vector("numeric", length(Q))
  y <- vector("numeric", length(Q))
  Spill <- vector("numeric", length(Q))
  Power <- vector("numeric", length(Q))
  for (t in 1:length(Q)) {
    S_state <- round(1 + ( (S[t] / capacity) * (length(S_states) - 1)))
    R[t] <- R_disc_x[R_policy[S_state, t]]
    E[t] <- GetEvap(s = S[t], q = Q[t], r = R[t], ev = evap[t])
    y[t] <- GetLevel(c, S[t] * (10 ^ 6))
    
    if ( (S[t] - R[t] + Q[t] - E[t]) > capacity) {
      S[t + 1] <- capacity
      Spill[t] <- S[t] - R[t] + Q[t] - E[t] - capacity
    } else {
      S[t + 1] <- max(0, S[t] - R[t] + Q[t] - E[t])
    }
    Power[t] <- efficiency * 1000 * 9.81 * (GetLevel(c, mean(S[t:(t+1)]) * (10 ^ 6)) + yconst) * R[t] / (365.25 / frequency(Q) * 24 * 60 * 60)
  }

  S <- ts(S[1:(length(S) - 1)], start = start(Q), frequency = frequency(Q))
  R <- ts(R, start = start(Q), frequency = frequency(Q))
  E <- ts(E, start = start(Q), frequency = frequency(Q))
  y <- ts(y, start = start(Q), frequency = frequency(Q))
  Spill <- ts(Spill, start = start(Q), frequency = frequency(Q))
  Power <- ts(Power, start = start(Q), frequency = frequency(Q))
  Energy_MWh <- sum(Power * (365.25 / frequency(Q)) * 24)
  results <- list(S, R, E, y, Spill, Power, Energy_MWh)
  names(results) <- c("Storage_Mm3", "Release_Mm3", "Evap_loss_Mm3", "Water_level_m", "Spill_Mm3", "Power_MW", "Total_Energy_MWh")
  
  # ===================================================================================
  
  if (plot) {
    plot(R, ylab = "Turbined release [Mm3]", ylim = c(0, qmax), main = paste0("Total output = ", round(Energy_MWh/1000000, 3), " TWh"))
    plot(Power, ylab = "Power [MW]", ylim = c(0, installed_cap))
    plot(S, ylab = "Storage [Mm3]", ylim = c(0, capacity))
    plot(Spill, ylab = "Uncontrolled spill [Mm3]")
  }
  return(results)
}