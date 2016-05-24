#' @title Simulate a water supply reservoir with specified operating policy.
#' @description Simulates a reservoir for a given inflow time series and assuming Standard Operating Policy (meet target at all times, unless constrained by available water in reservoir plus incoming flows) or an optimised policy deived using \code{\link{sdp_supply}}. 
#' @param Q             vector or time series object. Net inflow totals to the reservoir. Mm^3 (Million cubic meters).
#' @param target        numerical constant, or a time series or vector of the target releases. Must be the same length as Q is given as a vector or time series. Mm^3 (Million cubic meters).
#' @param capacity      numerical. The reservoir capacity. Should be same volumetric unit as Q. Mm^3 (Million cubic meters).
#' @param surface_area  numerical. The reservoir surface area at full capacity. Must be in square kilometers (km^2), or Mm^2.
#' @param max_depth     numerical. The maximum water depth of the reservoir at maximum capacity. Must be in meters. If omitted, the depth-storage-area relationship will be estimated from surface area and capacity only.
#' @param evap          vector or time series object of length Q, or a numerical constant.  Evaporation from losses from reservoir surface. Varies with level if depth and surface_area parameters are specified. Recommended units: meters, or kg/m2 * 10 ^ -3.
#' @param S_initial     numerical. The initial storage as a ratio of capacity (0 <= S_initial <= 1). The default value is 1.
#' @param double_cycle  logical. If TRUE the Q and R time series will be replicated and placed end-to-end to double the simulation. Recommended if the critical period occurs at the end of the sequence.
#' @param plot          logical. If TRUE (the default) the storage and release time series are plotted.
#' @param policy        list. The output of the SDP function. If omitted, Standard Operating Policy is assumed.
#' @return Returns the no-fail storage capacity and corresponding storage behaviour time series.
#' @examples # simulate a reservoir assuming standard operating policy, then compare with SDP-derived policy
#' #trained on historical flows.
#' 
#' # DEFINE RESERVOIR SPECS AND MODEL INPUTS
#' res_cap <- 1500 #Mm3
#' targ <- 150 #Mm3
#' area <- 40 #km2
#' max_d <- 40 #m
#' ev = 0.2 #m
#' Q_pre1980 <- window(resX$Q_Mm3, end = c(1979, 12), frequency = 12)
#' Q_post1980 <- window(resX$Q_Mm3, start = c(1980, 1), frequency = 12)
#' 
#' # SIMULATE WITH SOP
#' layout(1:3)
#' simSOP <- simRes(Q_post1980, capacity = res_cap, target = targ,
#' surface_area = area, max_depth = max_d, evap = ev)
#' 
#' # TRAIN SDP POLICY ON HISTORICAL FLOWS
#' policy_x <- sdp_supply(Q_pre1980, capacity = res_cap, target = targ,
#' surface_area = area, max_depth = max_d, evap = ev, Markov = TRUE, plot = FALSE, S_disc = 100)
#' 
#' # SIMULATE WITH SDP-DERIVED POLICY
#' simSDP <- simRes(Q_post1980, capacity = res_cap, target = targ,
#' surface_area = area, max_depth = max_d, evap = ev, policy = policy_x)
#' @import stats
#' @importFrom graphics abline lines
#' @export
simRes <- function(Q, target, capacity, surface_area, max_depth, evap,
                   double_cycle = FALSE, plot = TRUE, S_initial = 1, policy) {
  
  

  frq <- frequency(Q)
  
  if(length(target)==1){
    target <- rep(target, length(Q))
  }
  if(length(target) != length(Q))
    stop("target must be numerical constant or time series of length Q")
  
  target <- ts(target, start = start(Q), frequency = frequency(Q))
  
  if (double_cycle) {
    Q <- ts(c(Q, Q), start = start(Q), frequency = frq)
    target <- c(target, target)
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
  
  S <- vector("numeric", length = length(Q) + 1); S[1] <- capacity * S_initial
  R <- vector("numeric", length = length(Q))
  E <- vector("numeric", length(Q))
  y <- vector("numeric", length(Q))
  Spill <- vector("numeric", length(Q))
  
  # SIMULATION--------------------------------------------------------------------
  
  if (missing(policy)){
    for (t in 1:length(Q)) {
      
      #x <- Q[t] - target[t] + S[t]
      
      E[t] <- GetEvap(s = S[t], q = Q[t], r = target[t], ev = evap[t])
      y[t] <- GetLevel(c, S[t] * 10 ^ 6)
      
      if ( (S[t] - target[t] + Q[t] - E[t]) > capacity) {
        S[t + 1] <- capacity
        Spill[t] <- S[t] - target[t] + Q[t] - capacity - E[t]
        R[t] <- target[t]
      } else {
        if (S[t] - target[t] + Q[t] - E[t] < 0) {
          S[t + 1] <- 0
          R[t] <- max(0, S[t] + Q[t] - E[t])
        } else {
          S[t + 1] <- S[t] +  Q[t] - target[t] - E[t]
          R[t] <- target[t]
        }
      }
    }
  }
  
  if (!missing(policy)){
    if (is.list(policy) == FALSE) stop("The policy must be the full output returned by the sdp_supply function")
    
    Q_month_mat <- matrix(Q, byrow = TRUE, ncol = frq)
    
    if (length(dim(policy$release_policy)) == 2){
      S_disc <- length(policy$release_policy[,1]) - 1
    } else {
      S_disc <- length(policy$release_policy[,1,1]) - 1
    }
    
    S_states <- seq(from = 0, to = capacity, by = capacity / S_disc)
    Q.probs <- diff(policy$flow_disc)
    Q_class_med <- apply(Q_month_mat, 2, quantile, type = 8,
                         probs = policy$flow_disc[-1] - (Q.probs / 2))
    
    for (yr in 1:nrow(Q_month_mat)) {
      for (month in 1:frq) {
        t_index <- (frq * (yr - 1)) + month   
        S_state <- which.min(abs(S_states - S[t_index]))
        Qx <- Q_month_mat[yr, month]
        
        if (length(dim(policy$release_policy)) == 2){
          Rx <- target[t_index] * policy$release_policy[S_state, month]
        } else {
          Q_class <- which.min(abs(as.vector(Q_class_med[,month] - Qx)))
          Rx <- target[t_index] * policy$release_policy[S_state, Q_class, month]
        }
        
        E[t_index] <- GetEvap(s = S[t_index], q = Qx, r = R, ev = evap[t_index])
        y[t_index] <- GetLevel(c, S[t_index] * 10 ^ 6)
        
        if ( (S[t_index] - Rx + Qx - E[t_index]) > capacity) {
          S[t_index + 1] <- capacity
          R[t_index] <- Rx
          Spill[t_index] <- S[t_index] - Rx + Qx - capacity - E[t_index]
        }else{
          if ( (S[t_index] - Rx + Qx - E[t_index]) < 0) {
            S[t_index + 1] <- 0
            R[t_index] <- max(0, S[t_index] + Qx - E[t_index])
          }else{
            S[t_index + 1] <- S[t_index] - Rx + Qx - E[t_index]
            R[t_index] <- Rx
          }
        }
      }
    }
  }
  
  S <- ts(S[1:length(S) - 1], start = start(Q), frequency = frequency(Q))
  R <- ts(R, start = start(Q), frequency = frequency(Q))
  E <- ts(E, start = start(Q), frequency = frequency(Q))
  y <- ts(y, start = start(Q), frequency = frequency(Q))
  Spill <- ts(Spill, start = start(Q), frequency = frequency(Q))
  
  results <- list(S, R, E, y, Spill)
  names(results) <- c("storage", "releases", "evaporation", "water_level", "spill")
  
  
  if (plot) {
    plot(R, ylab = "Release", ylim = c(0, max(target))); lines(target, lty = 3)
    plot(S, ylab = "Storage", ylim = c(0, capacity))
    plot(Spill, ylab = "Spill")
  }
  
  return(results)

}
