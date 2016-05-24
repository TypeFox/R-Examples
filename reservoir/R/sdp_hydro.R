#' @title Stochastic Dynamic Programming for hydropower reservoirs
#' @description Determines the optimal policy of turbined releases to maximise the total energy produced by the reservoir. The policy can be based on season and storage level, or season, storage level, and current-period inflow.
#' @param Q             time series object. Net inflows to the reservoir. Must be in volumetric units of Mm^3.
#' @param capacity      numerical. The total reservoir storage capacity (including unusable "dead" storage). Must be in Mm^3.
#' @param capacity_live numerical. The volume of usable water in the reservoir ("live capacity" or "active storage"). capacity_live <= capacity. Default capacity_live = capacity. Must be in Mm^3.
#' @param surface_area  numerical. The reservoir surface area at full capacity. Must be in square kilometers (km^2), or Mm^2.
#' @param max_depth     numerical. The maximum water depth of the reservoir at maximum capacity. If omitted, the depth-storage-area relationship will be estimated from surface area and capacity only. Recommended units: meters.
#' @param evap          vector or time series object of length Q, or a numerical constant, representing evaporation loss potential from reservoir surface. Varies with level if depth and surface_area parameters are specified. Must be in meters, or kg/m2 * 10 ^ -3.
#' @param installed_cap numerical. The hydropower plant electric capacity (MW).
#' @param efficiency    numerical. The hydropower plant efficiency. Default is 0.9, but, unless user specifies an efficiency, it will be automatically re-estimated if head and qmax are supplied.
#' @param head          numerical. The maximum hydraulic head of the hydropower plant (m). Can be omitted if qmax is supplied.
#' @param qmax          numerical. The maximum flow into the hydropower plant. Can be omitted and estimated if head is supplied. Must be in volumetric units of Mm^3.
#' @param S_disc        integer. Storage discretization--the number of equally-sized storage states. Default = 1000.
#' @param R_disc        integer. Release discretization. Default = 10 divisions.
#' @param Q_disc        vector. Inflow discretization bounding quantiles. Defaults to five inflow classes bounded by quantile vector c(0.0, 0.2375, 0.4750, 0.7125, 0.95, 1.0).
#' @param S_initial     numeric. The initial storage as a ratio of capacity (0 <= S_initial <= 1). The default value is 1. 
#' @param plot          logical. If TRUE (the default) the storage behavior diagram and release time series are plotted.
#' @param tol           numerical. The tolerance for policy convergence. The default value is 0.990.
#' @param Markov        logical. If TRUE the current period inflow is used as a hydrological state variable and inflow persistence is incorporated using a first-order, periodic Markov chain. The default is FALSE.
#' @return Returns the optimal release policy, associated Bellman function, simulated storage, release, evaporation, depth, uncontrolled spill, and power generated, and total energy generated.
#' @seealso \code{\link{dp_hydro}} for deterministic Dynamic Programming for hydropower reservoirs.
#' @examples layout(1:4)
#' sdp_hydro(resX$Q_Mm3, resX$cap_Mm3, surface_area = resX$A_km2,
#' installed_cap = resX$Inst_cap_MW, qmax = mean(resX$Q_Mm3))
#' sdp_hydro(resX$Q_Mm3, resX$cap_Mm3, surface_area = resX$A_km2,
#' installed_cap = resX$Inst_cap_MW, qmax = mean(resX$Q_Mm3), Markov = TRUE)
#' @import stats
#' @export
sdp_hydro <- function (Q, capacity, capacity_live = capacity,
                       surface_area, max_depth, evap, installed_cap, head, qmax,
                       efficiency = 0.9, S_disc = 1000, R_disc = 10,
                       Q_disc = c(0.0, 0.2375, 0.4750, 0.7125, 0.95, 1.0),
                       S_initial = 1, plot = TRUE, tol = 0.99, Markov = FALSE){
  
  frq <- frequency(Q)
  if (is.ts(Q)==FALSE) stop("Q must be seasonal time series object with frequency of 12 or 4")
  if (frq != 12 && frq != 4) stop("Q must have frequency of 4 or 12")  
  if (missing(evap)) {
    evap <- ts(rep(0, length(Q)), start = start(Q), frequency = frq)
  }
  if(length(evap) == 1) {
    evap <- ts(rep(evap, length(Q)), start = start(Q), frequency = frq)
  }
  if (length(evap) != length(Q) && length(evap) != frq){
    stop("Evaporation must be either a time series of length Q, a vector of length frequency(Q), or a single numeric constant")
  }
  if (start(Q)[2] != 1){
    message("NOTE: First incomplete year of time series removed")
    Q <- window(Q, start = c(start(Q)[1] + 1, 1), frequency = frq)
  }
  if(end(Q)[2] != frq){
    message("NOTE: Final incomplete year of time series removed")
    Q <- window(Q, end = c(end(Q)[1] - 1, frq), frequency = frq)
  }
  if (length(evap) == frq){
    evap <- ts(rep(evap, length(Q) / frq), start = start(Q), frequency = frq)
  } else {
    if(is.ts(evap)==FALSE) stop("Evaporation must be either a time series of length Q or a vector of length frequency(Q) for a seasonal evaporation profile")
    evap <- window(evap, start = start(Q), end = end(Q), frequency = frq)
  }
  if ((missing(head) || is.na(head)) && (missing(qmax) || is.na(qmax))) {
    stop("You must enter a value for either head or qmax")
  }
  if (!missing(head) && !missing(qmax) && missing(efficiency) && !is.na(head) && !is.na(qmax)) {
    efficiency <- installed_cap / (9.81 * 1000 * head * (qmax / ((365.25/frq) * 24 * 60 * 60)))
    if (efficiency > 1) {
      warning("Check head, qmax and installed_cap: calculated efficiency exceeds 100 %")
    }
  }
  if (missing(head) || is.na(head)) {
    head <- installed_cap / (efficiency * 9.81 * 1000 * (qmax / ((365.25/frq) * 24 * 60 * 60)))
  }
  if (missing(qmax) || is.na(qmax)){
    qmax <- (installed_cap / (efficiency * 9.81 * 1000 * head)) * ((365.25/frq) * 24 * 60 * 60)
  }
  
  evap_seas <- as.vector(tapply(evap, cycle(evap), FUN = mean))

  # SET UP (non-Markov)---------------------------------------------------------------------------
  
  if (Markov == FALSE){
    Q_month_mat <- matrix(Q, byrow = TRUE, ncol = frq)                                        
    Q.probs <- diff(Q_disc)
    Q_class_med <- apply(Q_month_mat, 2, quantile, type = 8,
                         probs = Q_disc[-1] - (Q.probs / 2))
    S_states <- seq(from = 0, to = capacity, by = capacity / S_disc)                   
    R_disc_x <- seq(from = 0, to = qmax, by = qmax / R_disc)
    Shell.array <- array(0,dim=c(length(S_states),length(R_disc_x),length(Q.probs)))
    #R.star <- aperm(apply(Shell.array, c(1, 3), "+", R_disc_x), c(2, 1, 3))             
    Rev_to_go <- vector("numeric",length=length(S_states))
    Results_mat <- matrix(0,nrow=length(S_states),ncol=frq)
    R_policy <- matrix(0,nrow=length(S_states),ncol=frq)
    Bellman <- R_policy
    R_policy_test <- R_policy
    
    
    # SET UP (Markov)-------------------------------------------------------------------------------  
    
  } else if (Markov == TRUE){
    Q_month_mat <- matrix(Q, byrow = TRUE, ncol = frq)                                        
    n_Qcl <- length(Q_disc) - 1
    Q.probs <- diff(Q_disc)
    Q_class_med <- apply(Q_month_mat, 2, quantile, type = 8,
                         probs = Q_disc[-1] - (Q.probs / 2))
    S_states <- seq(from = 0, to = capacity, by = capacity / S_disc)                   
    R_disc_x <- seq(from = 0, to = qmax, by = qmax / R_disc)
    Shell.array <- array(0, dim = c(length(S_states), length(R_disc_x),
                                    length(Q.probs)))
    #R.star <- aperm(apply(Shell.array, c(1, 3), "+", R_disc_x), c(2, 1, 3))             
    Q_class.mat <- matrix(nrow=length(Q_month_mat[,1]),ncol=frq)
    for (m in 1:frq){
      Q_disc_x <- gtools::quantcut(Q_month_mat[,m], Q_disc)
      Q_class.mat[,m] <- as.numeric(as.vector(factor(Q_disc_x,
                                                     labels = c(1:n_Qcl))))
    }
    Q_trans_probs <- array(0, c(length(Q_disc) - 1, length(Q_disc) - 1, frq))             
    for (m in 1 : frq){
      for (cl in 1 : n_Qcl){
        if (m == frq){
          Tr.count <- table(factor(Q_class.mat[which(Q_class.mat[1:(length(Q_month_mat[,1]) - 1),
                                                                 frq] == cl) + 1, 1], 1:n_Qcl))
        }else{
          Tr.count <- table(factor(Q_class.mat[which(Q_class.mat[,m] == cl),
                                               m + 1], 1:n_Qcl)) 
        }
        Tr.freq <-  Tr.count / sum(Tr.count)
        Q_trans_probs[cl,,m] <- Tr.freq
      }}
    Rev_to_go <- matrix(0, nrow = (length(S_states)), ncol = n_Qcl)
    R_policy <- array(0,dim = c(length(S_states), n_Qcl, frq))
    Bellman <- R_policy
    R_policy_test <- R_policy
  }
  
  
  # SET UP (storage-depth-area relationships)----------------------------------------------------- 
  
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
    if (yconst <0){
      capacity_live <- min(capacity_live, capacity - (-yconst) ^ 3 * c ^ 2 / 6 / 10 ^ 6)
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
    yconst <- head - max_depth
    if (yconst <0){
      capacity_live <- min(capacity_live, capacity - (-yconst / max_depth) ^ (2 / c) * capacity )
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
  
  
  message(paste0("policy converging... (>", tol,")"))
  
  # POLICY OPTIMIZATION (non-Markov)-------------------------------------------------------------
  
  if (Markov == FALSE){
    repeat{
      for (t in frq:1){
        R.cstr <- sweep(Shell.array, 3, Q_class_med[,t], "+") +
          sweep(Shell.array, 1, S_states, "+") - 
          sweep(Shell.array, 1, evap_seas[t] * S_area_rel / 10 ^ 6, "+")
        R.star <- aperm(apply(Shell.array, c(1, 3), "+", R_disc_x), c(2, 1, 3))
        R.star[,2:(R_disc + 1),][which(R.star[,2:(R_disc + 1),] > R.cstr[,2 : (R_disc + 1),] - (capacity - capacity_live))] <- NaN
        S.t_plus_1 <- R.cstr - R.star
        S.t_plus_1[which(S.t_plus_1 < 0)] <- 0
        S.t_plus_1[which(S.t_plus_1 > capacity)] <- capacity
        
        H_arr <- GetLevel(c, ((S.t_plus_1 + S_states) * (10 ^ 6))  / 2) + yconst
        Rev_arr <- R.star * H_arr
        Implied_S_state <- round(1 + (S.t_plus_1 / capacity)
                                 * (length(S_states) - 1))
        Rev_to_go.arr <- array(Rev_to_go[Implied_S_state],
                                dim = c(length(S_states), length(R_disc_x) , length(Q.probs)))       
        Max_rev_arr <- Rev_arr + Rev_to_go.arr
        Max_rev_arr_weighted <- sweep(Max_rev_arr, 3, Q.probs, "*")
        Max_rev_expected <- apply(Max_rev_arr_weighted, c(1, 2), sum)
        Bellman[,t] <- Rev_to_go
        Rev_to_go <- apply(Max_rev_expected, 1, max, na.rm = TRUE)
        Results_mat[,t] <- Rev_to_go
        R_policy[,t] <- apply(Max_rev_expected, 1, which.max)
      }
      message(sum(R_policy == R_policy_test) / (frq * length(S_states)))   
      if (sum(R_policy == R_policy_test) / (frq * length(S_states)) > tol){
        break
      }
      R_policy_test <- R_policy
    }
    
    # POLICY OPTIMIZATION (Markov)-------------------------------------------------------------
    
  } else  if (Markov == TRUE){
    repeat{
      for (t in frq:1){
        R.cstr <- sweep(Shell.array, 3, Q_class_med[,t], "+") +
          sweep(Shell.array, 1, S_states, "+") -
          sweep(Shell.array, 1, evap_seas[t] * S_area_rel / 10 ^ 6, "+")
        R.star <- aperm(apply(Shell.array, c(1, 3), "+", R_disc_x), c(2, 1, 3))
        R.star[,2:(R_disc + 1),][which(R.star[,2:(R_disc + 1),] > R.cstr[,2 : (R_disc + 1),] - (capacity - capacity_live))] <- NaN
        S.t_plus_1 <- R.cstr - R.star
        S.t_plus_1[which(S.t_plus_1 < 0)] <- 0
        S.t_plus_1[which(S.t_plus_1 > capacity)] <- capacity
        
        H_arr <- GetLevel(c, ((S.t_plus_1 + S_states) * (10 ^ 6))  / 2) + yconst
        Rev_arr <- R.star * H_arr
        
        Implied_S_state <- round(1 + (S.t_plus_1 / capacity)
                                 * (length(S_states) - 1))
        Rev_to_go.arr <- array(Rev_to_go,
                                dim = c(length(S_states), n_Qcl, n_Qcl))       
        Expectation <- apply(sweep(Rev_to_go.arr, c(2,3),
                                   t(Q_trans_probs[,,t]), "*"), c(1,3), sum)
        Exp.arr <- Shell.array
        for (Qt in 1:n_Qcl){
          Exp.arr[,,Qt] <- matrix(Expectation[,Qt][Implied_S_state[,,Qt]], 
                                  ncol = length(R_disc_x))
        }
        R_policy[,,t] <- apply( (Rev_arr + Exp.arr), c(1,3), which.max)
        Rev_to_go <- apply( (Rev_arr + Exp.arr), c(1,3), max, na.rm = TRUE)
        Bellman[,,t] <- Rev_to_go
      }
      message(sum(R_policy == R_policy_test) / (frq * length(S_states) * n_Qcl))   
      if (sum(R_policy == R_policy_test) / (frq * length(S_states) * n_Qcl) > tol){
        break
      }
      R_policy_test <- R_policy
    }
  }
  
  # ===================================================================================
  
  # POLICY SIMULATION------------------------------------------------------------------
  
  S <- vector("numeric",length(Q) + 1); S[1] <- S_initial * capacity    
  R_rec <- vector("numeric",length(Q))
  E <- vector("numeric", length(Q))
  y <- vector("numeric", length(Q))
  Spill <- vector("numeric", length(Q))
  Power <- vector("numeric", length(Q))
  for (yr in 1:nrow(Q_month_mat)) {
    for (month in 1:frq) {
      t_index <- (frq * (yr - 1)) + month   
      S_state <- which.min(abs(S_states - S[t_index]))
      Qx <- Q_month_mat[yr,month]
      if (Markov == FALSE){
        R <- R_disc_x[R_policy[S_state,month]]
      } else if (Markov == TRUE){
        Q_class <- which.min(abs(as.vector(Q_class_med[,month] - Qx)))
        R <- R_disc_x[R_policy[S_state,Q_class,month]]
      }
      R <- min(R, S[t_index] + Qx - (capacity - capacity_live))
      R_rec[t_index] <- R
      E[t_index] <- GetEvap(s = S[t_index], q = Qx, r = R, ev = evap[t_index])
      y[t_index] <- GetLevel(c, S[t_index] * 10 ^ 6)
      
      
      if ( (S[t_index] - R + Qx - E[t_index]) > capacity) {
        S[t_index + 1] <- capacity
        Spill[t_index] <- S[t_index] - R + Qx - capacity - E[t_index]
      }else{
        if ( (S[t_index] - R + Qx - E[t_index]) < 0) {
          S[t_index + 1] <- 0
          R_rec[t_index] <- max(0, S[t_index] + Qx - E[t_index])
        }else{
          S[t_index + 1] <- S[t_index] - R + Qx - E[t_index]
        }
      }
      Power[t_index] <- max(efficiency * 1000 * 9.81 * (GetLevel(c,mean(S[t_index:(t_index + 1)]) * (10 ^ 6)) + yconst) * 
                              R_rec[t_index] / (365.25 / frq * 24 * 60 * 60), 0)
    }
  }
  R_policy <- (R_policy - 1) / R_disc
  S <- ts(S[1:(length(S) - 1)],start = start(Q),frequency = frq)
  R_rec <- ts(R_rec, start = start(Q), frequency = frq)
  E <- ts(E, start = start(Q), frequency = frq)
  y <- ts(y, start = start(Q), frequency = frq)
  Spill <- ts(Spill, start = start(Q), frequency = frq)
  Power <- ts(Power, start = start(Q), frequency = frq)
  Energy_MWh <- sum(Power * (365.25 / frq) * 24)
  
  if(plot) {
    plot(R_rec, ylab = "Turbined release [Mm3]", ylim = c(0, qmax), main = paste0("Total output = ", round(Energy_MWh/1000000, 3), " TWh"))
    plot(Power, ylab = "Power [MW]", ylim = c(0, installed_cap))
    plot(S, ylab = "Storage [Mm3]", ylim = c(0, capacity))
    plot(Spill, ylab = "Uncontrolled spill [Mm3]")
  }
  

  results <- list(R_policy, Bellman, S, R_rec, E, y, Spill, Power, Q_disc, Energy_MWh)
  names(results) <- c("release_policy", "Bellman", "storage",
                      "releases", "evap_loss", "water_level",
                      "spill", "power", "flow_disc", "Energy_MWh")
  
  
  return(results)
}