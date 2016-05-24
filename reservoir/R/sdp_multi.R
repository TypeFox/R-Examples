#' @title Stochastic Dynamic Programming with multiple objectives (supply, flood control, amenity)
#' @description Determines the optimal sequence of releases from the reservoir to minimise a penalty cost function based on water supply, spill, and water level. For water supply: Cost[t] = ((target - release[t]) / target) ^ loss_exp[1]). For flood control: Cost[t] = (Spill[t] / quantile(Q, spill_targ)) ^ loss_exp[2]. For amenity: Cost[t] = abs(((storage[t] - (vol_targ * capacity)) / (vol_targ * capacity))) ^ loss_exp[3].  
#' @param Q             time series object. Net inflow to the reservoir.
#' @param capacity      numerical. The reservoir storage capacity (must be the same volumetric unit as Q and the target release).
#' @param target        numerical. The target release constant. Recommended units: Mm^3 (Million cubic meters).
#' @param surface_area  numerical. The reservoir water surface area at maximum capacity. Recommended units: km^2 (square kilometers).
#' @param max_depth     numerical. The maximum water depth of the reservoir at maximum capacity. If omitted, the depth-storage-area relationship will be estimated from surface area and capacity only. Recommended units: meters.
#' @param evap          vector or time series object of length Q, or a numerical constant.  Evaporation from losses from reservoir surface. Varies with level if depth and surface_area parameters are specified. Recommended units: meters, or kg/m2 * 10 ^ -3.
#' @param spill_targ    numerical. The quantile of the inflow time series used to standardise the "minimise spill" objective.
#' @param vol_targ      numerical. The target storage volume constant (as proportion of capacity).
#' @param R_max         numerical. The maximum controlled release.
#' @param weights       vector of length 3 indicating weighting to be applied to release, spill and water level objectives respectively.
#' @param S_disc        integer. Storage discretization--the number of equally-sized storage states. Default = 1000.
#' @param R_disc        integer. Release discretization. Default = 10 divisions.
#' @param Q_disc        vector. Inflow discretization bounding quantiles. Defaults to five inflow classes bounded by quantile vector c(0.0, 0.2375, 0.4750, 0.7125, 0.95, 1.0).
#' @param loss_exp      vector of length 3 indicating the exponents on release, spill and water level deviations from target. Default exponents are c(2,2,2).
#' @param S_initial     numeric. The initial storage as a ratio of capacity (0 <= S_initial <= 1). The default value is 1. 
#' @param plot          logical. If TRUE (the default) the storage behavior diagram and release time series are plotted.
#' @param tol           numerical. The tolerance for policy convergence. The default value is 0.990.
#' @param Markov        logical. If TRUE the current period inflow is used as a hydrological state variable and inflow persistence is incorporated using a first-order, periodic Markov chain. The default is FALSE.
#' @param rep_rrv       logical. If TRUE then reliability, resilience and vulnerability metrics are computed and returned.
#' @return Returns a list that includes: the optimal policy as an array of release decisions dependent on storage state, month/season, and current-period inflow class; the Bellman cost function based on storage state, month/season, and inflow class; the optimized release and storage time series through the training inflow data; the flow discretization (which is required if the output is to be implemented in the rrv function); and, if requested, the reliability, resilience, and vulnerability of the system under the optimized policy. 
#' @seealso \code{\link{dp_multi}} for deterministic Dynamic Programming.
#' @examples layout(1:3)
#' sdp_multi(resX$Q_Mm3, cap = resX$cap_Mm3, target = 0.2 * mean(resX$Q_Mm3))
#' @import stats
#' @importFrom graphics abline lines
#' @export
sdp_multi <- function (Q, capacity, target, surface_area, max_depth, evap,
                       R_max = 2 * target, spill_targ = 0.95, vol_targ = 0.75, Markov = FALSE,
                       weights = c(0.7, 0.2, 0.1), S_disc = 1000, R_disc = 10,
                       Q_disc = c(0.0, 0.2375, 0.4750, 0.7125, 0.95, 1.0),
                       loss_exp = c(2,2,2), S_initial = 1, plot = TRUE, tol = 0.99, rep_rrv = FALSE){
  
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
  if (missing(surface_area)) {
    surface_area <- 0
  }
  
  evap_seas <- as.vector(tapply(evap, cycle(evap), FUN = mean))
  
  # SET UP (non-Markov)---------------------------------------------------------------------------
  
  if (Markov == FALSE){
    Q_month_mat <- matrix(Q, byrow = TRUE, ncol = frq)                                        
    Q.probs <- diff(Q_disc)
    Q_class_med <- apply(Q_month_mat, 2, quantile, type = 8,
                         probs = Q_disc[-1] - (Q.probs / 2))
    S_states <- seq(from = 0, to = capacity, by = capacity / S_disc)                   
    R_disc_x <- seq(from = 0, to = R_max, by = R_max / R_disc)
    Shell.array <- array(0,dim=c(length(S_states),length(R_disc_x),length(Q.probs)))
    #R.star <- aperm(apply(Shell.array, c(1, 3), "+", R_disc_x), c(2, 1, 3))             
    Cost_to_go <- vector("numeric",length=length(S_states))
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
    R_disc_x <- seq(from = 0, to = R_max, by = R_max / R_disc)
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
    Cost_to_go <- matrix(0, nrow = (length(S_states)), ncol = n_Qcl)
    R_policy <- array(0,dim = c(length(S_states), n_Qcl, frq))
    Bellman <- R_policy
    R_policy_test <- R_policy
  }
  

  
  # SET UP (storage-depth-area relationships)----------------------------------------------------- 
  
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
  
  message(paste0("policy converging... (>", tol,")"))
  
  # POLICY OPTIMIZATION----------------------------------------------------------------
  
  if (Markov == FALSE){
    repeat{
      for (t in frq:1){
        R.cstr <- sweep(Shell.array, 3, Q_class_med[,t], "+") +
          sweep(Shell.array, 1, S_states, "+") -
          sweep(Shell.array, 1, evap_seas[t] * S_area_rel / 10 ^ 6, "+")
        R.star <- aperm(apply(Shell.array, c(1, 3), "+", R_disc_x), c(2, 1, 3))
        R.star[,2:(R_disc + 1),][which(R.star[,2:(R_disc + 1),] > R.cstr[,2 : (R_disc + 1),])] <- NaN
        Deficit.arr <- (target - R.star) / target                                           
        Deficit.arr[which(Deficit.arr < 0)] <- 0
        Cost_arr <- ( (abs(Deficit.arr)) ^ loss_exp[1])                          
        S.t_plus_1 <- R.cstr - R.star
        S.t_plus_1[which(S.t_plus_1 < 0)] <- 0
        
        Spill_costs <- S.t_plus_1 - capacity
        Spill_costs[which(Spill_costs < 0)] <- 0
        Spill_costs <- (Spill_costs / quantile(Q, spill_targ)) ^ loss_exp[2]
        
        S.t_plus_1[which(S.t_plus_1 > capacity)] <- capacity
        Vol_costs <- abs(((S.t_plus_1 - vol_targ * capacity) / (vol_targ * capacity))) ^ loss_exp[3]
        Implied_S_state <- round(1 + (S.t_plus_1 / capacity)
                                 * (length(S_states) - 1))
        
        Cost_arr <- weights[1] * Cost_arr + weights[2] * Spill_costs  + weights[3] * Vol_costs
        
        Cost_to_go.arr <- array(Cost_to_go[Implied_S_state],
                                dim = c(length(S_states), length(R_disc_x) , length(Q.probs)))       
        
        Min_cost_arr <- Cost_arr + Cost_to_go.arr
        Min_cost_arr_weighted <- sweep(Min_cost_arr, 3, Q.probs, "*")
        Min_cost_expected <- apply(Min_cost_arr_weighted, c(1, 2), sum)
        
        Bellman[,t] <- Cost_to_go
        Cost_to_go <- apply(Min_cost_expected, 1, min, na.rm = TRUE)
        Results_mat[,t] <- Cost_to_go
        R_policy[,t] <- apply(Min_cost_expected, 1, which.min)
      }
      message(sum(R_policy == R_policy_test) /
                (frq * length(S_states)))   
      if (sum(R_policy == R_policy_test) /
          (frq * length(S_states)) > tol){
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
        R.star[,2:(R_disc + 1),][which(R.star[,2:(R_disc + 1),] > R.cstr[,2 : (R_disc + 1),])] <- NaN
        Deficit.arr <- (target - R.star) / target
        Deficit.arr[which(Deficit.arr < 0)] <- 0
        Cost_arr <- ( (abs(Deficit.arr)) ^ loss_exp[1])                          
        S.t_plus_1 <- R.cstr - R.star
        S.t_plus_1[which(S.t_plus_1 < 0)] <- 0
        
        Spill_costs <- S.t_plus_1 - capacity
        Spill_costs[which(Spill_costs < 0)] <- 0
        Spill_costs <- (Spill_costs / quantile(Q, spill_targ)) ^ loss_exp[2]
        
        S.t_plus_1[which(S.t_plus_1 > capacity)] <- capacity
        Vol_costs <- abs(((S.t_plus_1 - vol_targ * capacity) / (vol_targ * capacity))) ^ loss_exp[3]
        
        Implied_S_state <- round(1 + (S.t_plus_1 / capacity)
                                 * (length(S_states) - 1))
        
        Cost_arr <- weights[1] * Cost_arr + weights[2] * Spill_costs  + weights[3] * Vol_costs
        
        Cost_to_go.arr <- array(Cost_to_go,
                                dim = c(length(S_states), n_Qcl, n_Qcl))       
        Expectation <- apply(sweep(Cost_to_go.arr, c(2,3),
                                   t(Q_trans_probs[,,t]), "*"), c(1,3), sum)
        Exp.arr <- Shell.array
        for (Qt in 1:n_Qcl){
          Exp.arr[,,Qt] <- matrix(Expectation[,Qt][Implied_S_state[,,Qt]], 
                                  ncol = length(R_disc_x))
        }
        R_policy[,,t] <- apply( (Cost_arr + Exp.arr), c(1,3), which.min)
        Cost_to_go <- apply( (Cost_arr + Exp.arr), c(1,3), min, na.rm = TRUE)
        Bellman[,,t] <- Cost_to_go
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
  Spill <- vector("numeric",length(Q))
  
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
      
      R_rec[t_index] <- R
      E[t_index] <- GetEvap(s = S[t_index], q = Qx, r = R, ev = evap[t_index])
      y[t_index] <- GetLevel(c, S[t_index] * 10 ^ 6)
      if ( (S[t_index] - R + Qx - E[t_index]) > capacity) {
        S[t_index + 1] <- capacity
        Spill[t_index] <- S[t_index] - R + Qx - capacity - E[t_index]
      }else{
        if ( (S[t_index] - R + Qx) < 0) {
          S[t_index + 1] <- 0
          R_rec[t_index] <- max(0, S[t_index] + Qx - E[t_index])
        }else{
          S[t_index + 1] <- S[t_index] - R + Qx - E[t_index]
        }
      }
    }
  }
  
  R_policy <- (R_policy - 1) / R_disc
  S <- ts(S[1:(length(S) - 1)],start = start(Q),frequency = frq)
  R_rec <- ts(R_rec, start = start(Q), frequency = frq)
  E <- ts(E, start = start(Q), frequency = frequency(Q))
  y <- ts(y, start = start(Q), frequency = frequency(Q))
  Spill <- ts(Spill, start = start(Q), frequency = frq)
  
  if(plot) {
    plot(R_rec, ylab = "Controlled release", ylim = c(0, R_max)); abline(h = target, lty = 2)
    plot(S, ylab = "Storage", ylim = c(0, capacity)); abline(h = vol_targ * capacity, lty = 2)
    plot(Spill, ylab = "Uncontrolled spill")
  }
  
  total_release_cost <- sum((1 - R_rec/target)[which((R_rec/target) <  1)] ^ loss_exp[1])
  total_spill_cost <- sum((Spill / quantile(Q, spill_targ)) ^ loss_exp[2])
  total_volume_cost <- sum(((S - vol_targ * capacity) / (vol_targ * capacity)) ^ loss_exp[3])
  total_weighted_cost <- weights[1] * total_release_cost + weights[2] * total_spill_cost + weights[3] * total_volume_cost 
  costs <- list(total_release_cost, total_spill_cost, total_volume_cost, total_weighted_cost)
  names(costs) <- c("total_release_cost", "total_spill_cost", "total_volume_cost", "total_weighted_cost")
  
  if (rep_rrv == TRUE){
    
    # COMPUTE RRV METRICS FROM SIMULATION RESULTS---------------------------------------
    
    deficit <- ts(round(1 - (R_rec / target),5), start = start(Q), frequency = frequency(Q))
    rel_ann <- sum(aggregate(deficit, FUN = mean) == 0) /
      length(aggregate(deficit, FUN = mean))
    rel_time <- sum(deficit == 0) / length(deficit)
    rel_vol <- sum(R_rec) / (target * length(deficit))
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
    
    #===============================================================================
    
    results <- list(R_policy, Bellman, S, R_rec, E, y, Spill, rel_ann, rel_time, rel_vol, resilience, vulnerability, Q_disc, costs)
    names(results) <- c("release_policy", "Bellman", "storage", "releases", "evap_loss", "water_level", "spill", "annual_reliability",
                        "time_based_reliability", "volumetric_reliability",
                        "resilience", "vulnerability", "flow_disc", "costs")

  } else {
    results <- list(R_policy, Bellman, S, R_rec, E, y, Spill, Q_disc, costs)
    names(results) <- c("release_policy", "Bellman", "storage", "releases", "evap_loss", "water_level", "spill", "flow_disc", "total_costs")
  }
  
  return(results)
}