#' @title Stochastic Dynamic Programming (Deprecated function; use 'sdp_supply' instead)
#' @description Derives the optimal release policy based on storage state, inflow class and within-year period.
#' @param Q             time series object. Net inflows to the reservoir.
#' @param capacity      numerical. The reservoir storage capacity (must be the same volumetric unit as Q and the target release).
#' @param target        numerical. The target release constant.  
#' @param S_disc        integer. Storage discretization--the number of equally-sized storage states. Default = 1000.
#' @param R_disc        integer. Release discretization. Default = 10 divisions.
#' @param Q_disc        vector. Inflow discretization bounding quantiles. Defaults to five inflow classes bounded by quantile vector c(0.0, 0.2375, 0.4750, 0.7125, 0.95, 1.0).
#' @param loss_exp      numeric. The exponent of the penalty cost function--i.e., Cost[t] <- ((target - release[t]) / target) ^ **loss_exp**). Default value is 2.
#' @param S_initial     numeric. The initial storage as a ratio of capacity (0 <= S_initial <= 1). The default value is 1. 
#' @param plot          logical. If TRUE (the default) the storage behavior diagram and release time series are plotted.
#' @param tol           numerical. The tolerance for policy convergence. The default value is 0.990.
#' @param rep_rrv       logical. If TRUE then reliability, resilience and vulnerability metrics are computed and returned.
#' @return Returns a list that includes: the optimal policy as an array of release decisions dependent on storage state, month/season, and current-period inflow class; the Bellman cost function based on storage state, month/season, and inflow class; the optimized release and storage time series through the training inflow data; the flow discretization (which is required if the output is to be implemented in the rrv function); and, if requested, the reliability, resilience, and vulnerability of the system under the optimized policy. 
#' @references Loucks, D.P., van Beek, E., Stedinger, J.R., Dijkman, J.P.M. and Villars, M.T. (2005) Water resources systems planning and management: An introduction to methods, models and applications. Unesco publishing, Paris, France.
#' @references Gregory R. Warnes, Ben Bolker and Thomas Lumley (2014). gtools: Various R programming tools. R package version 3.4.1. http://CRAN.R-project.org/package=gtools
#' @seealso \code{\link{sdp}} for deterministic Dynamic Programming
#' @import stats
#' @export
sdp <- function (Q, capacity, target, S_disc = 1000, R_disc = 10,
                 Q_disc = c(0.0, 0.2375, 0.4750, 0.7125, 0.95, 1.0),
                 loss_exp = 2, S_initial = 1, plot = TRUE, tol = 0.99, rep_rrv = FALSE){
  
  .Deprecated("dp_supply")
  
  frq <- frequency(Q)
  if (is.ts(Q)==FALSE) stop("Q must be seasonal time series object with frequency of 12 or 4")
  if (frq != 12 && frq != 4) stop("Q must have frequency of 4 or 12")  
  if (start(Q)[2] != 1){
    message("NOTE: First incomplete year of time series removed")
    Q <- window(Q, start = c(start(Q)[1] + 1, 1), frequency = frq)
  }
  if(end(Q)[2] != frq){
    message("NOTE: Final incomplete year of time series removed")
    Q <- window(Q, end = c(end(Q)[1] - 1, frq), frequency = frq)
  }
  Q_month_mat <- matrix(Q, byrow = TRUE, ncol = frq)                                        
  n_Qcl <- length(Q_disc) - 1
  Q.probs <- diff(Q_disc)
  Q_class_med <- apply(Q_month_mat, 2, quantile, type = 8,
                       probs = Q_disc[-1] - (Q.probs / 2))
  S_states <- seq(from = 0, to = capacity, by = capacity / S_disc)                   
  R_disc_x <- seq(from = 0, to = target, by = target / R_disc)
  Shell.array <- array(0, dim = c(length(S_states), length(R_disc_x),
                                  length(Q.probs)))
  R.star <- aperm(apply(Shell.array, c(1, 3), "+", R_disc_x), c(2, 1, 3))             
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
  message(paste0("policy converging... (>", tol,")"))
  
  # POLICY OPTIMIZATION----------------------------------------------------------------
  
  repeat{
    for (t in frq:1){
      R.cstr <- sweep(Shell.array, 3, Q_class_med[,t], "+") +
        sweep(Shell.array, 1, S_states, "+")                               
      R.star[which(R.star > R.cstr)] <- NaN                              
      Deficit.arr <- (R.star - target) / target                                           
      Cost_arr <- ( (abs(Deficit.arr)) ^ loss_exp)                          
      S.t_plus_1 <- R.cstr - R.star
      Implied_S_state <- round(1 + (S.t_plus_1 / capacity)
                               * (length(S_states) - 1))
      Implied_S_state[which(Implied_S_state > length(S_states))] <- length(S_states)
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
    message(sum(R_policy == R_policy_test) /
              (frq * length(S_states) * n_Qcl))   
    if (sum(R_policy == R_policy_test) /
        (frq * length(S_states) * n_Qcl) > tol){
      break
    }
    
    R_policy_test <- R_policy
  }
  
  # ===================================================================================
  
  # POLICY SIMULATION------------------------------------------------------------------
  
  
  S <- vector("numeric",length(Q) + 1); S[1] <- S_initial * capacity    
  R_rec <- vector("numeric",length(Q))                      
  for (yr in 1:nrow(Q_month_mat)) {
    for (month in 1:frq) {
      t_index <- (frq * (yr - 1)) + month   
      S_state <- which.min(abs(S_states - S[t_index]))
      Qx <- Q_month_mat[yr,month]
      Q_class <- which.min(abs(as.vector(Q_class_med[,month] - Qx)))
      R <- R_disc_x[R_policy[S_state,Q_class,month]]
      R_rec[t_index] <- R
      if ( (S[t_index] - R + Qx) > capacity) {
        S[t_index + 1] <- capacity
      }else{
        if ( (S[t_index] - R + Qx) < 0) {
          S[t_index + 1] <- 0
          R_rec[t_index] <- S[t_index] + Qx
        }else{
          S[t_index + 1] <- S[t_index] - R + Qx
        }
      }
    }
  }
  R_policy <- (R_policy - 1) / (max(R_policy) - 1)
  S <- ts(S[2:length(S)],start = start(Q),frequency = frq)
  R_rec <- ts(R_rec, start = start(Q), frequency = frq)
  if(plot) {
    plot(S, ylab = "storage", ylim = c(0, capacity))
    plot(R_rec, ylab = "release", ylim = c(0, target))
  }
  
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
    
    results <- list(R_policy, Bellman, S, R_rec, rel_ann, rel_time, rel_vol, resilience, vulnerability, Q_disc)
    names(results) <- c("release_policy", "Bellman", "storage", "releases", "annual_reliability",
                        "time_based_reliability", "volumetric_reliability",
                        "resilience", "vulnerability", "flow_disc")
    
    
    
  } else {
    results <- list(R_policy, Bellman, S, R_rec, Q_disc)
    names(results) <- c("release_policy", "Bellman", "storage", "releases", "flow_disc")
  }
  
  return(results)
}