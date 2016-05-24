
#########################################################
# Count the number of transisitions in a state sequence #
######################################################### 

state.trans <- function(path, no_states){
  # matrix/varbiables for the transition overview, the breakpoint series,
  # and the sojourn times
  m           <- matrix(rep(0, no_states^2), no_states)
  trans       <- rep(0, length(path))
  # calculation of the transition overview and the breakpoints
  # ----------------------------------------------------------
  # loop over the sequence
  for (u in 1:(length(path)-1)){
    # number of transitions
    # analyze type of transition and write it into a matrix
    m[path[u], path[u+1]] <- m[path[u], path[u+1]] + 1
    # location of breakpoints
    # if a breakpoint occurs, write a "1" into the series, else 0
    # the last one is always equal to zero
    if (path[u] != path[u+1]){
      trans[u] <- 1
      }
    }
  # calculation of the sojourn times
  # --------------------------------
  # determine the postitions of the breakpoints
  pos_bp    <- which(trans == 1)
  # calculate all states at breakpoints (before state-switch)
  sj_states <- path[c(pos_bp)]           
  # calculate the sojourn times (exept for the last visited state)
  temp      <- diff(pos_bp)
  sj_times  <- c(pos_bp[1], temp)
  # special treatment necessary for the last visited state
  # last visited state and sojourn in the last visited state
  t_state   <- path[max(pos_bp)+1]
  t_sj_time <- length(path)-max(pos_bp)
  # join the information of the last and the other states, write everything into a matrix
  sj_states <- c(sj_states, t_state)
  sj_times  <- c(sj_times, t_sj_time)
  sojourns  <- rbind(sj_states, sj_times)
  rownames(sojourns) <- c("state", "time")
  # return values
  return(list(transitions = m, breakpoints = trans, sojourns = sojourns))
  }
