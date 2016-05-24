activate.interventions <-
function(P, star, active) {
  if (P$recursive) {
    for (i in 1:length(P$children)) {
      P$children[[i]] <- activate.interventions(P$children[[i]], P$star & P$children[[i]]$star, union(P$do, active))
    }
  } else {
    if (!star) P$star <- FALSE
    P$cond <- setdiff(P$cond, active)
    P$do <- active
    return(P)
  }
  return(P)
}
