deconstruct <- function(P, P.res) {
  if (!P$recursive) {
    init <- length(P.res$children) 
    if (length(P$sumset) == 0) {
      P.temp <- probability(var = P$var[1], cond = union(P$cond, P$var[-1]))
      P.res$children[[init + 1]] <- P.temp
      if (length(P$var) > 1) {
        for (i in 2:length(P$var)) {
          P.temp <- probability(var = P$var[i], cond = union(P$cond, P$var[-(1:i)]))
          P.res$children[[init + i]] <- P.temp
        }
      }
    } else {
      P.temp <- probability(recursive = TRUE, children = list(), sumset = P$sumset) 
      P.temp$children[[1]] <- probability(var = P$var[1], cond = union(P$cond, P$var[-1]), sumset = P$sumset)
      if (length(P$var) > 1) {
        for (i in 2:length(P$var)) {
          P.temp$children[[i]] <- probability(var = P$var[i], cond = union(P$cond, P$var[-(1:i)]))
        }
      }
      P.res$children[[init + 1]] <- P.temp
    }
  } 
  if (P$recursive & length(P$sumset) > 0) {
    P.temp <- probability(recursive = TRUE, children = P$children)
    P.temp <- deconstruct(P.temp, probability(recursive = TRUE, children = list()))
    P.temp$sumset <- P$sumset
    P.res$children[[length(P.res$children)+1]] <- P.temp
  }
  if (P$fraction) {
    P.res$divisor <- deconstruct(P$divisor, P.res$divisor)
  }
  if (P$recursive & length(P$sumset) == 0) {
    for (i in 1:length(P$children)) {
      P.res <- deconstruct(P$children[[i]], P.res)
    }
  }
  return(P.res)
}