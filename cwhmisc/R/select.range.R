"select.range" <- function (data, groupvec, min, max) {
  min.cond <- (groupvec >= min)
  max.cond <- (groupvec < max)
  cond <- min.cond & max.cond
  return(data[cond])
}
