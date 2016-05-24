"autoThreshold" <-
function(d.m, est=0.5) {
  est.old <- 0
  while (est.old != est) {
    est.old <- est
    t1 <- mean(d.m[d.m < est], na.rm=TRUE)
    t2 <- mean(d.m[d.m > est], na.rm=TRUE)
    est <- mean(c(t1, t2), na.rm=TRUE)
  }
  return(c(t1, mean(c(t1, est)), est, mean(c(t2,est)), t2))
}
