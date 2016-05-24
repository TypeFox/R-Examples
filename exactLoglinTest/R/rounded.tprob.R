rounded.tprob <- function(y, m, s, df) {
  t1 <- (y + .5 - m) / sqrt(s)
  t2 <- (y - .5 - m) / sqrt(s)
  sum(log(pt(t1, df) - pt(t2, df)))
}
