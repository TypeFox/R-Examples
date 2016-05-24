valid_amp <- function(x) {
  tres <- t.test(head(x), tail(x), alternative="less")$p.value < 0.01
  sigres <- mean(tail(x)) > mean(x[5:15]) + 3*sd(x[5:15])
  as.logical(tres * sigres)
}
