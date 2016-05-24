LogRank <- function (y, tr, tl, event, block, ...) {
  ## Log-rank statistic
  if (is.null(block)) {
    lr <- coin::surv_test(Surv(y, event) ~ factor(tr==tl))
  } else {
    lr <- coin::surv_test(Surv(y, event) ~ factor(tr==tl) | factor(block))
  }
  return(coin::statistic(lr, "test"))
}
