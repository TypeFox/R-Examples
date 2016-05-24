#' @import stats

Bi.test <- function(tot, act)
{
  lower <- 0
  upper <- 1
  threshold <- 1
  if(act == 0)
    test <- 0
  else {
    if(act == tot)
      test <- 0.05^(1/act)
    else {
      while(threshold > 0.0001) {
        test <- (lower + upper)/2
        prop <- pbinom((act - 1),tot,test)
        if(prop < 0.95)
          upper <- test
        else lower <- test
        threshold <- abs(prop - 0.95)
      }
    }
  }
  return(test)
}
