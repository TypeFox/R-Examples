summary.learnIQ2 <-
function (object, ...){

  res <- list (s2Reg=object$s2Fit);
  class (res) <- "summary.learnIQ2"
  res
}
