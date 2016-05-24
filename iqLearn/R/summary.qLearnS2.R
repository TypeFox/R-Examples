summary.qLearnS2 <-
function (object, ...){

  res <- list (reg2=object$s2Fit);
  class (res) <- "summary.qLearnS2"
  res
}
