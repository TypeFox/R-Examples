summary.qLearnS1 <-
function (object, ...){

  res <- list (reg1=object$s1Fit);
  class (res) <- "summary.qLearnS1"
  res
}
