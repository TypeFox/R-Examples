summary.learnIQ1var <-
function (object, ...){

  if (!object$homo){
    res <- list (s1Reg=object$s1VarFit, varType=object$homo);
    class (res) <- "summary.learnIQ1var"
    res
  }
  else{
    res <- list (s1Reg=object$stdDev, varType=object$homo);
    class (res) <- "summary.learnIQ1var"
    res
  }
}
