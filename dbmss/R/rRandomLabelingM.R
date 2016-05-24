rRandomLabelingM <-
function(X, CheckArguments = TRUE) {
  
  if (CheckArguments)
    CheckdbmssArguments()
  
  # Randomize marks
  RandomizedX <- rlabel(X)
  # Restore weights
  RandomizedX$marks$PointWeight <- X$marks$PointWeight
  
  class(RandomizedX) <- c("wmppp", "ppp")
  return (RandomizedX)
}
