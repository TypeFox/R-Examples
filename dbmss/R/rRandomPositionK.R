rRandomPositionK <-
function(X, CheckArguments = TRUE) {
  
  if (CheckArguments)
    CheckdbmssArguments()
  
  # Draw in a binomial process
  RandomizedX <- runifpoint(X$n, win=X$window)
  # Apply original marks to new points
  marks(RandomizedX) <- data.frame(PointWeight=X$marks$PointWeight, PointType=X$marks$PointType)
  class(RandomizedX) <- c("wmppp", "ppp")
  return (RandomizedX)
}
