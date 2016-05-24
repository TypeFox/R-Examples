rPopulationIndependenceK <-
function (X, ReferenceType, NeighborType, CheckArguments = TRUE) {

  if (CheckArguments)
    CheckdbmssArguments()

  # Eliminate useless points
  X.reduced <- X[X$marks$PointType==ReferenceType | X$marks$PointType==NeighborType]
  RandomizedX <- X.reduced
  # Reduce the factor levels to two (factor eliminates the levels with no points)
  Marks <- factor(X.reduced$marks$PointType)
  # The new point pattern has classical spatstat marks
  RandomizedX <- RandomizedX %mark% Marks
  # Split reference and neighbor points
  X.split <- split(RandomizedX)
  # Randomly shift the neighbors
  rshift(X.split, which=NeighborType) -> RandomizedX.split
  # Reunify the split point pattern
  RandomizedX.split ->  split(RandomizedX)
  # Reorganize the marks (add weight)
  PointWeight <- rep(1, RandomizedX$n)
  PointType   <- marks(RandomizedX)
  marks(RandomizedX) <- data.frame(PointWeight, PointType)  
  class(RandomizedX) <- c("wmppp", "ppp")
  return (RandomizedX)
}
