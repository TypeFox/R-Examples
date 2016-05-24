rPopulationIndependenceM <-
function(X, ReferenceType, CheckArguments = TRUE) {
  
  if (CheckArguments)
    CheckdbmssArguments()
  
  ReferencePP <- X[X$marks$PointType==ReferenceType]
  OtherPointsPP <- X[X$marks$PointType!=ReferenceType]
  RandomizedX <- superimpose(ReferencePP, rlabel(OtherPointsPP))
  class(RandomizedX) <- c("wmppp", "ppp")
  return (RandomizedX)
}
