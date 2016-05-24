Khat <-
function(X, r = NULL, ReferenceType = "", NeighborType = ReferenceType, CheckArguments = TRUE) {

  if (CheckArguments) {
    CheckdbmssArguments()
    # Eliminate erroneous configurations
    if ((ReferenceType == "" | NeighborType == "") & (ReferenceType != NeighborType)) {
      stop("Either two or no point type must be specified.")
    }
  }
  
  # K intra
  if (ReferenceType == "" & NeighborType == "") {
    return (Kest(X, r=r, correction="best"))
  }
  # K intra for a single point type
  if (ReferenceType == NeighborType) {
    X.reduced <- X[X$marks$PointType == ReferenceType]
    return (Kest(X.reduced,  r=r, correction="best"))
  }  
  # K inter calls Kcross. The marks must contain the type, with no weight.
  if (ReferenceType != NeighborType) {
    X.cross <- X
    X.cross$marks <- X$marks$PointType
    return (Kcross(X.cross, i=ReferenceType , j=NeighborType, r=r, correction="best"))
  }  
}
