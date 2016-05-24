Mhat <-
function(X, r = NULL, ReferenceType, NeighborType = ReferenceType, CaseControl = FALSE, CheckArguments = TRUE) {
  # Eliminate erroneous configurations
  if (CheckArguments) {
    CheckdbmssArguments()
    if (CaseControl & (ReferenceType==NeighborType)) {
      warning("Cases and controls are identical.")
      return(rep(1,length(r)))
    }
  }
  
  # Default r values: 64 values up to half the max distance
  if (is.null(r)) {
    rMax <- diameter(X$window)
    r <- rMax*c(0, 1:20, seq(22, 40, 2), seq(45, 100,5), seq(110, 200, 10), seq(220, 400, 20))/800
  }
  
  # Vectors to recognize point types
  IsReferenceType <- X$marks$PointType==ReferenceType
  IsNeighborType <- X$marks$PointType==NeighborType
  
  # Global ratio
  if (ReferenceType==NeighborType | CaseControl) {
    WrMinusReferencePoint <- sum(X$marks$PointWeight[IsReferenceType])-X$marks$PointWeight
    Wn <- WrMinusReferencePoint[IsReferenceType]
  } else {
    Wn <- sum(X$marks$PointWeight[IsNeighborType])
  }
  if (CaseControl) {
    Wa <- sum(X$marks$PointWeight[IsNeighborType]) 
  } else {
    WaMinusReferencePoint <- sum(X$marks$PointWeight)-X$marks$PointWeight
    Wa <- WaMinusReferencePoint[IsReferenceType]
  }
  GlobalRatio <- Wn/Wa
  
  Nr <- length(r)
  # Neighborhoods (i.e. all neighbors of a point less than a distance apart)
  # Prepare matrix (serial version only), one line for each point, one column for each distance
  # Store weights of neighbors of interest in first Nr columns, all points from Nr+1 to 2*Nr
  # Nbd <- matrix(0.0, nrow=X$n, ncol=2*Nr)
  
  # Call C routine to fill Nbd
  if (CaseControl) {
    Nbd <- parallelCountNbdCC(r, X$x, X$y, X$marks$PointWeight, IsReferenceType, IsNeighborType)
    # Serial version (returns nothing but modifies Nbd)
    #CountNbdCC(r, X$x, X$y, X$marks$PointWeight, Nbd, IsReferenceType, IsNeighborType)    
  } else {
    Nbd <- parallelCountNbd(r, X$x, X$y, X$marks$PointWeight, IsReferenceType, IsNeighborType)
    # Serial version (returns nothing but modifies Nbd)
    # CountNbd(r, X$x, X$y, X$marks$PointWeight, Nbd, IsReferenceType, IsNeighborType)
  }
  
  # Keep the lines of the matrix corresponding to reference points (cases).
  # Other lines are useless and have not been filled by the loops
  NbdInt <- Nbd[IsReferenceType, 1:Nr]
  NbdAll <- Nbd[IsReferenceType, (Nr+1):(2*Nr)]
  # Cumulate weights up to each distance
  NbdInt <- t(apply(NbdInt, 1, cumsum))
  NbdAll <- t(apply(NbdAll, 1, cumsum))
  
  # Calulate the ratio of points of interest around each point
  LocalRatio <- NbdInt/NbdAll
  # Divide it by the global ratio. Ignore points with no neighbor at all.
  Mvalues <- apply(LocalRatio, 2, function(x) sum(x[is.finite(x)])/sum(GlobalRatio[is.finite(x)]))
  # Put the results into an fv object
  MEstimate <- data.frame(r, rep(1, length(r)), Mvalues)
  colnames(MEstimate) <- c("r", "theo", "M")
  
  # Return the values of M(r)
  return (fv(MEstimate, argu="r", ylab=quote(M(r)), valu="M", fmla= ". ~ r", alim=c(0, max(r)), labl=c("r", "%s[ind](r)", paste("hat(%s)(r)", sep="")), desc=c("distance argument r", "theoretical independent M(r)", "Estimated M(r)"), unitname=X$window$unit, fname="M"))
  
}
