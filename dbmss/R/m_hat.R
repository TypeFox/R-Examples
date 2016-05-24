mhat <-
function(X, r = NULL, ReferenceType, NeighborType = ReferenceType, CaseControl = FALSE, 
         Original = TRUE, Approximate = ifelse(X$n < 10000, 0, 1), Adjust = 1, MaxRange = "ThirdW", 
         CheckArguments = TRUE) {

  # Eliminate erroneous configurations
  if (CheckArguments) {
    CheckdbmssArguments()
    if (CaseControl & (ReferenceType==NeighborType)) {
      warning("Cases and controls are identical.")
      return(rep(1,length(r)))
    }
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

  # Roughly estimate max distances
  if(is.null(r)) {
    # Default interval for R: between the min distance
    rmin <- 0
    # and the diameter of the window /2 or /4. DO2005 is ignored at this stage.
    rmax <- switch(MaxRange,
                   HalfW = diameter(X$win)/2,
                   ThirdW =  diameter(X$win)/3,
                   QuarterW = diameter(X$win)/4)
    if(is.null(rmax)) rmax <- diameter(X$win)/3
  } else {
    rmin <- 0
    rmax <- max(r)
  }
  
  if (Approximate) {
    # Round distances to save memory
    # Prepare steps so that 1024*Approximate steps are between 0 and rmax. Pairs further than 2*rmax apart are dropped.
    rseq <- seq(from = rmin, to = rmax*2, length.out = 2048*Approximate)
    Nr <- length(rseq)
    # Prepare a matrix (serial version only), one line for each point, one column for each distance
    # Store weights of neighbors of interest in first Nr columns, all points from Nr+1 to 2*Nr
    # Nbd <- matrix(0.0, nrow=X$n, ncol=2*Nr)
    
    # Call C routine to fill Nbd
    if (CaseControl) {
      Nbd <- parallelCountNbdCC(rseq, X$x, X$y, X$marks$PointWeight, IsReferenceType, IsNeighborType)
      # Serial version (returns nothing but modifies Nbd)
      #CountNbdCC(r, X$x, X$y, X$marks$PointWeight, Nbd, IsReferenceType, IsNeighborType)    
    } else {
      Nbd <- parallelCountNbd(rseq, X$x, X$y, X$marks$PointWeight, IsReferenceType, IsNeighborType)
      # Serial version (returns nothing but modifies Nbd)
      # CountNbd(r, X$x, X$y, X$marks$PointWeight, Nbd, IsReferenceType, IsNeighborType)
    }
    
    # Keep the lines of the matrix corresponding to reference points (cases).
    # Other lines are useless and have not been filled by the loops
    Nbd <- Nbd[IsReferenceType, ]
    
    # Adjust distances: values are the centers of intervals
    rseq <- c(0, (rseq[2:Nr]+rseq[1:Nr-1])/2)
    
    # Estimate the bandwith according to adjust if requested.
    if (Original) {
      h <- stats::bw.nrd0(rseq)*Adjust
    } else {
      h <- stats::bw.SJ(rseq)*Adjust
    }

    # Calculate densities of neighbors (with unnormalized weights so suppress warnings)
    Djc <- t(apply(Nbd[, 1:Nr], 1, function(x) suppressWarnings(stats::density(rseq, bw=h, weights=x, cut=0, from=rmin, to=rmax, na.rm=TRUE))$y))
    Dj <- t(apply(Nbd[, (Nr+1):(2*Nr)], 1, function(x) suppressWarnings(stats::density(rseq, bw=h, weights=x, cut=0, from=rmin, to=rmax, na.rm=TRUE))$y))
    # Get the x values of the density estimation: estimate one vector
    x <- stats::density(rseq, bw=h, cut=0, from=rmin, to=rmax, na.rm=TRUE)$x
    
    
  } else {
    # Classical estimation
    
    # Call C routine to fill Nbd. Send x, y, and the list of reference points (indexed from 0 in C instead of 1 in R)
    Nbd <- parallelCountNbdm(X$x, X$y, which(IsReferenceType)-1)
    # Negative values are actually NA
    Nbd[Nbd < 0] <- NA
    
    # Choose the bandwith based on all distance pairs between reference and neighbor points
    # Prepare the data
    RefDistances <- Nbd[, IsNeighborType]
    if (ReferenceType==NeighborType) {
      # Make a square matrix and keep the upper half as a vector
      RefDistances <- RefDistances[upper.tri(RefDistances)]
    }
    
    if (Original) {
      h <- stats::bw.nrd0(RefDistances) * Adjust
    } else {
      h <- stats::bw.SJ(RefDistances) * Adjust
    }

    if(is.null(r)) {
      # Min distance obtained from the data rather than 0
      rmin <- min(Nbd, na.rm=TRUE)
      # Max distance may be obtained from the data rather than from the window
      if (MaxRange == "DO2005") rmax <- stats::median(Nbd, na.rm = TRUE)
    }

    # Calculate densities of neighbors (with unnormalized weights so suppress warnings)
    Djc <- t(apply(Nbd[, IsNeighborType], 1, function(x) suppressWarnings(stats::density(x, bw=h, weights=X$marks$PointWeight[IsNeighborType][!is.na(x)], cut=0, from=rmin, to=rmax, na.rm=TRUE))$y))
    Dj <- t(apply(Nbd, 1, function(x) suppressWarnings(stats::density(x, bw=h, weights=X$marks$PointWeight[!is.na(x)], cut=0, from=rmin, to=rmax, na.rm=TRUE))$y))
    # Get the x values of the density estimation: estimate one vector
    x <- stats::density(Nbd[1, IsNeighborType], bw=h, cut=0, from=rmin, to=rmax, na.rm=TRUE)$x
  }
  
  
  # Calculate the local ratio (at distance r)
  LocalRatio <- Djc/Dj
  # Divide it by the global ratio. Ignore points with no neighbor at all.
  mvalues <- colSums(LocalRatio)/sum(GlobalRatio)
  
  # Interpolate if necessary
  if (is.null(r)) {
    r <- x
  } else {
    mvalues <- stats::approx(x, mvalues, xout=r)$y
  }
  # Put the results into an fv object
  mEstimate <- data.frame(r, rep(1, length(r)), mvalues)
  colnames(mEstimate) <- c("r", "theo", "m")
  
  # Return the values of M(r)
  return (fv(mEstimate, argu="r", ylab=quote(m(r)), valu="m", fmla= ". ~ r", alim=c(0, max(r)), labl=c("r", "%s[ind](r)", paste("hat(%s)(r)", sep="")), desc=c("distance argument r", "theoretical independent m(r)", "Estimated m(r)"), unitname=X$window$unit, fname="m")) 
}
