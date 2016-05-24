Kdhat <-
function(X, r = NULL, ReferenceType, NeighborType = ReferenceType, Weighted = FALSE,
         Original = TRUE, Approximate = ifelse(X$n < 10000, 0, 1), Adjust = 1,
         MaxRange = "ThirdW", CheckArguments = TRUE) {
  
  if (CheckArguments) {
    CheckdbmssArguments()
  }
  
  # Select the bandwith: original choice by Duranton and Overman or optimized one.
  if (Original) {
    bw <- "nrd0"
  } else {
    bw <- "sj"
  }
  
  # Vectors to recognize point types
  if (ReferenceType == "") {
    # All points (reference value as the center of the confidence interval)
    IsReferenceType <- IsNeighborType <- rep(TRUE, X$n)
    Y <- X
  } else {
    # Current use
    IsReferenceType <- X$marks$PointType==ReferenceType
    IsNeighborType <- X$marks$PointType==NeighborType    
    # Eliminate useless points
    Y <- X[IsReferenceType | IsNeighborType]
    # Update for Y
    IsReferenceType <- Y$marks$PointType==ReferenceType
    IsNeighborType <- Y$marks$PointType==NeighborType
  }

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
    # Prepare steps so that 1024*Approximate steps are between 0 and rmax.
    # Pairs further than 2*rmax apart will be stored in an extra element.
    rseq <- seq(from = rmin, to = rmax*2, length.out = 2048*Approximate)
    # Number of distances
    Nr <- length(rseq)
    # Prepare a matrix, single line, one value for each distance + 1 extra for pairs far away.
    NeighborWeight <- matrix(0.0, nrow=1, ncol=Nr+1)
    # Weights
    if (Weighted) {
      Weight <- Y$marks$PointWeight
    } else {
      Weight <- rep(1, Y$n)
    }
    
    # Call C routine to fill NeighborWeights
    CountNbdKd(rseq, Y$x, Y$y, Weight, NeighborWeight, IsReferenceType, IsNeighborType)
    
    # Adjust distances: values are the centers of intervals
    rseq <- c(0, (rseq[2:Nr]+rseq[1:Nr-1])/2)

    # Estimate the density. Change the bandwith according to adjust if requested.
    if (Adjust != 1) {
      if (Original) {
        bw <- stats::bw.nrd0(rseq) * Adjust
      } else {
        bw <- stats::bw.SJ(rseq) * Adjust
      }
    }
    
    # The last element of the vector NeighborWeight contains the weight of pairs farther than 2 rmax
    FarWeight <- NeighborWeight[Nr+1]
    # Keep the other elements
    NeighborWeight <- NeighborWeight[-(Nr+1)]
    
    # Estimate density taking into account far pairs. Suppress warnings because density does not sum to 1.
    Density <- suppressWarnings(stats::density(rseq, Weight=NeighborWeight/(sum(NeighborWeight)+FarWeight), cut=0, to=rmax, bw=bw))
    
  } else {
    # Classical estimation
    # Prepare a vector for distances between all point pairs.
    if (ReferenceType == NeighborType) {
      # Univariate Kd: n(n-1)/2 pairs
      NbDist <- sum(IsReferenceType)*(sum(IsReferenceType)-1)/2
    } else {
      # Bivariate Kd: n1*n2 pairs
      NbDist <- sum(IsReferenceType)*sum(IsNeighborType)
    }
    Dist <- vector(mode="double", length=NbDist)
    
    # Prepare a vector for weights if Weighted. Else, set a single value.
    if (Weighted) {
      Weight <- vector(mode="double", length=NbDist)
    } else {
      Weight <- 1
    }
    
    # C++ routine to fill distances and weights
    DistKd(Y$x, Y$y, Y$marks$PointWeight, Weight, Dist, IsReferenceType, IsNeighborType)
    
    if(is.null(r)) {
      # Min distance obtained from the data rather than 0
      rmin <- min(Dist)
      # Max distance may be obtained from the data rather than from the window
      if (MaxRange == "DO2005") rmax <- stats::median(Dist)
    }
    
    # Estimate the density. Change the bandwith according to adjust if requested.
    if (Adjust != 1) {
      if (Original) {
        bw <- stats::bw.nrd0(Dist) * Adjust
      } else {
        bw <- stats::bw.SJ(Dist) * Adjust
      }
    }
    # Increase the number of estimation points if necessary. 
    # If rmax << Dist, too few estimation points may be below rmax
    # Try to have at least 128 of them. Limit the total number of points to 4096 (i.e. rmax must be >1/32 max(Dist))
    nDensity <- max(min(512, max(Dist)/rmax*128), 4096)
    # Estimate density
    if (Weighted) {
      Density <- stats::density(Dist, Weight=Weight/sum(Weight), cut=0, from=rmin, to=rmax, bw=bw, n=nDensity)
    } else {
      Density <- stats::density(Dist, cut=0, from=rmin, to=rmax, bw=bw, n=nDensity)
    }
  }
  

  if(is.null(r)) {
    # Return estimated values
    r <- Density$x
    Kd <- Density$y
  } else {
    # Interpolate results at the chosen R
    Kd <- stats::approx(Density$x, Density$y, xout=r)$y    
  }
  KdEstimate <- data.frame(r, Kd)
  colnames(KdEstimate) <- c("r", "Kd")
  
  # Return the values of Kd(r)
  return (fv(KdEstimate, argu="r", ylab=quote(Kd(r)), valu="Kd", fmla= ". ~ r", alim=c(0, max(r)), labl=c("r", paste("hat(%s)(r)", sep="")), desc=c("distance argument r", "Estimated Kd(r)"), unitname=X$window$unit, fname="Kd"))
}
