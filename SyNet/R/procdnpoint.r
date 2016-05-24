procdnpoint <- function (dnpoint, tolerance = 1e-03) {
    if (is.null(class(dnpoint)) | class(dnpoint) != "dnpoint") {
        cat("Argument is not of class 'dnpoint' \n")
        return(invisible())
    }
    stopifnot(is.vector(tolerance) & is.numeric(tolerance)) 
    if(length(tolerance) != dnpoint$Numpoints) tolerance <- rep(tolerance, length.out = dnpoint$Numpoints) 
    o <- order(dnpoint$Points[, 2], dnpoint$Points[, 3], dnpoint$Points[, 1])
    dnpoint$Points <- dnpoint$Points[o, ]
    tolerance <- tolerance[o]
    ident <- apply(dnpoint$Points, 2, function(x) return(abs(diff(x)) > 1e-10))
    difpt <- c(TRUE, apply(ident[,2:3], 1, any)) #TRUE for unique points
    coords <- dnpoint$Points[difpt,2:3]
    #Matrix of distances between points
    if(dnpoint$Type == "geographical") {
      latrad <- coords[,2]*pi/180
      lonrad <- coords[,1]*pi/180
      R <- 6378.388
      coslat1 <- cos(latrad)
      sinlat1 <- sin(latrad)
      coslon1 <- cos(lonrad)
      sinlon1 <- sin(lonrad)
      pp <- cbind(coslat1 * coslon1, coslat1 * sinlon1, sinlat1) %*%
              t(cbind(coslat1 * coslon1, coslat1 * sinlon1, sinlat1))
      PtPtDist <- R * acos(ifelse(abs(pp) > 1, 1 * sign(pp), pp))
    }
    if(dnpoint$Type == "cartesian") PtPtDist <- as.matrix(dist(coords))
    diag(PtPtDist) <- 0
    nsp <- length(dnpoint$Label)
    npt <- sum(difpt)
    tol <- abs(tolerance[difpt])
    dntable <- matrix(FALSE, nrow = nsp, ncol = npt)
    rownames(dntable) <- dnpoint$Label
    j <- 0
    #Boolean table of species x actual records 
    for(i in 1:dnpoint$Numpoints) {
      if(difpt[i]) j <- j + 1
      dntable[dnpoint$Points[i,1], j] <- TRUE
    }
    PtSpDist <- c() # Hausdorff distance between single points and species set points 
    #List of species records (plus aggregated points by tolerance radius) and respective MST information
    occupancy <- MSTsp <- vector("list", nsp)
    for (i in 1:nsp) {
      ref <- apply(PtPtDist[, which(dntable[i,]), drop = FALSE], 1, min)
      pts <- which(ref <= tol) #Should be actual points in addition to those found at the 
                               #tolerated immediacy
      PtSpDist <- cbind(PtSpDist, ref)
      occupancy[[i]] <- pts 
      MSTsp[[i]] <- mst(as.matrix(PtPtDist[pts, pts, drop = FALSE]))
      MSTsp[[i]][[3]] <- pts[MSTsp[[i]][[3]]]
      MSTsp[[i]][[4]] <- pts[MSTsp[[i]][[4]]]
    }
    names(occupancy) <- dnpoint$Label
    rslts <- list(Call = match.call(), Label = dnpoint$Label, dntable = dntable, occupancy = occupancy, 
             coords = coords, PtPtDist = PtPtDist, PtSpDist = PtSpDist, MSTsp = MSTsp)
    class(rslts) <- "dotdata"
    rslts
}

