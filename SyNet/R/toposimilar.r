toposimilar <- function (dotdata)
{
    if (is.null(class(dotdata)) | class(dotdata) != "dotdata") {
        cat("Argument is not of class 'dotdata' \n")
        return(invisible())
    }
    nsp <- length(dotdata$Label)
    MSTlength <- TopoSim <- matrix(0, nrow = nsp, ncol = nsp)
    diag(MSTlength) <- unlist(lapply(dotdata$MSTsp, function(x) return(x$path)))
    diag(TopoSim) <- 1
    for (i in 1:(nsp-1)) {
      pt1 <- dotdata$occupancy[[i]]
      for (j in (i+1):nsp) {
      		pt2 <- dotdata$occupancy[[j]]
          if(all(pt1%in%pt2) & all(pt2%in%pt1)) {TopoSim[i, j] <- TopoSim[j, i] <- 1; next}
          pts <- union(pt1, pt2)
          MSTlength[i,j] <- MSTlength[j,i] <- mst(as.matrix(dotdata$PtPtDist[pts, pts]))$path
          a <- MSTlength[i,i] + MSTlength[j, j] - MSTlength[i, j]
          if (a > 0) {
            b <- ifelse((MSTlength[i, i] - a) < 0, 0, MSTlength[i, i] - a)
            ce <- ifelse((MSTlength[j, j] - a) < 0, 0, MSTlength[j, j] - a)
            if(abs(b - ce) < 1e-7) {TopoSim[i, j] <- TopoSim[j, i] <- a/(a + b)}
            else TopoSim[i, j] <- TopoSim[j, i] <- a/abs(b - ce)*log((a + max(b, ce))/(a + min(b, ce)))
          } else TopoSim[i, j] <- TopoSim[j, i] <- a
       }
    }
    return(TopoSim)
}

