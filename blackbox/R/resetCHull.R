## possible formats are vertices, vectices001, and constraints that returns both the 'Hrep' and 'a' and 'b'
resetCHull <- function(object,
                     formats=c("vertices"), ## by default returns vertices...
                     redundant.mode="double", ## ...after redundant elimination...
                     quickndirty = FALSE,
                     returnRational=FALSE,
                     verboseThreshold=if (.Platform$r_arch=="x64") {1000} else {600}) {
  if (quickndirty) { ## tries to avoid two bottlenecks: (1) unique.matrix; (2) rational arithm
    vertices <- object[unique(as.numeric(convhulln(object, "Pp"))), ]
    resu <- list()
    if ("vertices" %in% formats ) resu <- c(resu, list(vertices=vertices))
    if ("constraints" %in% formats) {
      Hdec <- cbind(0, 1, vertices)
      Hdec <- scdd(Hdec, representation = "V", roworder="maxcutoff")$output # H repr from V repr.
      if (!is.null(colNames <- colnames(object))) colnames(Hdec) <- c("eq", "b", colNames) ## col names checked by isPointInCHull
      resu <- c(resu, list(Hrep=Hdec, a=-Hdec[, -c(1, 2), drop=FALSE], b=Hdec[, 2]))
    }
    return(resu)
  } ## ELSE
  if("OKriglistplus" %in% class(object)) stop.redef("(!) no 'purefn' method for object of class 'Kriglistplus'")
  if("OKriglist" %in% class(object)) { ## if Kriging by blocs
    fullrange <- object[[1]]$x
    for (jj in 2:(length(object$Kriglist))) {fullrange <- rbind(fullrange, object[[jj]]$x)}
    redundantpts <- unique.matrix(fullrange) ## FR->FR a potential faster code is Unique(as.matrix(fullrange)) using geometry::Unique
  } else if("OKrig" %in% class(object)) {
    redundantpts <- unique.matrix(object$x) ## or Unique(as.matrix(object$x))
  } else redundantpts <- as.matrix(unique(object)) ## or Unique(as.matrix(object)) ## possible SLOW !
  onr <- nrow(redundantpts) ## points, but with redundant rows
  colNames <- colnames(redundantpts)
  if ( onr>2L ) { ## else cannot be redundant
    #    if (onr>verboseThreshold) message("       [* Computing convex hull can be very slow!")
    if (nrow(redundantpts)<=ncol(redundantpts)) { ## convhulln crashes!
      vertices <- redundantpts ## stupid case, should never occur
    } else if (ncol(redundantpts)==1L) {
      vertices <- array(c(min(redundantpts), max(redundantpts)), dim=c(2, 1))
    } else {
      trycn <- try(convhulln(redundantpts, "Pp"),silent=TRUE) ## possibly SLOW
      #if (inherits(trycn,"try-error")) trycn <- try(convhulln_allowing_lowerdim(redundantpts))
      if (inherits(trycn,"try-error")) {
        amknots <- as.matrix(redundantpts) ## the sweep syntax does not work on data.frame's 
        colvecs <- t(sweep(redundantpts[-1,,drop=FALSE],2,redundantpts[1,],`-`))
        qrvecs <- qr(colvecs)
        if (qrvecs$rank<ncol(redundantpts)) {
          message.redef("resetCHull: qhull failed because input points are in a linear subspace\n of lower dimension than the number of coordinates.\n The calling function may handle this problem.")
        } else {
          message.redef("resetCHull: qhull failed for an unidentified reason (input points are not\n in a linear subspace of lower dimension than the number of coordinates.")
        }
        return(trycn)
      } else vertices <- redundantpts[unique(as.numeric(trycn)), ] ## removes redundant vertices
    }
    #    if (onr>verboseThreshold) message("       ... done. *]")
    colnames(vertices) <- colNames
    nr <- nrow(vertices)
  } else if ( onr==2L ) {
    vertices <- redundantpts
  } else if ( onr==1L ) {
    vertices <- redundantpts
    message("convex hull requested for a single 'unique()' point.")
  }
  resu <- list()
  if ("vertices" %in% formats ) resu <- c(resu, list(vertices=vertices))
  if ("vertices001" %in% formats) { ## MINUS vertices, for lpcdd ! May become a rarely used option...
    if ( returnRational) {
      resu <- c(resu, list(vertices001=cbind("0", "0", "1", -vertices)))
    } else {
      resu <- c(resu, list(vertices001=cbind(0, 0, 1, -vertices)))
    }
  }
  if ("constraints" %in% formats) {
    ## vertices are already minimal...
    if (onr > verboseThreshold) { ## test bc system.time itself wastes a comparatively large among of time when nrow(vertices) is ~ 3
      toto <- system.time(Hdec <- Hreprrational(vertices)) ## return value always rational
      if (toto["user.self"]>100) message.redef(paste(" Convex hull (V->H) computation of ",
                                                     paste(dim(vertices),collapse="x"),"vertices\n took ",
                                                     signif(toto["user.self"], 2), " s.",sep=""))
    } else Hdec <- Hreprrational(vertices)
    if ( returnRational ) {
      resu <- c(resu, list(Hrep=Hdec, a=qneg(Hdec[, -c(1, 2), drop=FALSE]), b=Hdec[, 2])) ## redundant info... but Hrep directly useful for subchull...
    } else resu <- c(resu, list(Hrep=q2d(Hdec), a=q2d(qneg(Hdec[, -c(1, 2), drop=FALSE])), b=q2d(Hdec[, 2])))
  }
  return(resu)
}
