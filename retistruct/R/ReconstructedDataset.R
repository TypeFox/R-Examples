##' This function infers the coordinates of datapoints \code{Ds }and
##' landmarks \code{Ss} in  spherical coordinates.
##'
##' @title Constructor for RecontructedDataset object
##' @param r Object that of clases \code{reconstructedOutline} and
##' \code{dataset}.
##' @param report Function used to report progress.
##' @return \code{\link{ReconstructedDataset}} object containing the input
##' information and the following modified and extra information:
##' \item{\code{Dsb}}{Datapoints in barycentric coordinates}
##' \item{\code{Dsc}}{Datapoints on reconstructed sphere in cartesian coordinates}
##' \item{\code{Dss}}{Datapoints on reconstructed sphere in spherical coordinates}
##' \item{\code{Ssb}}{Landmarks in barycentric coordinates}
##' \item{\code{Ssc}}{Landmarks on reconstructed sphere in cartesian coordinates}
##' \item{\code{Sss}}{Landmarks on reconstructed sphere in spherical coordinates}
##' @author David Sterratt
##' @import geometry
##' @export
ReconstructedDataset <- function(r, report=message) {
  report("Inferring coordinates of datapoints")
  Dsb <- list() # Datapoints in barycentric coordinates
  Dsc <- list() # Datapoints on reconstructed sphere in cartesian coordinates
  Dss <- list() # Datapoints on reconstructed sphere in spherical coordinates
  if (!is.null(r$Ds) & (length(r$Ds) > 0)) {
    for (name in names(r$Ds)) {
      Dsb[[name]] <- tsearchn(r$P, r$T, r$Ds[[name]])
      oo <- is.na(Dsb[[name]]$idx)     # Points outwith outline
      if (any(oo)) {
        warning(paste(sum(oo), name, "datapoints outwith the outline will be ignored."))
      }
      Dsb[[name]]$p   <- Dsb[[name]]$p[!oo,,drop=FALSE]
      Dsb[[name]]$idx <- Dsb[[name]]$idx[!oo]
      Dsc[[name]] <- bary.to.sphere.cart(r$phi, r$lambda, r$R, r$Tt, Dsb[[name]])
      Dss[[name]] <- sphere.cart.to.sphere.spherical(Dsc[[name]], r$R)
    }
  }
  
  report("Inferring coordinates of grouped datapoints")
  Gsb <- list() # Datapoints in barycentric coordinates
  Gsc <- list() # Datapoints on reconstructed sphere in cartesian coordinates
  Gss <- list() # Datapoints on reconstructed sphere in spherical coordinates
  if (!is.null(r$Gs) & (length(r$Gs) > 0)) {
    for (name in names(r$Gs)) {
      Gsb[[name]] <- tsearchn(r$P, r$T, r$Gs[[name]][,1:2,drop=FALSE])
      oo <- is.na(Gsb[[name]]$idx)     # Points outwith outline
      if (any(oo)) {
        warning(paste(sum(oo), name, "datapoints outwith the outline will be ignored."))
      }
      Gsb[[name]]$p   <- Gsb[[name]]$p[!oo,,drop=FALSE]
      Gsb[[name]]$idx <- Gsb[[name]]$idx[!oo]
      Gsc[[name]] <- bary.to.sphere.cart(r$phi, r$lambda, r$R, r$Tt, Gsb[[name]])
      Gss[[name]] <- cbind(sphere.cart.to.sphere.spherical(Gsc[[name]], r$R),
                           r$Gs[[name]][!oo,3])
    }
  }
  
  report("Inferring coordinates of landmarks")
  Ssb <- list() # Landmarks in barycentric coordinates
  Ssc <- list() # Landmarks on reconstructed sphere in cartesian coordinates
  Sss <- list() # Landmarks on reconstructed sphere in spherical coordinates
  if (!is.null(r$Ss) & (length(r$Ss) > 0)) {
    for (i in 1:length(r$Ss)) {
      Ssb[[i]] <- with(r, tsearchn(P, T, r$Ss[[i]]))
      Ssc[[i]] <- bary.to.sphere.cart(r$phi, r$lambda, r$R, r$Tt, Ssb[[i]])
      Sss[[i]] <- sphere.cart.to.sphere.spherical(Ssc[[i]], r$R)
    }
    names(Ssb) <- names(r$Ss)
    names(Ssc) <- names(r$Ss)
    names(Sss) <- names(r$Ss)
  }

  d <- merge(list(Dsb=Dsb, Dsc=Dsc, Dss=Dss,
                  Gsb=Gsb, Gsc=Gsc, Gss=Gss,
                  Ssb=Ssb, Ssc=Ssc, Sss=Sss), r)
  ## Trigger recomputation of Kernel density estimates
  d$KDE <- NULL
  d$KR <- NULL

  class(d) <- addClass("reconstructedDataset", r)
  return(d)
}

##' @title Get IDs of groups of data within a ReconstructedDataset
##' @param r \code{\link{ReconstructedDataset}} object
##' @return Array of IDs
##' @author David Sterratt
##' @method getIDs reconstructedDataset
##' @export
getIDs.reconstructedDataset <- function(r) {
  return(names(r$Dss))
}

##' Get spherical coordinates of datapoints.
##'
##' @title Get transformed spherical coordinates of datapoints
##' @param r \code{ReconstructedDataset} object.
##' @return \code{Dss}
##' @method getDss reconstructedDataset
##' @author David Sterratt
##' @export
getDss.reconstructedDataset <- function(r) {
  return(r$Dss)
}

##' @title Get grouped variable with locations in spherical coordinates.
##' @param r \code{ReconstructedDataset} object.
##' @return \code{Gss}
##' @method getGss reconstructedDataset
##' @author David Sterratt
##' @export
getGss.reconstructedDataset <- function(r) {
  return(r$Gss)
}

##' Get Karcher mean of datapoints in spherical coordinates.
##'
##' @title Karcher mean of datapoints in spherical coordinates
##' @param r \code{\link{ReconstructedDataset}} or \code{\link{RetinalReconstructedDataset}} object.
##' @return \code{Dss.mean}
##' @method getDssMean reconstructedDataset
##' @author David Sterratt
##' @export
getDssMean.reconstructedDataset <- function(r) {
  Dss.mean <- list()
  if (length(r$Dss)) {
    for (i in 1:length(r$Dss)) {
      km <- karcher.mean.sphere(r$Dss[[i]], na.rm=TRUE)
      Dss.mean[[i]] <- cbind(phi=km["phi"], lambda=km["lambda"])
    }
  }
  names(Dss.mean) <- names(r$Dss)
  return(Dss.mean)
}

##' @title Get area of convex hull around data points on sphere
##' @param r code{\link{ReconstructedDataset}} or \code{\link{RetinalReconstructedDataset}} object.
##' @return Area in degress squared
##' @author David Sterratt
##' @export
getDssHullarea <- function(r) {
  Dss.hullarea <- list()
  if (length(r$Dss)) {
    for (i in 1:length(r$Dss)) {
      if (nrow(r$Dss[[i]]) >= 3) {
        Dsp <- sphere.spherical.to.polar.cart(r$Dss[[i]], pa=TRUE)
        Dspt <- suppressMessages(delaunayn(Dsp))
        Dss.hullarea[[i]] <- sum(sphere.tri.area(r$Dss[[i]], Dspt))*(180/pi)^2
      } else {
        Dss.hullarea[[i]] <- NA
      }
    }
  }
  names(Dss.hullarea) <- names(r$Dss)
  return(Dss.hullarea)
}

##' Get spherical coordinates of landmarks.
##'
##' @title Get transformed spherical coordinates of landmarks.
##' @param r \code{ReconstructedDataset} object.
##' @return \code{Sss}
##' @method getSss reconstructedDataset
##' @author David Sterratt
##' @export
getSss.reconstructedDataset <- function(r) {
  return(r$Sss)
}

##' @export
getSssMean.reconstructedDataset <- function(r) {
  Sss.mean <- list()
  if (length(r$Sss)) {
    for (i in 1:length(r$Sss)) {
      km <- karcher.mean.sphere(r$Sss[[i]], na.rm=TRUE)
      Sss.mean[[i]] <- cbind(phi=km["phi"], lambda=km["lambda"])
    }
  }
  names(Sss.mean) <- names(r$Sss)
  return(Sss.mean)
}

##' @title Get kernel density estimate of data points
##' @param r \code{\link{ReconstructedDataset}} object
##' @return See \code{\link{compute.kernel.estimate}}
##' @author David Sterratt
##' @export
getKDE <- function(r) {
  if (is.null(r$KDE)) {
    return(compute.kernel.estimate(getDss(r), r$phi0, kde.fhat, kde.compute.concentration))
  }
  return(r$KDE)
}

##' Compute a kernel estimate over a grid and do a contour analsysis
##' of this estimate. The contour heights the determined by finding
##' heights that exclude a certain fraction of the probability. For
##' example, the 95% contour is excludes 95% of the probability mass,
##' and it should enclose about 5% of the points. The contour levels
##' are specified by  the \code{contour.levels} option; by default
##' they are \code{c(5, 25, 50, 75, 95)}.
##' 
##' @title Kernel estimate over grid
##' @param Dss List of datasets. The first two columns of each datasets
##' are coordinates of points on the sphere in spherical polar
##' (lattitude, \code{phi}, and longitude, \code{lambda})
##' coordinates. In the case kernel smoothing, there is a third column
##' of values of dependent variables at those points.
##' @param phi0 Rim angle in radians
##' @param fhat Function such as \code{\link{kde.fhat}} or
##' \code{\link{kr.yhat}} to compute the density given data and a
##' value of the concentration parameter \code{kappa} of the Fisher
##' density.
##' @param compute.conc Function to return the optimal value of the
##' concentration parameter kappa given the data.
##' @return A list containing
##' \item{\code{kappa}}{The concentration parameter}
##' \item{\code{h}}{A pseudo-bandwidth parameter, the inverse of the square root of \code{kappa}. Units of degrees.}
##' \item{\code{flevels}}{Contour levels.}
##' \item{\code{labels}}{Labels of the contours.}
##' \item{\code{g}}{Raw density estimate drawn on non-area-preserving projection. Comprises locations of gridlines in Cartesian coordinates (\code{xs} and \code{ys}), density estimates at these points, \code{f} and location of maximum in Cartesian coordinates (\code{max}).}
##' \item{\code{gpa}}{Raw density estimate drawn on area-preserving projection. Comprises same elements as above.}
##' \item{\code{contour.areas}}{Area of each individual contour. One level may ahave more than one contour; this shows the areas of all such contours.}
##' \item{\code{tot.contour.areas}}{Data frame containing the total area within the contours at each level.}
##' @author David Sterratt
##' @export
compute.kernel.estimate <- function(Dss, phi0, fhat, compute.conc) {
  vols <- getOption("contour.levels")
  res <- 100

  ## Helper function to get kde as locations gs in spherical coordinates
  get.kde <- function(gs, mu, kappa, res, crop=TRUE) {
    ## Make space for the kernel density estimates
    gk <- fhat(gs, mu, kappa)

    ## If we're cropping, were going to eliminate pixels whos centres
    ## lie outwith the outline altogether. If crop isn't set, there's
    ## a wee margin, that ought to allow all contours to lie outwith
    ## the outline.
    if (crop) {
      gk[gs[,"phi"] > phi0] <- NA
    }
    ## else {
    ##  gk[gs[,"phi"] > phi0 + 2/res*(phi0+pi/2)] <- 0
    ## }

    ## Put the estimates back into a matrix. The matrix is filled up
    ## column-wise, so the matrix elements should match the elements of
    ## gxs and gys
    k <- matrix(gk, res, res)
    k[is.na(k)] <- 0
    return(k)
  }
  
  ## Get data points
  KDE <- list()
  if (length(Dss) > 0) {
    ## First create a grid in Cartesian coordinates with
    ## area-preserving coords. The extra margin on the grid is needed
    ## so that the contour lines are continuous round the edge.
    gpa <- create.polar.cart.grid(TRUE, res, min(phi0 + 0.2*(phi0+pi/2), pi/2))
    ## And one without area-preserving coords
    g   <- create.polar.cart.grid(FALSE, res, min(phi0 + 0.2*(phi0+pi/2), pi/2))

    ## Check conversion
    ## gcb <- sphere.spherical.to.polar.cart(gs, pa)
    ## points(rho.to.degrees(gcb, phi0, pa), pch='.')
    
    for (i in names(Dss)) {
      if (nrow(Dss[[i]]) > 2) {
        ## Find the optimal concentration of the kernel density
        ## estimator
        kappa <- compute.conc(Dss[[i]])
        
        ## Now we've found the concentration, let's try to estimate
        ## and display the density at grid points on an azimuthal
        ## equidistant projection (f) and on an aziumuthal equal-area
        ## projection (fpa).
        fpa <- get.kde(gpa$s, Dss[[i]], kappa, res)
        f  <-  get.kde(g$s,   Dss[[i]], kappa, res)
        maxs   <- gpa$s[which.max(fpa),,drop=FALSE]

        ## The above estimates are set to NA outwith the outline. For
        ## the purposes of computing smooth contours, this causes
        ## jagged edges, so we also get uncropped versions there the
        ## density spreads outwith the outline.
        fpau <- get.kde(gpa$s, Dss[[i]], kappa, res, crop=FALSE)
        fu  <-  get.kde(g$s,   Dss[[i]], kappa, res, crop=FALSE)

        ## Determine the value of gk that encloses 0.95 of the
        ## density.  To compute the density, we need to know the
        ## area of each little square, which is why we have used the
        ## are-preserving projection. FIXME: I think this method of
        ## "integration" could be improved.
        vol.contours <- TRUE
        if (vol.contours) {
          f.sort <- sort(as.vector(fpa))
          js <- findInterval(vols/100, cumsum(f.sort)/sum(f.sort))
          flevels <- f.sort[js]
          if (length(unique(flevels)) < length(flevels)) {
            warning("The volume contours method has found duplicated levels - there is probably something wrong with the data.")
            KDE[[i]] <- NULL
            next
          }
        } else {
          flevels <- vols/100*max(fpa)
        }

        ## Store full kde matrices
        KDE[[i]] <- list(kappa=kappa,
                         h=180/pi/sqrt(kappa),
                         flevels=flevels,
                         maxs=maxs,
                         g=  list(xs=g$xs,   ys=g$ys,   f=f  , fu=fu),
                         gpa=list(xs=gpa$xs, ys=gpa$ys, f=fpa, fu=fpau))

        ## Get contours in Cartesian space
        cc <- contourLines(gpa$xs, gpa$ys, fpau, levels=flevels)
        cs <- list()
        ## Must be careful, as there is a core function called labels
        labels <- rep(NA, length(cc))
        contour.areas <- rep(NA, length(cc))
        if (length(cc) > 0) {
          for (j in 1:length(cc)) {
            cs[[j]] <- list()
            contour.areas[j] <- polyarea(cc[[j]]$x, cc[[j]]$y) * (180/pi)^2
            ccj <- cbind(x=cc[[j]]$x, y=cc[[j]]$y)
            ## Convert back to Spherical coordinates
            cs[[j]] <- polar.cart.to.sphere.spherical(ccj, TRUE)
            ## Push any points outwith the outline back into it
            cs[[j]][cs[[j]][,"phi"] > phi0, "phi"] <- phi0
            labels[j] <- vols[which(flevels==cc[[j]]$level)]
          }
        }
        KDE[[i]]$contours <- cs
        KDE[[i]]$labels <- labels
        names(contour.areas) <- c()
        KDE[[i]]$contour.areas <- contour.areas
        KDE[[i]]$tot.contour.areas <- aggregate(contour.areas ~ labels,
                                                data.frame(labels, contour.areas), sum)
      }
    }
  }
  return(KDE)
}

##' @title Get kernel regression estimate of grouped data points
##' @param r \code{\link{ReconstructedDataset}} object
##' @return See \code{\link{compute.kernel.estimate}}
##' @author David Sterratt
##' @export
getKR <- function(r) {
  if (is.null(r$KR)) {
    yhat <- function(r, mu, kappa) {
      kr.yhat(r, mu[,1:2], mu[,3], kappa)
    }
    compute.conc <- function(mu) {
      kr.compute.concentration(mu[,1:2], mu[,3])
    }
    Gss <- getGss(r)
    ## Get rid of datasets with very little data
    for (n in names(Gss)) {
      if (sum(Gss[[n]][,3]) <= 2) {
        Gss[[n]] <- NULL
      }
    }
    return(compute.kernel.estimate(Gss, r$phi0, yhat, compute.conc))
  }
  return(r$KR)
}

##' @title Plot projection of reconstructed dataset
##' @param r \code{\link{ReconstructedDataset}} object
##' @param transform Transform function to apply to spherical coordinates
##' before rotation
##' @param axisdir Direction of axis (North pole) of sphere in external space
##' @param projection Projection in which to display object,
##' e.g. \code{\link{azimuthal.equalarea}} or \code{\link{sinusoidal}}
##' @param proj.centre Location of centre of projection as matrix with
##' column names \code{phi} (elevation) and \code{lambda} (longitude).
##' @param lambdalim Limits of longitude (in degrees) to display
##' @param datapoints If \code{TRUE}, display data points
##' @param datapoint.means If \code{TRUE}, display Karcher mean of data points.
##' @param datapoint.contours If \code{TRUE}, display contours around
##' the data points generated using Kernel Density Estimation.
##' @param grouped If \code{TRUE}, dipslay grouped data.
##' @param grouped.contours If \code{TRUE}, display contours around
##' the grouped data generated using Kernel Regression.
##' @param landmarks If \code{TRUE}, dipslay landmarks.
##' @param ids IDs of groups of data within a dataset, returned using
##' \code{\link{getIDs}}.
##' @param ... Graphical parameters to pass to plotting functions
##' @method projection reconstructedDataset
##' @author David Sterratt
##' @export
projection.reconstructedDataset <-
  function(r, 
           transform=identity.transform,
           axisdir=cbind(phi=90, lambda=0),
           projection=azimuthal.equalarea,
           proj.centre=cbind(phi=0, lambda=0),
           lambdalim=c(-180, 180),
           datapoints=TRUE, 
           datapoint.means=TRUE,
           datapoint.contours=TRUE, 
           grouped=FALSE,
           grouped.contours=FALSE,
           landmarks=TRUE,
           ids=getIDs(r),
           ...)
{
  ## This will call projection.reconstructedOutline(), but without
  ## drawing a grid.  The grid will be drawn later, after all the data
  ## has appeared..
  NextMethod(grid=FALSE, pole=FALSE)

  ## Datapoints
  if (datapoints) {
    Dss <- getDss(r)
    for (id in ids) {
      if (!is.null(Dss[[id]])) {
        suppressWarnings(points(projection(rotate.axis(transform(Dss[[id]],
                                                                 phi0=r$phi0),
                                                       axisdir*pi/180),
                                           proj.centre=pi/180*proj.centre),
                                col=r$cols[[id]],
                                pch=20, ...))
      }
    }
  }

  ## Mean datapoints
  if (datapoint.means) {
    Dss.mean <- getDssMean(r)
    for (id in ids) {
      if (!is.null(Dss.mean[[id]])) {
        suppressWarnings(points(projection(rotate.axis(transform(Dss.mean[[id]],
                                                                 phi0=r$phi0),
                                                       axisdir*pi/180),
                                           proj.centre=pi/180*proj.centre),
                                bg=r$cols[[id]], col="black",
                                pch=23, cex=1.5))
      }
    }
  }

  ## Groups
  if (grouped) {
    Gss <- getGss(r)
    for (id in ids) {
      if (!is.null(Gss[[id]])) {
        if (nrow(Gss[[id]]) > 0) {
          rc <- projection(rotate.axis(transform(Gss[[id]][,c("phi", "lambda")],
                                                 phi0=r$phi0),
                                       axisdir*pi/180),
                           proj.centre=pi/180*proj.centre)
        
          suppressWarnings(text(rc[,"x"], rc[,"y"], Gss[[id]][,3],
                                col=r$cols[[id]],
                                ...))
        }
      }
    }
  }
  
  ## KDE
  if (datapoint.contours) {
    k <- getKDE(r)
    for (id in ids) {
      if (!is.null(k[[id]])) {
        css <- k[[id]]$contours
        for(cs in css) {
          suppressWarnings(lines(projection(rotate.axis(transform(cs,
                                                                  phi0=r$phi0),
                                                        axisdir*pi/180),
                                            lambdalim=lambdalim*pi/180,
                                            lines=TRUE,
                                            proj.centre=pi/180*proj.centre),
                                 col=r$cols[[id]]))
        }
        ## FIXME: contours need to be labelled
      }
    }

    ## Plot locations of highest contours
    for (id in ids) {
      if (!is.null(k[[id]])) {
        suppressWarnings(points(projection(rotate.axis(transform(k[[id]]$maxs,
                                                                 phi0=r$phi0),
                                                       axisdir*pi/180),
                                           proj.centre=pi/180*proj.centre),
                                pch=22, cex=1, lwd=1,
                                col="black", bg=r$cols[[id]]))
      }
    }
  }
  
  ## KR
  if (grouped.contours) {
    k <- getKR(r)
    for (id in ids) {
      if (!is.null(k[[id]])) {
        css <- k[[id]]$contours
        for(cs in css) {
          suppressWarnings(lines(projection(rotate.axis(transform(cs,
                                                                  phi0=r$phi0),
                                                        axisdir*pi/180),
                                            lambdalim=lambdalim*pi/180,
                                            lines=TRUE,
                                            proj.centre=pi/180*proj.centre),
                                 col=r$cols[[id]]))
        }
        ## FIXME: contours need to be labelled
      }
    }
    ## Plot locations of highest contours
    for (id in ids) {
      if (!is.null(k[[id]])) {
        suppressWarnings(points(projection(rotate.axis(transform(k[[id]]$maxs,
                                                                 phi0=r$phi0),
                                                       axisdir*pi/180),
                                           proj.centre=pi/180*proj.centre),
                                pch=23, cex=1, lwd=1,
                                col="black", bg=r$cols[[id]]))
      }
    }
  }

  ## Landmarks
  if (landmarks) {
    Sss <- getSss(r)
    if (length(Sss)) {
      for (i in 1:length(Sss)) {
        name <- names(Sss)[i]
        col <- ifelse(is.null(name) || (name==""), "default", name)
        suppressWarnings(lines(projection(rotate.axis(transform(Sss[[i]],
                                                                phi0=r$phi0),
                                                      axisdir*pi/180),
                                          lines=TRUE,
                                          lambdalim=lambdalim*pi/180,
                                          proj.centre=pi/180*proj.centre),
                               col=r$cols[[col]], ...))
        ## FIXME: Might want to put lines argument here
      }
    }
  }

  ## This will call projection.reconstructedOutline() and will draw a
  ## grid but not add an image. Thus the grid appears over the data.
  NextMethod(add=TRUE,
             image=FALSE)
}

##' Draw a spherical plot of datapoints.
##'
##' @title Spherical plot of reconstructed outline
##' @param r \code{reconstructedOutline} object
##' @param datapoints If \code{TRUE}, display data points
##' @param ... Other graphics parameters -- not used at present
##' @method sphericalplot reconstructedDataset
##' @author David Sterratt
##' @export
sphericalplot.reconstructedDataset <- function(r,
                                               datapoints=TRUE,
                                               ...) {
  NextMethod()

  size <- r$R/10
  
  if (datapoints) {
    Dsc <- r$Dsc

    for (i in 1:length(Dsc)) {
      Dc <- Dsc[[i]]
      if (nrow(Dc) > 0) {
        
        ## Find axis in z=0 plane that is orthogonal to projection of
        ## datapoint onto that plane
        ax1 <- 1/sqrt(apply(Dc[,1:2,drop=FALSE]^2, 1, sum)) * cbind(-Dc[,2], Dc[,1], 0)

        ax1 <- as.matrix(ax1, ncol=3)

        ## Find axis that is orthogonal to the plane of axis 1 and the
        ## datapoint
        ax2 <- extprod3d(Dc, ax1)
        ax2 <- matrix(ax2, ncol=3)

        ax2 <- ax2/sqrt(apply(ax2^2, 1, sum))

        ## Create the verticies of an equillateral triangle to plot
        v1 <- Dc + size *  ax1/2
        v2 <- Dc + size * (-ax1/4 + sqrt(3)/4*ax2)
        v3 <- Dc + size * (-ax1/4 - sqrt(3)/4*ax2)

        ## Plot the triangle inside and outside the sphere
        inmag <- 0.99
        outmag <- 1.02
        
        x <- rbind(v2[,1], v1[,1], v3[,1])
        y <- rbind(v2[,2], v1[,2], v3[,2])
        z <- rbind(v2[,3], v1[,3], v3[,3])
        triangles3d(inmag*x, inmag*y, inmag*z, color=r$cols[[names(Dsc)[i]]])

        x <- rbind(v1[,1], v2[,1], v3[,1])
        y <- rbind(v1[,2], v2[,2], v3[,2])
        z <- rbind(v1[,3], v2[,3], v3[,3])
        triangles3d(outmag*x, outmag*y, outmag*z, color=r$cols[[names(Dsc)[i]]],
                    pch=20)
      }
    }
  }
}  

