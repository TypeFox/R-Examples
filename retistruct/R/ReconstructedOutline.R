##' Reconstruct outline into spherical surface. Reconstruction
##' proceeds in a number of stages:
##'
##' \enumerate{
##' 
##' \item The flat object is triangulated with at least \code{n}
##' triangles. This can introduce new vertices in the rim. 
##'
##' \item The triangulated object is stitched.
##'
##' \item The stitched object is triangulated again, but this time it
##' is not permitted to add extra vertices to the rim.
##'
##' \item The corresponding points determined by the stitching process
##' are merged to form a new set of merged points and a new
##' triangulation.
##'
##' \item The merged points are projected roughly to a sphere.
##'
##' \item The locations of the points on the sphere are moved so as to
##' minimise the energy function.
##' }
##'
##' @title Reconstruct outline into spherical surface
##' @param o \code{\link{AnnotatedOutline}} object, containing the following information:\describe{
##' \item{\code{P}}{outline points as N-by-2 matrix}
##' \item{\code{V0}}{indicies of the apex of each tear}
##' \item{\code{VF}}{indicies of the forward vertex of each tear}
##' \item{\code{VB}}{indicies of the backward vertex of each tear}
##' \item{\code{i0}}{index of the landmark on the rim}
##' \item{\code{phi0}}{lattitude of rim of partial sphere}
##' \item{\code{lambda0}}{longitude of landmark on rim}
##' }
##' @param n Number of points in triangulation.
##' @param alpha Area scaling coefficient
##' @param x0 Area cutoff coefficient
##' @param report Function used to report progress.
##' @param plot.3d Whether to show 3D picture during optimisation.
##' @param dev.flat Device to plot grid onto. Value of \code{NA} (default)
##' means no plotting.
##' @param dev.polar Device display projection. Value of NA
##' (default) means no plotting.
##' @return \code{reconstructedOutline} object containing the input
##' information and the following modified and extra information:
##' \item{\code{P}}{New set of points in flattened object}
##' \item{\code{gf}}{New set of forward pointers in flattened object}
##' \item{\code{gb}}{New set of backward pointers in flattened object}
##' \item{\code{phi}}{lattitude of new points on sphere}
##' \item{\code{lambda}}{longitude of new points on sphere}
##' \item{\code{Tt}}{New triangulation}
##' @author David Sterratt
ReconstructedOutline <- function(o, 
                                 n=500, alpha=8, x0=0.5,
                                 report=print,
                                 plot.3d=FALSE, dev.flat=NA, dev.polar=NA) {
  ## Clear polar plot, if it's required
  if (!is.na(dev.polar)) {
    dev.set(dev.polar)
    projection(o$phi0)
  }
  
  report("Triangulating...")
  t <- TriangulatedOutline(o, n=n)
  if (!is.na(dev.flat)) {
    dev.set(dev.flat)
    flatplot(t)
  }
    
  report("Stitching...")
  s <- StitchedOutline(t)
  if (!is.na(dev.flat)) {
    dev.set(dev.flat)
    flatplot(s, datapoints=FALSE)
  }

  report("Triangulating...")  
  r <- TriangulatedOutline(s, n=n,
                           suppress.external.steiner=TRUE)
  
  if (!is.na(dev.flat)) {
    dev.set(dev.flat)
    flatplot(r, datapoints=FALSE)
  }

  report("Merging points...")
  r <- mergePointsEdges(r)
  
  report("Projecting to sphere...")
  r <- projectToSphere(r)
  
  if (!is.na(dev.flat)) {
    ## Plot of initial gridlines
    dev.set(dev.flat)
      flatplot(r, grid=TRUE, strain=TRUE,
                datapoints=FALSE, landmarks=FALSE, mesh=FALSE, markup=FALSE)
    
    ## Initial plot in 3D space
    if (plot.3d) {
      sphericalplot(r)
    }
  }

  ## Check for flipped triangles and record initial number
  ft <- with(r, flipped.triangles(phi, lambda, Tt, R))
  r$nflip0 <- sum(ft$flipped)
  
  report("Optimising mapping with no area constratint using BFGS...")
  r <- optimiseMapping(r, alpha=0, x0=0, nu=1,
                        plot.3d=plot.3d, 
                        dev.flat=dev.flat, dev.polar=dev.polar)
  report("Optimising mapping with area constraint using FIRE...")
  r <- solveMappingCart(r, alpha=alpha, x0=x0, nu=1,
                          dtmax=500, maxmove=1E2, tol=1e-5,
                          plot.3d=plot.3d,
                          dev.flat=dev.flat, dev.polar=dev.polar)
  report("Optimising mapping with strong area constratint using BFGS...")
  r <- optimiseMapping(r, alpha=alpha, x0=x0, nu=1,
                        plot.3d=plot.3d,
                        dev.flat=dev.flat, dev.polar=dev.polar)
  report("Optimising mapping with weak area constratint using BFGS...")
  r <- optimiseMapping(r, alpha=alpha, x0=x0, nu=0.5,
                        plot.3d=plot.3d, 
                        dev.flat=dev.flat, dev.polar=dev.polar)
  
  report("Transforming image...")
  r <- transform.image.reconstructedOutline(r)
  
  report(paste("Mapping optimised. Deformation energy E:", format(r$opt$value, 5),
               ";", r$nflip, "flipped triangles."))
  class(r) <- addClass("reconstructedOutline", r)
  return(r)
}

##' Try a range of values of phi0s in the reconstruction, recording the
##' energy of the mapping in each case.
##'
##' @title Titrate values of phi0
##' @param r \code{\link{ReconstructedOutline}} object
##' @param alpha Area penalty scaling coefficient
##' @param x0 Area cutoff coefficient
##' @param byd Increments in degrees
##' @param len.up How many increments to go up from starting value of
##' \code{phi0} in \code{r}.
##' @param len.down How many increments to go up from starting value
##' of \code{phi0} in \code{r}.
##' @return dat Output data frame
##' @author David Sterratt
##' @export
titrate.reconstructedOutline <- function(r, alpha=8, x0=0.5, byd=1,
                                         len.up=5, len.down=20) {
  dat <- data.frame(phi0=r$phi0, sqrt.E=sqrt(r$E.l))

  by <- byd*pi/180

  ## Going up from phi0
  message("Going up from phi0")
  s <- r
  sqrt.E.min <- sqrt(r$E.l)
  r.opt <- r
  phi0s <- r$phi0 + seq(by, by=by, len=len.up)
  for (phi0 in phi0s)  {
    message(paste("phi0 =", phi0*180/pi))
    s$phi0 <- phi0
    ## Stretch the mapping to help with optimisation
    s$phi <- -pi/2 + (s$phi + pi/2)*(phi0+pi/2)/(s$phi0+pi/2)
    s <- optimiseMapping(s, alpha=alpha, x0=x0, nu=0.5,
                          plot.3d=FALSE)
    sqrt.E <- sqrt(s$E.l)
    dat <- rbind(dat, data.frame(phi0=s$phi0, sqrt.E=sqrt.E))
    if (sqrt.E < sqrt.E.min) {
      r.opt <- s
    }
  }

  ## Going down from phi0
  message("Going down from phi0")
  s <- r
  phi0s <- r$phi0 - seq(by, by=by, len=len.down)
  for (phi0 in phi0s)  {
    message(paste("phi0 =", phi0*180/pi))
    s$phi0 <- phi0
    ## Stretch the mapping to help with optimisation
    s$phi <- -pi/2 + (s$phi + pi/2)*(phi0+pi/2)/(s$phi0+pi/2)
    s <- optimiseMapping(s, alpha=alpha, x0=x0, nu=0.5,
                          plot.3d=FALSE)
    sqrt.E <- sqrt(s$E.l)
    dat <- rbind(dat, data.frame(phi0=s$phi0, sqrt.E=sqrt(s$E.l)))
    if (sqrt.E < sqrt.E.min) {
      r.opt <- s
    }
  }
  dat$phi0d <- dat$phi0*180/pi
  dat <- dat[order(dat$phi0d),]
  phi0d.opt <- dat[which.min(dat$sqrt.E),"phi0d"]

  ## Find mean difference between grid points
  ## First map range of original positions onto 
  phi.adj <- -pi/2 + (r$phi + pi/2)*(phi0d.opt*pi/180+pi/2)/(r$phi0+pi/2)
  Dtheta.mean <- mean(central.angle(phi.adj, r$lambda, r.opt$phi, r.opt$lambda)) * 180/pi
  
  return(list(dat=dat, phi0d.orig=r$phi0*180/pi,
              phi0d.opt=phi0d.opt,
              r.opt=r.opt,
              Dtheta.mean=Dtheta.mean))
}

##' This function returns information about how edges on the sphere
##' have been deformed from their flat state.
##'
##' @title Return strains edges are under in spherical retina
##' @param r A \code{\link{ReconstructedOutline}} object
##' @return A list containing two data frames \code{flat} and \code{spherical}. 
##' Each data frame contains for each edge in the flat or spherical meshes:
##' \item{\code{L}}{Length of the edge in the flat outline }
##' \item{\code{l}}{Length of the corresponding edge on the sphere}
##' \item{\code{strain}}{The strain of each connection}
##' \item{\code{logstrain}}{The logarithmic strain of each connection}
##' @author David Sterratt
getStrains <- function(r) {
  ## Original lengths in flattened outline is a vector with
  ## M elements, the number of rows of Cu
  L <- r$L
  ## New lengths in reconstructed object is a vector wtih Mt < M
  ## elements, the number of rows of Cut
  lt <- compute.lengths(r$phi, r$lambda, r$Cut, r$R)
  ## For each connection in the flattened object, we want the length of
  ## the corresponding connection in the reconstructed object
  ## The mapping Ht achieves this
  l <- lt[r$Ht]
  stretch <- l/L
  strain <- stretch - 1
  logstrain <- log(stretch)

  ## Compute quantities in spherical retina too
  Lt <- r$Lt
  stretcht <- lt/Lt
  straint <- stretcht - 1
  logstraint <- log(stretcht)

  return(list(flat=
              data.frame(L=L,  l=l,
                         strain=strain,  logstrain=logstrain),
              spherical=
              data.frame(L=Lt, l=lt,
                         strain=straint, logstrain=logstraint)))
}

##' Transform an image into the reconstructed space. The four corner
##' coordinates of each pixel are transformed into spherical
##' coordinates and a mask matrix with the same dimensions as
##' \code{im} is created. This has \code{TRUE} for pixels that should
##' be displayed and \code{FALSE} for ones that should not.
##'
##' @title Transform an image into the reconstructed space
##' @param r \code{reconstructedOutline} object
##' @return \code{reconstructedOutline} object with extra elements
##' \item{\code{ims}}{Coordinates of corners of pixes in spherical coordinates}
##' \item{\code{immask}}{Mask matrix with same dimensions as image \code{im}}
##' @author David Sterratt
transform.image.reconstructedOutline <- function(r) {
  if (!is.null(r$im)) {
    ## Need to find the *boundaries* of pixels
    N <- ncol(r$im)
    M <- nrow(r$im)

    ## Create grid coords of corners of pixels.  These run from the
    ## top left of the image down each column of the image.
    xs <- 0:N
    ys <- M:0
    ## x-coords of pixel corners, arranged in (N+1) by (M+1) grid 
    Ix <- outer(ys*0, xs, FUN="+")
    ## Ditto for y-coords
    Iy <- outer(ys, xs*0, FUN="+")
    ## Join to give (x, y) coordinates of all corners
    I  <- cbind(as.vector(Ix), as.vector(Iy))
    
    ## Find Barycentric coordinates of corners of pixels
    Ib <- tsearchn(r$P, r$T, I)
    
    ## Create mask depending on whether corners are in outline
    idx <- matrix(Ib$idx, M+1, N+1)
    r$immask <- (!is.na(idx[1:M    , 1:N    ]) & 
                 !is.na(idx[1:M    , 2:(N+1)]) &
                 !is.na(idx[2:(M+1), 1:N    ]) &
                 !is.na(idx[2:(M+1), 2:(N+1)]))
    
    ## Convert to spherical coordinates
    Ic <- with(r, bary.to.sphere.cart(phi, lambda, R, Tt, Ib))
    r$ims <- with(r, sphere.cart.to.sphere.spherical(Ic, R))
  }
  return(r)
}
##' @title Get coordinates of corners of pixels of image in spherical
##' coordinates 
##' @param r \code{\link{ReconstructedOutline}} object
##' @return Coordinates of corners of pixels in spherical coordinates 
##' @author David Sterratt
##' @method getIms reconstructedOutline
##' @export
getIms.reconstructedOutline <- function(r) {
  return(r$ims)
}

##' @export
getTss.reconstructedOutline <- function(r) {
  Tss <- list()
  for (TF in r$TFset) {
    ## Convert indicies to the spherical frame of reference
    j <- r$ht[TF]
    Tss <- with(r, c(Tss, list(cbind(phi=phi[j], lambda=lambda[j]))))
  }
  return(Tss)
}

##' Plot \code{\link{ReconstructedOutline}} object. This adds a mesh
##' of gridlines from the spherical retina (described by points
##' \code{phi}, \code{lambda} and triangulation \code{Tt} and cutoff
##' point \code{phi0}) onto a flattened retina (described by points
##' \code{P} and triangulation \code{T}).
##'
##' @title Flat plot of reconstructed outline
##' @param x \code{\link{ReconstructedOutline}} object
##' @param axt whether to plot axes
##' @param ylim y-limits
##' @param grid Whether or not to show the grid lines of
##' lattitude and longitude
##' @param strain Whether or not to show the strain
##' @param ... Other plotting parameters
##' @method flatplot reconstructedOutline
##' @author David Sterratt
##' @export
flatplot.reconstructedOutline <- function(x, axt="n", ylim=NULL,
                                          grid=TRUE,
                                          strain=FALSE,
                                          ...) {
  NextMethod()

  if (strain) {
    o <- getStrains(x)
    palette(rainbow(100))
    scols <- strain.colours(o$flat$logstrain)
    with(x, 
         segments(P[Cu[,1],1], P[Cu[,1],2],
                  P[Cu[,2],1], P[Cu[,2],2], col=round(scols)))
  }
  
  ## Plot a gridline from the spherical retina (described by points phi,
  ## lambda and triangulation Tt) onto a flattened retina (described by
  ## points P and triangulation T). The gridline is described by a
  ## normal n to a plane and a distance to the plane. The intersection of
  ## the plane and the spehere is the gridline.
  get.gridline.flat <- function(P, T, phi, lambda, Tt, n, d, ...) {
    mu <- compute.intersections.sphere(phi, lambda, Tt, n, d)

    ## Take out rows that are not intersections. If a plane intersects
    ## one edge of a triangle and the opposing vertex, in the row
    ## corresponding to the triangle, there will be a 0, a 1 and a
    ## value between 0 and 1. We get rid of the 1 in the
    ## following. Triangles in which one line is in the plane have mu
    ## values 0, 1 and NaN; we want to include these.
    tri.int <- (rowSums((mu >= 0) & (mu <= 1), na.rm=TRUE) == 2) 
    ## | apply(mu, 1, function(x) setequal(x, c(0, 1, NaN))))
    if (any(tri.int)) {
      T  <- T[tri.int,,drop=FALSE]
      mu <- mu[tri.int,,drop=FALSE]

      ## Create a logical matrix of which points are involved in lines
      ## that interscect the plane.
      line.int <- (mu >= 0) & (mu < 1)
      ## If any element of mu contained a NaN, due to a line being in
      ## the plane, this should be set to false as the point opposite
      ## the NaN is not in the plane
      line.int[is.na(line.int)] <- FALSE

      ## Order rows so that the false indicator is in the third column
      T[!line.int[,2] ,] <- T[!line.int[,2], c(3,1,2)]
      mu[!line.int[,2],] <- mu[!line.int[,2],c(3,1,2)]
      T[!line.int[,1] ,] <- T[!line.int[,1], c(2,3,1)]
      mu[!line.int[,1],] <- mu[!line.int[,1],c(2,3,1)]

      P <- cbind(mu[,1] * P[T[,3],] + (1-mu[,1]) * P[T[,2],], 
                 mu[,2] * P[T[,1],] + (1-mu[,2]) * P[T[,3],])
      # suppressWarnings(segments(P1[,1], P1[,2], P2[,1], P2[,2], ...))
    } else {
      P <- matrix(0, nrow=1, ncol=4)
    }
    colnames(P) <- c("X1", "Y1", "X2", "Y2")
    return(P)
  }
  
  if (grid) {
    grid.int.minor <- 15
    grid.int.major <- 45
    grid.maj.col <- getOption("grid.maj.col")
    grid.min.col <- getOption("grid.min.col")

    phi0d <- x$phi0 * 180/pi
    
    P <- matrix(0, nrow=0, ncol=4)
    cols <- NULL
    Phis <- seq(-90, phi0d, by=grid.int.minor)
    Lambdas <- seq(0, 180-grid.int.minor, by=grid.int.minor)
    for (Phi in Phis) {
      if ((!(Phi %% grid.int.major) || Phi == phi0d)) {
        col <- grid.maj.col
      } else {
        col <- grid.min.col
      }
      P1 <- get.gridline.flat(x$P, x$T, x$phi, x$lambda, x$Tt,
                              c(0,0,1), sin(Phi*pi/180))
      cols <- c(cols, rep(col, nrow(P1)))
      P <- rbind(P, P1)
    }
    for (Lambda in Lambdas) {
      if (!(Lambda %% grid.int.major)) {
        col <- grid.maj.col
      } else {
        col <- grid.min.col
      }
      Lambda <- Lambda * pi/180
      P1 <- get.gridline.flat(x$P, x$T, x$phi, x$lambda, x$Tt,
                              c(sin(Lambda),cos(Lambda),0), 0)
      cols <- c(cols, rep(col, nrow(P1)))
      P <- rbind(P, P1)
    }
    if (nrow(P) > 0) {
      segments(P[,"X1"], P[,"Y1"], P[,"X2"], P[,"Y2"], col=cols)
    }
  }
}

##' Draw a projection of a \code{\link{ReconstructedOutline}}. This method sets up
##' the grid lines and the angular labels and draws the image.
##'
##' @title  Projection of a reconstructed outline
##' @param r \code{ReconstructedOutline} object
##' @param transform Transform function to apply to spherical coordinates
##' before rotation
##' @param axisdir Direction of axis (North pole) of sphere in
##' external space as matrix with column names \code{phi} (elevation)
##' and \code{lambda} (longitude).
##' @param projection Projection in which to display object,
##' e.g. \code{\link{azimuthal.equalarea}} or \code{\link{sinusoidal}}
##' @param proj.centre Location of centre of projection as matrix with
##' column names \code{phi} (elevation) and \code{lambda} (longitude).
##' @param lambdalim Limits of longitude (in degrees) to display
##' @param philim Limits of latitude (in degrees) to display
##' @param labels Vector of 4 labels to plot at 0, 90, 180 and 270 degrees 
##' @param grid Whether or not to show the grid lines of
##' lattitude and longitude
##' @param grid.bg Background colour of the grid
##' @param grid.int.minor Interval between minor grid lines in degrees
##' @param grid.int.major Interval between major grid lines in degrees
##' @param colatitude If \code{TRUE} have radial labels plotted with
##' respect to colatitude rather than latitude
##' @param pole If \code{TRUE} indicate the pole with a "*"
##' @param image If \code{TRUE}, show the image
##' @param add If \code{TRUE}, don't draw axes; add to existing plot.
##' @param ... Graphical parameters to pass to plotting functions
##' @method projection reconstructedOutline
##' @export
projection.reconstructedOutline <- function(r,
                                            transform=identity.transform,
                                            axisdir=cbind(phi=90, lambda=0),
                                            projection=azimuthal.equalarea,
                                            proj.centre=cbind(phi=0, lambda=0),
                                            lambdalim=c(-180, 180),
                                            philim=c(-90, 90),          
                                            labels=c(0, 90, 180, 270),
                                            grid=TRUE,
                                            grid.bg="transparent", 
                                            grid.int.minor=15,
                                            grid.int.major=45,
                                            colatitude=TRUE,
                                            pole=FALSE,
                                            image=TRUE,
                                            add=FALSE,
                                            ...) {
  plot.image <- image
  ## Compute grid lines

  ## Lines of latitude (parallels)

  ## Determine the major and minor parallels
  phis.maj <- seq(-90, 90, by=grid.int.major)
  phis.maj <- c(philim[1],
                phis.maj[(phis.maj > philim[1]) & (phis.maj < philim[2])],
                philim[2])
  phis.min <- seq(-90, 90, by=grid.int.minor)
  phis.min <- c(philim[0],
                phis.min[(phis.min > philim[1]) & (phis.min < philim[2])],
                philim[2])
  phis.min <- setdiff(phis.min, phis.maj)

  ## Longitudes at which to draw lines; the smaller the by interval,
  ## the smoother
  lambdas <- seq(lambdalim[1], lambdalim[2], by=1)

  ## Compute the minor and and major parallels to draw
  paras.min <- projection(pi/180*cbind(phi   =as.vector(outer(c(lambdas, NA)*0, phis.min, FUN="+")),
                                       lambda=as.vector(outer(c(lambdas, NA), phis.min*0, FUN="+"))),
                          proj.centre=pi/180*proj.centre)
  paras.maj <- projection(pi/180*cbind(phi   =as.vector(outer(c(lambdas, NA)*0, phis.maj, FUN="+")),
                                       lambda=as.vector(outer(c(lambdas, NA), phis.maj*0, FUN="+"))),
                          proj.centre=pi/180*proj.centre)

  ## Lines of longitude (meridians)

  ## Determine the major and minor parallels
  lambdas.maj <- c(rev(seq(0         , lambdalim[1], by=-grid.int.major)),
                   seq(grid.int.major, lambdalim[2], by= grid.int.major))
  lambdas.min <- c(rev(seq(0         , lambdalim[1], by=-grid.int.minor)),
                   seq(grid.int.minor, lambdalim[2], by= grid.int.minor))
  lambdas.min <- setdiff(lambdas.min, lambdas.maj)

  ## Lattidues at which to draw lines; the smaller the by interval,
  ## the smoother
  phis <- seq(philim[1], philim[2], by=1)

  ## Compute the minor and and major meridians to draw
  merids.min <- projection(pi/180*cbind(phi   =as.vector(outer(c(phis, NA), lambdas.min*0, FUN="+")),
                                        lambda=as.vector(outer(c(phis*0, NA), lambdas.min, FUN="+"))),
                           proj.centre=pi/180*proj.centre)
  merids.maj <- projection(pi/180*cbind(phi   =as.vector(outer(c(phis, NA), lambdas.maj*0, FUN="+")),
                                        lambda=as.vector(outer(c(phis*0, NA), lambdas.maj, FUN="+"))),
                           proj.centre=pi/180*proj.centre)
  
  ## Set up the plot region
  xlim <- range(na.omit(rbind(paras.min, paras.maj))[,"x"])
  ylim <- range(na.omit(rbind(paras.min, paras.maj))[,"y"])

  if (!add) {
    plot(NA, NA, xlim=xlim, ylim=ylim,# xaxs="i", yaxs="i",
         type = "n", axes = FALSE, xlab = "", ylab = "", asp=1)
  }

  ## Plot an image

  ## Get the spherical coordinates of the corners of pixels
  ims <- getIms(r)
  if (plot.image && !is.null(ims)) {
    ## Reconstitute image from stored values of phi and lambda
    ## coordinates of corners of pixels

    ## Get the size of the image
    M <- nrow(r$im)
    N <- ncol(r$im)

    ## Downsample the image by first selecting rows and columns to
    ## look at
    max.proj.dim <- getOption("max.proj.dim")
    by <- ceiling(max(N, M)/max.proj.dim) # Number of pixels to merge
    Ms <- seq(1, M - (M %% by), by=by)
    Ns <- seq(1, N - (N %% by), by=by)

    ## Downsample the image
    im <- r$im[Ms, Ns]
    
    ## Now need to do the more complex job of downsampling the matrix
    ## containing the coordinates of the corners of pixels
    imsmask <- matrix(FALSE, M+1, N+1)
    imsmask[c(Ms, (max(Ms) + by)), c(Ns, (max(Ns) + by))] <- TRUE
    ims <- ims[imsmask,]

    ## Convenience variables for the new image sizes
    M <- nrow(im)
    N <- ncol(im)
    
    ## Transform the pixel coordinates and compute x and y positions
    ## of corners of pixels.
    rc <- projection(rotate.axis(transform(ims, phi0=r$phi0),
                                 axisdir*pi/180),
                     lambdalim=lambdalim*pi/180,
                     proj.centre=pi/180*proj.centre)
    xpos <- matrix(rc[,"x"], M+1, N+1)
    ypos <- matrix(rc[,"y"], M+1, N+1)
    
    ## Convert these to format suitable for polygon()
    impx <- rbind(as.vector(xpos[1:M    , 1:N    ]),
                  as.vector(xpos[1:M    , 2:(N+1)]),
                  as.vector(xpos[2:(M+1), 2:(N+1)]),
                  as.vector(xpos[2:(M+1), 1:N    ]),
                  NA)
    impy <- rbind(as.vector(ypos[1:M    , 1:N    ]),
                  as.vector(ypos[1:M    , 2:(N+1)]),
                  as.vector(ypos[2:(M+1), 2:(N+1)]),
                  as.vector(ypos[2:(M+1), 1:N    ]),
                  NA)

    ## Pixels outside the image should be masked. The mask is a matrix
    ## the same size as the image, containing TRUE for pixels that
    ## should be displayed and FALSE for those that should be
    ## masked. It is computed by finding the corners of the poly-pixel
    ## lie outwith the outline. These corners will have the coordinate
    ## NA.  print("sum(!is.na(colSums(impx[1:4,])))")
    ## print(sum(!is.na(colSums(impx[1:4,]))))
    immask <- matrix(!is.na(colSums(impx[1:4,])), M, N)
    
    ## We want to get rid of any poly-pixels that cross either end of
    ## the longitude range in a pseudocylindrical projection. A simple
    ## way of doing this is to say that if a pixel is very large,
    ## don't plot it.
    
    ## This code is chronically slow, hence replacing with the
    ## pmin/pmax version
    ## bigpx <- which(apply(impx[1:4,], 2,
    ##                      function(x) {max(x) - min(x)}) > 0.1*abs(diff(xlim)) |
    ##                apply(impy[1:4,], 2,
    ##                      function(y) {max(y) - min(y)}) > 0.1*abs(diff(ylim)))
    bigpx <- which((pmax(impx[1,], impx[2,], impx[3,], impx[4,]) -
                    pmin(impx[1,], impx[2,], impx[3,], impx[4,]) >
                    0.1*abs(diff(xlim))) |
                   (pmax(impy[1,], impy[2,], impy[3,], impy[4,]) -
                    pmin(impy[1,], impy[2,], impy[3,], impy[4,]) >
                    0.1*abs(diff(xlim))))
    immask[bigpx] <- FALSE
    
    ## Plot the polygon, masking as we go
    polygon(impx[,immask], impy[,immask],
            col=im[immask], border=im[immask])
  }

  grid.maj.col <- getOption("grid.maj.col")
  grid.min.col <- getOption("grid.min.col")
  
  ## Plot the grid
  if (grid) {
    ## Minor paralells and meridians
    lines(paras.min,  col=grid.min.col)
    lines(merids.min, col=grid.min.col)

    ## Major lines of latitude on top of all minor lines
    lines(paras.maj,  col=grid.maj.col)
    lines(merids.maj, col=grid.maj.col)

    ## Boundary of projection
    boundary <- projection("boundary")
    polygon(boundary[,"x"], boundary[,"y"], border="black")
  }

  ## Plot rim in visutopic space
  rs <- cbind(phi=r$phi0, lambda=seq(0, 2*pi, len=360))
  rs.rot <- rotate.axis(transform(rs, phi0=r$phi0), axisdir*pi/180)
  ## "Home" position for a cyclops looking ahead
  ## r$axisdir = cbind(phi=0, lambda=0)
  
  lines(projection(rs.rot, lambdalim=lambdalim*pi/180, lines=TRUE,
                   proj.centre=pi/180*proj.centre),
        col=getOption("TF.col"))

  ## Projection of pole
  if (pole) {
    oa.rot <- rotate.axis(transform(cbind(phi=-pi/2, lambda=0), phi0=r$phi0),
                          axisdir*pi/180)
    points(projection(oa.rot, proj.centre=pi/180*proj.centre),
           pch="*", col=getOption("TF.col"), cex=2)
  }

  ## Plot outline
  Tss <- getTss(r)
  for (Ts in Tss) {
    ## Plot
    suppressWarnings(lines(projection(rotate.axis(transform(Ts, phi0=r$phi0),
                                                  axisdir*pi/180),
                                      lines=TRUE,
                                      lambdalim=lambdalim*pi/180,
                                      proj.centre=pi/180*proj.centre),
                           col=getOption("TF.col"), ...))
  }
  
  ## Longitude labels around rim - not on actual frame of reference!
  if (!is.null(labels)) {
    ## Longitudes (meridians) at which to plot at
    angles <- seq(0, by=2*pi/length(labels), len=length(labels))

    ## This is a nasty hack: we want to find out how far to plot the
    ## labels from the rim. We can't do this simply in terms of
    ## angles, so we have to find the right angular distance to give
    ## the desired fraction of the axes at which to plot the
    ## labels. This is done by this optimisation function.
    label.fax <- 0.02                   # Fraction of axis length from axes to plot labels
    opt <- optimise(function(a) {
      rs0 <- cbind(phi=r$phi0,     lambda=angles[1])
      rs  <- cbind(phi=r$phi0 + a, lambda=angles[1])
      rc0 <- projection(rotate.axis(transform(rs0, phi0=r$phi0),
                                    axisdir*pi/180),
                        proj.centre=pi/180*proj.centre)
      rc  <- projection(rotate.axis(transform(rs, phi0=r$phi0),
                                    axisdir*pi/180),
                        proj.centre=pi/180*proj.centre)
      return((vecnorm(rc - rc0) - label.fax*abs(diff(xlim)))^2)
    }
                                ,interval=c(1, 20)*pi/180)
    lambda.label.off <- opt$minimum

    ## Now plot the labels themselves. Phew!!
    rs <- cbind(phi=r$phi0 + lambda.label.off, lambda=angles)
    rc <- projection(rotate.axis(transform(rs, phi0=r$phi0), axisdir*pi/180),
                     proj.centre=pi/180*proj.centre)
    text(rc[,"x"], rc[,"y"], labels, xpd=TRUE)
  }

  ## Lattitude Labels
  ## rlabels <- c(seq(philim[1], philim[2], by=grid.int.major))
  rlabels <- phis.maj
  rs <- cbind(phi=rlabels*pi/180, lambda=proj.centre[1,"lambda"])
  rc <- projection(rs, proj.centre=pi/180*proj.centre)
  text(rc[,"x"], rc[,"y"], rlabels + ifelse(colatitude, 90, 0),
       xpd=TRUE, adj=c(1, 1), col=grid.maj.col)
}

##' Draw a spherical plot of reconstructed outline. This method just
##' draws the mesh.
##'
##' @title Spherical plot of reconstructed outline
##' @param r \code{\link{ReconstructedOutline}} object
##' @param strain If \code{TRUE}, plot the strain
##' @param surf If \code{TRUE}, plot the surface
##' @param ... Other graphics parameters -- not used at present
##' @method sphericalplot reconstructedOutline
##' @author David Sterratt
##' @import rgl
##' @export
sphericalplot.reconstructedOutline <- function(r,
                                               strain=FALSE,
                                               surf=TRUE, ...) {
  NextMethod()
  
  ## FIXME: This needs to be looked at with a view to replacing
  ## functions in plots.R
  with(r, {
    ## Obtain Cartesian coordinates of points
    P <- sphere.spherical.to.sphere.cart(phi, lambda, R)

    if (surf) {
      ## Outer triangles
      fac <- 1.005
      triangles3d(matrix(fac*P[t(Tt[,c(2,1,3)]),1], nrow=3),
                  matrix(fac*P[t(Tt[,c(2,1,3)]),2], nrow=3),
                  matrix(fac*P[t(Tt[,c(2,1,3)]),3], nrow=3),
                  color="darkgrey", alpha=1)
      
      ## Inner triangles
      triangles3d(matrix(P[t(Tt),1], nrow=3),
                  matrix(P[t(Tt),2], nrow=3),
                  matrix(P[t(Tt),3], nrow=3),
                  color="white", alpha=1)
    }
    
    ## Plot any flipped triangles
    ft <- flipped.triangles(phi, lambda, Tt, R)
    with(ft, points3d(cents[flipped,1], cents[flipped,2], cents[flipped,3],
                      col="blue", size=5))

    ## Shrink so that they appear inside the hemisphere
    fac <- 0.997
    rgl.lines(fac*rbind(P[ht[gb[gb]],1], P[ht[gb],1]),
              fac*rbind(P[ht[gb[gb]],2], P[ht[gb],2]),
              fac*rbind(P[ht[gb[gb]],3], P[ht[gb],3]),
              lwd=3, color=getOption("TF.col"))
    
    fac <- 1.006
    rgl.lines(fac*rbind(P[ht[gb[gb]],1], P[ht[gb],1]),
              fac*rbind(P[ht[gb[gb]],2], P[ht[gb],2]),
              fac*rbind(P[ht[gb[gb]],3], P[ht[gb],3]),
              lwd=3, color=getOption("TF.col"))

    if (strain) {
      o <- getStrains(r)
      palette(rainbow(100))
      scols <- strain.colours(o$spherical$logstrain)

      fac <- 0.999
      P1 <- fac*P[Cut[,1],]
      P2 <- fac*P[Cut[,2],]
      
      width <- 40
      ## Compute displacement vector to make sure that strips are
      ## parallel to surface of sphere
      d <- extprod3d(P1, P2-P1)
      d <- width/2*d/vecnorm(d)
      PA <- P1 - d
      PB <- P1 + d
      PC <- P2 + d
      PD <- P2 - d

      ## This is a ridiculously inefficient way of drawing the strain,
      ## but if you try presenting a color vector, it makes each line
      ## multi-coloured. It has taking HOURS of fiddling round to
      ## discover this! GRRRRRRRRRRRRRRRRR!
      for (i in 1:nrow(PA)) {
        quads3d(rbind(PA[i,1], PB[i,1], PC[i,1], PD[i,1]),
                rbind(PA[i,2], PB[i,2], PC[i,2], PD[i,2]),
                rbind(PA[i,3], PB[i,3], PC[i,3], PD[i,3]),
                color=round(scols[i]), alpha=1)
      }
    }
  })
}

##' @export
lvsLplot.reconstructedOutline <- function(r) {
  o <- getStrains(r)$spherical
  palette(rainbow(100)) ## Green is about 35; dark blue about 70
  cols <- strain.colours(o$logstrain)
  L <- o$L
  l <- o$l
  
  if (!is.null(r$scale)) {
    scale <- r$scale
    if (!is.na(scale)) {
      L <- L*scale
      l <- l*scale
    }
  }

  plot(L, l, col=cols, pch=20,
       xlim=c(0, max(L, l)), ylim=c(0, max(L, l)),
       xlab=ifelse(!is.na(scale),
         expression(paste("Flat ", italic(L), " (", mu, "m)")),
         expression(paste("Flat ", italic(L)))),
       ylab=ifelse(!is.na(scale),
         expression(paste("Reconstructed ", italic(l), " (", mu, "m)")),
         expression(paste("Reconstructed ", italic(l)))),
       asp=1)
  par(xpd=FALSE)
  abline(0, 1)
  abline(0, 0.75, col="blue")
  abline(0, 1.25, col="red")
  text(0.2*max(L), 0.2*max(L)*0.5, "25% compressed", col="blue",
               pos=4)
  text(0.75*max(L), 0.75*max(L)*1.25, "25% expanded", col="red",
               pos=2)
}
