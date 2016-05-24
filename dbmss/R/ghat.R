ghat <-
function(X, r = NULL, ReferenceType = "", NeighborType = "", CheckArguments = TRUE) {

  if (CheckArguments) {
    CheckdbmssArguments()
    # Eliminate erroneous configurations
    if ((ReferenceType=="" | NeighborType=="") & (ReferenceType!=NeighborType)) {
      stop("Either two or no point type must be specified.")
    }
  }

  area <- area.owin(X$window)
  autor <- is.null(r)

  # Calculate densities
  # g intra
  if (ReferenceType=="" & NeighborType=="") {
    lambdaI <- lambdaJ <- X$n/area
  } else {
    # g intra for a single point type
    if (ReferenceType==NeighborType) {
      X.reduced <- X[X$marks$PointType==ReferenceType]
      lambdaI <- lambdaJ <- X.reduced$n/area
    }
    # g inter
    if (ReferenceType!=NeighborType) {
      X.cross <- X[X$marks$PointType==ReferenceType]
      Y.cross <- X[X$marks$PointType==NeighborType]
      lambdaI <- X.cross$n/area
      lambdaJ <- Y.cross$n/area
    }
  }

  # Find the best breaks
  rmaxdefault <- rmax.rule("K", X$window, lambdaJ)
  breaks <- handle.r.b.args(window = X$window, rmaxdefault = rmaxdefault)
  rBest <- breaks$r
  denargs <- resolve.defaults(list(kernel = "epanechnikov",  n = length(rBest), from = 0, to = max(rBest)))
  if (is.null(r)) {
    r <- rBest
  }
  
  # Find pairs
  # g intra
  if (ReferenceType=="" & NeighborType=="") {
    Pairs <- closepairs(X, max(r))
  } else {
    # g intra for a single point type
    if (ReferenceType==NeighborType) {
      Pairs <- closepairs(X.reduced, max(r))
    }
    # g inter
    if (ReferenceType!=NeighborType) {
      Pairs <- crosspairs(X.cross, Y.cross, max(r))
    }
  }

  # Adapted from pcf.ppp {spatstat}

  # Geometry
  XI <- ppp(Pairs$xi, Pairs$yi, window = X$window, check = FALSE)
  XJ <- ppp(Pairs$xj, Pairs$yj, window = X$window, check = FALSE)
  
  # Edge-effect correction
  if (is.rectangle(X$window) |  is.polygonal(X$window)) {
    edgewt <- edge.Ripley(XI, matrix(Pairs$d, ncol = 1))
    valu <- "iso"
    desc <- "Ripley isotropic correction estimate of %s"
  } else {
    edgewt <- edge.Trans(XI, XJ, paired = TRUE)
    valu <- "trans"    
    desc <- "translation-corrected estimate of %s"
  }

  # Estimate g  
  gEstimate <- sewpcf(Pairs$d, edgewt, denargs, lambdaI*lambdaJ*area)
  # Calculate values for r if specified
  if (!autor) {
    g <- approx(gEstimate$r, gEstimate$g, xout=r)$y
    gEstimate <- data.frame(r, g)
  }
  # Add theoretical value
  theo <- rep(1, length(r))
  gEstimate <- data.frame(gEstimate[, 1], theo, gEstimate[, 2])
  colnames(gEstimate) <- c("r", "theo", valu)
  
  # Return the values of g(r)
  return (fv(gEstimate, argu="r", ylab=quote(g(r)), valu=valu, fmla= ". ~ r", alim=c(0, max(r)), labl=c("r", "%s[pois](r)", paste("hat(%s)[", valu, "](r)", sep="")), desc=c("distance argument r", "theoretical Poisson g(r)", desc), unitname=X$window$unit, fname="g"))
}
