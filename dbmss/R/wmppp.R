wmppp <-
function(df, window = NULL, unitname = NULL) {
  
  if (!is.data.frame(df))
    stop("The data used to create a wmppp must be a dataframe.")
  if (ncol(df) < 2)
    stop("The data used to create a wmppp must have at least two columns for coordinates X and Y.")
  names(df) <- tolower(names(df))
  # Read X and Y
  if ("x" %in% names(df) & "y" %in% names(df)) {
    X <- as.numeric(df[, "x"])
    Y <- as.numeric(df[, "y"])
  } else {
    warning("No columns named X and Y have been found. Columns #1 and #2 have been used for coordinates.")
    X <- as.numeric(df[, 1])
    Y <- as.numeric(df[, 2])
  }
  if (!is.numeric(c(X, Y)))
    stop("Point coordinates X and Y must be numeric.")
  
  # Read Point Types
  if ("pointtype" %in% names(df)) {
    PointType <- df[, "pointtype"]
  } else {
    if (ncol(df) < 3) {
      warning("No column has been found for PointType. All point types have been set to All.")      
      PointType <- as.factor(rep("All", length(X)))
    } else {
      warning("No column named PointType has been found. Columns #3 has been used for labels.")
      PointType <- df[, 3]
    }
  }
  if (!is.factor(PointType))
    stop("Point types must be factors.")
  
  # Read Point Weights
  if ("pointweight" %in% names(df)) {
    PointWeight <- df[, "pointweight"]
  } else {
    if (ncol(df) < 4) {
      warning("No column has been found for PointWeight. All point weights have been set to 1.")      
      PointWeight <- rep(1, length(X))
    } else {
      warning("No column named PointWeight has been found. Columns #4 has been used for weights.")
      PointWeight <- df[, 4]
    }
  }
  if (!is.numeric(PointWeight))
    stop("Point weights must be numeric.")
  if (any(PointWeight < 0))
    stop("Point weights must be positive.")
  
  # Check the window
  if (is.null(window)) {
    window <- owin(xrange=c(min(X), max(X)), yrange=c(min(Y), max(Y)), unitname=unitname)        
    warning("No window has been specified. A rectangle window containing all points has been used.")
  }
  if (!is.owin(window))
    stop("window must be an object of class owin.")
  wmpppX <- ppp(X, Y, window=window, marks=data.frame(PointWeight, PointType))
  class(wmpppX) <- c("wmppp", "ppp")
  return (wmpppX)
}
