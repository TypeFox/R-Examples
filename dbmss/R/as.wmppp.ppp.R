as.wmppp.ppp <-
  function(X, ...) {

    Marks <- marks(X)
    
    # Marks are a dataframe with colums PointType and PointWeight
    if (is.data.frame(Marks)) {
      names(Marks) <- tolower(names(Marks))
      # Read Point Weights
      if ("pointtype" %in% names(Marks)) {
        PointType <- Marks[, "pointtype"]
      }
      if (!is.factor(PointType))  {
        PointType <- as.factor(PointType)
      }    
      if ("pointweight" %in% names(Marks)) {
        PointWeight <- Marks[, "pointweight"]
      }      
      if (!is.numeric(PointWeight))  {
        stop("Point weights have been found in the marks of the point pattern but they must be numeric.")
      }    
    } else {
      # Marks are types, set types to "All"  
      if (is.factor(Marks)) {
        PointType <- Marks
      } else {
        PointType <- rep("All", X$n)
      }
      # Marks are numbers: use them as weights, else set weights to 1
      if (is.numeric(Marks)) {
        PointWeight <- Marks
      } else {
        PointWeight <- rep(1, X$n)
      }
    }
    # Check weights are positive
    if (any(PointWeight < 0)) {
      stop("Point weights must be positive.")
    }
    
    wmpppX <- X
    marks(wmpppX) <- data.frame(PointWeight, PointType)
    class(wmpppX) <- c("wmppp", "ppp")
    return (wmpppX)
  }
