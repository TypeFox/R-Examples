estThetaRow <-
function(X, f, theta)
  {
    if (any(floor(X) != X))
      {
        warning("The elements of X should be integers !")
      }
    J <- length(f)
    if (nrow(X) != J)
      {
        stop("The number of rows in X must be equal to the length of f !")
      }
    X <- as.data.frame(X)
    ff <- as.factor(f)
    inds <- apply(X = X, MARGIN = 1, FUN = theta)
    indsdat <- data.frame(x = inds, y = ff)
    names(indsdat) <- c("inds", "fac")
    return(dat = indsdat)
  }

