ymatrix.mnp <- function(data, base=NULL, extra=FALSE, verbose=verbose) { 
  ## checking and formatting Y
  Y <- model.response(data)
  if (is.matrix(Y)) { # Multinomial ordered Probit model
    for (i in 1:nrow(Y))
      Y[i,] <- match(Y[i,], sort(unique(Y[i,]))) - 1
    p <- ncol(Y)
    lev <- colnames(Y)
    MoP <- TRUE
    if(!is.null(base))
      stop("Error: The last column of the response matrix must be the base category.\n No need to specify `base.'") 
    base <- lev[p]
  } else { # standard Multinomial Probit model        
    Y <- as.factor(Y)
    lev <- levels(Y)
    if (!is.null(base))
      if (base %in% lev) {
        Y <- relevel(Y, ref = base)
        lev <- levels(Y)
      } else {
        stop(paste("Error: `base' does not exist in the response variable."))
      }
    base <- lev[1]
    counts <- table(Y)
    if (any(counts == 0)) {
      warning(paste("group(s)", paste(lev[counts == 0], collapse = " "), "are empty"))
      Y <- factor(Y, levels  = lev[counts > 0])
      lev <- lev[counts > 0]
    }
    p <- length(lev)
    Y <- as.matrix(unclass(Y)) - 1
    MoP <- FALSE
  }
  if(extra)
    return(list(Y=Y, MoP=MoP, lev=lev, p=p, base=base))
  else
    return(Y)
}
