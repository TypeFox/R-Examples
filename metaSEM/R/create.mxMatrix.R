create.mxMatrix <- function(x, type=c("Full","Symm","Diag","Stand"), ncol=NA, nrow=NA, as.mxMatrix=TRUE, byrow=FALSE, ...) {
  my.mx <- function(y) {
    # suppress warnings
    ## warn <- options()$warn
    ## options(warn=-1)
    values <- suppressWarnings(as.numeric(y))  # They are NA for characters 
    free <- is.na(values)    # They are TRUE for parameters with labels 
    freePara1 <- y[free]     # Extract free parameters

    # check if there are any free parameters
    if (length(freePara1)>0) {
      freePara2 <- strsplit(freePara1, "*", fixed=TRUE)
      # replace NA with starting values 0.5 before "0.5*a"
      values[free] <- sapply(freePara2, function(y){ as.numeric(y[1]) })
      labels <- rep(NA, length(y))
      labels[free] <- sapply(freePara2, function(y){ y[2] })
      out <- list(values=values, free=free, labels=labels)
    } else {
      out <- list(values=values, free=free, labels=NA)
    }
    ## options(warn=warn)
    out
  }

  type <- match.arg(type)
  
  if (as.mxMatrix) {
    
    switch(type,
           Full = { if (length(x) != ncol*nrow)
                      stop("Length of \"x\" does not match the dimensions of \"ncol\" and \"nrow\".\n")
                    untitled1 <- matrix(x, ncol=ncol, nrow=nrow, byrow=byrow)
                    untitled1 <- as.mxMatrix(untitled1, ...)
                  },
           Symm = { no.var <- (sqrt(1 + 8 * length(x)) - 1)/2
                    if (abs(no.var - round(no.var)) > .Machine$double.eps^0.5)
                      stop("Length of \"x\" does not match the no. of elements for a symmetric matrix.\n")
                    my.x <- my.mx(x)
                    untitled1 <- mxMatrix(type="Symm", nrow=no.var, ncol=no.var, free=my.x$free,
                                          values=my.x$values, labels=my.x$labels, byrow=byrow, ...)   
                  },
           Diag = { no.var <- length(x)
                    my.x <- my.mx(x)
                    untitled1 <- mxMatrix(type="Diag", nrow=no.var, ncol=no.var, free=my.x$free,
                                          values=my.x$values, labels=my.x$labels, ...)                                         
                   },
           Stand = { no.var <- (sqrt(1 + 8 * length(x)) + 1)/2
                    if (abs(no.var - round(no.var)) > .Machine$double.eps^0.5)
                      stop("Length of \"x\" does not match the no. of elements for a standardized matrix.\n")
                    my.x <- my.mx(x)
                    untitled1 <- mxMatrix(type="Stand", nrow=no.var, ncol=no.var, free=my.x$free,
                                          values=my.x$values, labels=my.x$labels, byrow=byrow, ...)   
                  })
  } else {
        switch(type,
           Full = { if (length(x) != ncol*nrow)
                      stop("Length of \"x\" does not match the dimensions of \"ncol\" and \"nrow\".\n")
                    untitled1 <- matrix(x, ncol=ncol, nrow=nrow, byrow=byrow)                    
                  },
           Symm = { no.var <- (sqrt(1 + 8 * length(x)) - 1)/2
                    if (abs(no.var - round(no.var)) > .Machine$double.eps^0.5)
                      stop("Length of \"x\" does not match the no. of elements for a symmetric matrix.\n")
                    untitled1 <- matrix(0, ncol=no.var, nrow=no.var)
                    if (byrow) {
                        untitled1[upper.tri(untitled1, diag = TRUE)] <- x
                        untitled1[lower.tri(untitled1)] <- t(untitled1)[lower.tri(untitled1)]
                    } else {
                        untitled1[lower.tri(untitled1, diag = TRUE)] <- x
                        untitled1[upper.tri(untitled1)] <- t(untitled1)[upper.tri(untitled1)]
                    }
                  },
           Diag = { untitled1 <- Diag(x)                                         
                   },
           Stand = { no.var <- (sqrt(1 + 8 * length(x)) + 1)/2
                     if (abs(no.var - round(no.var)) > .Machine$double.eps^0.5)
                       stop("Length of \"x\" does not match the no. of elements for a standardized matrix.\n")
                     untitled1 <- Diag(rep(1,no.var))
                     if (byrow) {
                         untitled1[upper.tri(untitled1, diag = FALSE)] <- x
                         untitled1[lower.tri(untitled1)] <- t(untitled1)[lower.tri(untitled1)]
                     } else {
                         untitled1[lower.tri(untitled1, diag = FALSE)] <- x
                         untitled1[upper.tri(untitled1)] <- t(untitled1)[upper.tri(untitled1)]  
                     }
                  })    
  }
  untitled1
}

