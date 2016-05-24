preprocessX <- function(data.x) {

  #------------------------------------------------------------------#
  # Verify sufficient number of columns in datasets                  #
  #------------------------------------------------------------------#
  nc <- ncol(data.x)
  if( nc < 3L ) {
    stop("data.x must include {ID, time, measurement(s)}.", 
         call. = FALSE)
  }

  #------------------------------------------------------------------#
  # ensure that patient ids are integers                             #
  #------------------------------------------------------------------#
  if( !is.integer(data.x[,1L]) ) {
    data.x[,1L] <- as.integer(round(data.x[,1L],0))
    cat("Patient IDs in data.x were coerced to integer.\n")
  }

  intPresent <- FALSE
  if( nc > 2L ) {
    for( i in 3L:nc ) {
      if( isTRUE(all.equal(data.x[,i], 1.0)) ) {
        intPresent <- TRUE
        break
      }
    }
  }

  #------------------------------------------------------------------#
  # Remove any cases for which all covariates are NA                 #
  #------------------------------------------------------------------#
  rmRow <- apply(data.x, 1, function(x){all(is.na(x))})
  data.x <- data.x[!rmRow,]

  #------------------------------------------------------------------#
  # Determine total number of covariates.                            #
  # Note this includes an intercept                                  #
  #------------------------------------------------------------------#
  if( intPresent ){
    nCov <- ncol(data.x) - 2L
  } else {
    nCov <- ncol(data.x) - 2L + 1L
    #----------------------------------------------------------------#
    # Add intercept to data.x                                        #
    #----------------------------------------------------------------#
    nms <- colnames(data.x)
    if( nc > 2L ) {
      data.x <- cbind(data.x[,1L:2L], 1.0, data.x[,3L:nc])
      colnames(data.x) <- c(nms[1L:2L], "(Intercept)", nms[3L:nc])
    } else {
      data.x <- cbind(data.x[,1L:2L], 1.0)
      colnames(data.x) <- c(nms[1L:2L], "(Intercept)")
    }
  }

  #------------------------------------------------------------------#
  # Set missing cases to 0.0                                         #
  #------------------------------------------------------------------#
  tst <- is.na(data.x)
  data.x[tst] <- 0.0

  #------------------------------------------------------------------#
  # Determine if any times in data.x are < 0                         #
  #------------------------------------------------------------------#
  if( any(data.x[,2L] < {-1.5e-8}) ) {
    stop("Time is negative in data.x.", call. = FALSE)
  }

  return(data.x)

}
