# check if 'y' and 'networks' are compatible and preprocess them if necessary
checkDataTypes <- function(y, networks = NULL, lag = 0) {
  
  # error message if both are NULL
  if (is.null(y) && is.null(networks)) {
    stop("No 'y' and 'networks' arguments were provided.")
  }
  
  # make sure that 'y' is a list
  stopmsg <- paste0("The data type of the 'y' argument is '", class(y), 
      "'. 'y' must be provided as a numeric vector, a list of numeric ", 
      "vectors, or a data frame with numeric columns.")
  if (is.null(y)) {
    # leave NULL as is
  } else if (class(y) == "list") {
    # OK
  } else if (class(y) == "data.frame") {
    # convert to a list of vectors
    time.steps <- ncol(y)
    y.copy <- y
    y <- list()
    for (i in 1:ncol(y.copy)) {
      y[[i]] <- y.copy[, i]
      names(y[[i]]) <- rownames(y.copy)
    }
    rm(y.copy)
  } else if (class(y) == "matrix") {
    stop(stopmsg)
  } else if (class(y) == "integer") {
    current.names <- names(y)
    y <- as.numeric(y)
    names(y) <- current.names
    y <- list(y)  # wrap y in list
  } else if (class(y) == "numeric") {
    y <- list(y)  # wrap y in list
  } else if (class(y) == "character") {
    y <- list(y)  # wrap y in list
  } else {
    tryCatch({  # try to convert to numeric and wrap in list
      current.names <- names(y)
      y <- as.numeric(y)
      names(y) <- current.names
      y <- list(y)
    }, error = function(e) {
      stop(stopmsg)
    })
  }
  
  # check validity of elements of 'y' list and convert if necessary
  if (!is.null(y)) {
    for (i in 1:length(y)) {
      cl <- class(y[[i]])
      stopmsg <- paste0("At t=", i, ", 'y' contains ", cl, " objects. 'y' ", 
          "must be provided as a numeric vector, a list of numeric vectors, ", 
          "or a data frame with numeric columns.")
      if (is.integer(y[[i]])) {
        current.names <- names(y[[i]])
        y[[i]] <- as.numeric(y[[i]])
        names(y[[i]]) <- current.names
      } else if (cl == "character") {
        stop(stopmsg)
      } else {
        tryCatch({  # try to convert to numeric
          current.names <- names(y[[i]])
          y[[i]] <- as.numeric(y[[i]])
          names(y[[i]]) <- current.names
        }, error = function(e) {
          stop(stopmsg)
        })
      }
    }
  }

  # make sure network matrices are in a list
  if (is.null(networks)) {
    # leave NULL as is
  } else if (class(networks) == "matrix") {
    # OK, but wrap in list
    networks <- list(networks)
  } else if (class(networks) == "network") {
    # cast as matrix object and wrap in list
    networks <- list(as.matrix(networks))
  } else if (class(networks) == "list") {
    # OK
  } else {
    # wrap in list
    networks <- list(networks)
  }
  
  # check validity of network matrices and their cells and convert if necessary
  if (!is.null(networks)) {
    for (i in 1:length(networks)) {
      if (class(networks[[i]]) != "matrix") {
        # try to convert element to matrix
        tryCatch({
          networks[[i]] <- as.matrix(networks[[i]])
        }, error = function(e) {
          stop(paste0("At t=", i, ", the object in the 'networks' list could ", 
              "not be converted to a matrix object."))
        })
      }
      if (storage.mode(networks[[i]]) != "numeric") {
        # try to set numeric storage mode
        tryCatch({
          storage.mode(networks[[i]]) <- "numeric"
        }, error = function(e) {
          stop(paste0("At t=", i, ", the matrix in the 'networks' list does ", 
              "not contain numeric values."))
        })
      }
    }
  }
  
  # check if length of 'y' and 'networks' is compatible and adjust if necessary
  if (!is.null(y) && !is.null(networks) && length(y) != length(networks)) {
    if (length(y) == 1 && length(networks) > 1) {
      stop(paste("'y' has only one time step, but the network has multiple ", 
          "time steps."))
    } else if (length(y) > 1 && length(networks) == 1) {
      for (i in length(y)) {
        networks[[i]] <- networks[[1]]  # inflate networks list
      }
    } else {
      stop(paste0("There should be the same number of elements in 'y' and ", 
          "'networks'. There are ", length(y), " elements in 'y' and ", 
          length(networks), " elements in 'networks'."))
    }
  }
  
  # check if dimensions and labels match (if present); add labels if not present
  if (!is.null(y) && !is.null(networks)) {
    for (i in 1:length(y)) {
      # compare dimensions; mutually adjust if necessary
      if (length(y[[i]]) != nrow(networks[[i]])) {
        if (is.null(rownames(networks[[i]]))) {
          stop(paste0("The dimensions of 'y' and 'networks' differ at t=", 
              i, ", and the network (matrix) does not contain row names."))
        } else if (is.null(names(y[[i]]))) {
          stop(paste0("The dimensions of 'y' and 'networks' differ at t=", 
              i, ", and the elements in 'y' are not named."))
        } else {  # try to adjust dimensions of 'y' and 'networks'
          warning(paste0("Dimensions of 'y' and 'networks' do not match at ", 
              "t=", i, ". Trying to adjust them mutually."))
          y[[i]] <- xergm.common::adjust(y[[i]], networks[[i]], add = FALSE)
          networks[[i]] <- xergm.common::adjust(networks[[i]], y[[i]], 
              add = FALSE)
        }
      }
      # complain if labels do not match
      if (!is.null(names(y[[i]])) && !is.null(rownames(networks[[i]])) && 
          !all(names(y[[i]]) == rownames(networks[[i]]))) {
        warning(paste("At t=", i, "the names of 'y' and the row names of", 
            "'networks' do not match."))
      }
      # fix labels; use integers if no labels present at all
      if (is.null(rownames(networks[[i]])) && !is.null(names(y[[i]]))) {
        rownames(networks[[i]]) <- names(y[[i]])
      } else if (is.null(names(y[[i]])) && !is.null(rownames(networks[[i]]))) {
        names(y[[i]]) <- rownames(networks[[i]])
      } else if (is.null(rownames(networks[[i]])) && is.null(names(y[[i]]))) {
        names(y[[i]]) <- 1:length(y[[i]])
        rownames(networks[[i]]) <- 1:nrow(networks[[i]])
      }
    }
  }
  
  # if only one of them present, add labels if necessary
  if (is.null(y)) {
    for (i in 1:length(networks)) {
      if (is.null(rownames(networks[[i]]))) {
        rownames(networks[[i]]) <- 1:nrow(networks[[i]])
        if (i > 1) {
          if (nrow(networks[[i]]) != nrow(networks[[i - 1]])) {
            stop(paste0("No row names in 'networks' at t=", i, ". Tried to ", 
                "create custom row names, but the dimensions differ from ", 
                "the previous time point."))
          }
          if (!all(rownames(networks[[i]]) == rownames(networks[[i - 1]]))) {
            stop(paste0("At t=", i, ", the row names of the network matrix ", 
                "do not match the row names of the previous time step."))
          }
        }
      }
    }
  }
  if (is.null(networks)) {
    for (i in 1:length(y)) {
      if (is.null(names(y[[i]]))) {
        names(y[[i]]) <- 1:length(y[[i]])
        if (i > 1) {
          if (length(y[[i]]) != length(y[[i - 1]])) {
            stop(paste0("No names in 'y' vector at t=", i, ". Tried to ", 
                "create custom names, but the length differs from ", 
                "the previous time point."))
          }
          if (!all(names(y[[i]]) == names(y[[i - 1]]))) {
            stop(paste0("At t=", i, ", the names of the 'y' vector ", 
                "do not match the names of the previous time step."))
          }
        }
      }
    }
  }
  
  # check 'lag'
  if (!is.null(y)) {
    n <- length(y)
  } else {
    n <- length(networks)
  }
  if (!is.numeric(lag)) {
    stop("The 'lag' argument must be numeric.")
  } else if (length(lag) > 1) {
    stop("The 'lag' argument must be of length 1.")
  } else if (lag < 0) {
    stop("The 'lag' argument must be >= 0.")
  } else if (n - lag < 1) {
    if (n == 1) {
      stop(paste("A lag of", lag, "was specified, but there is only", 
          "one time step."))
    } else {
      stop(paste("A lag of", lag, "was specified, but there are only", 
          n, "time steps."))
    }
  }
  
  # save everything in an object and return
  objects <- list()
  objects$y <- y
  objects$networks <- networks
  objects$time.steps <- n
  if (!is.null(y)) {
    objects$n <- lapply(y, length)
    objects$nodelabels <- unlist(lapply(y, names))
  } else {
    objects$n <- lapply(networks, nrow)
    objects$nodelabels <- unlist(lapply(networks, rownames))
  }
  names(objects$nodelabels) <- NULL
  
  return(objects)
}

