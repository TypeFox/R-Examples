#' Stratification and categorization of continuous data
#' 
#' Compute break points and marginal strata proportions, stratify and convert
#' continuous data (numeric) into categorical data (factor or integer).
#' 
#' @param x Vector, data frame or matrix; the continuous data to be processed.
#' 
#' @param breaks Vector or list; the lower and upper limits that should be used 
#' to break the continuous data into categories. See \sQuote{Details} for more 
#' information.
#' 
#' @param n Integer value; the number of strata that should be created.
#' 
#' @param type Character value; the type of strata, with options \code{"area"},
#' for equal-area, and \code{"range"}, for equal-range strata. Defaults to
#' \code{type = "area"}.
#' 
#' @param integer Logical value; should the categorical data be returned as
#' integers? Defaults to \code{integer = FALSE}.
#' 
#' @param prop Logical value; should the marginal strata proportions be 
#' returned? Defaults to \code{prop = FALSE}.
#' 
#' @details
#' Breaks must be a vector if \code{x} is a vector, but a list if \code{x} is a
#' data frame or matrix. Using a list allows breaking the data into a different
#' number of classes.
#' 
#' @return
#' A vector, data frame, or matrix, depending on the class of \code{x}.
#' 
#' @author Alessandro Samuel-Rosa \email{alessandrosamuelrosa@@gmail.com}
#' @seealso \code{\link[Hmisc]{cut2}}
#' @aliases cont2cat breakPoints stratify
#' @examples
#'
#'## Compute the break points of marginal strata
#'x <- data.frame(x = round(rnorm(10), 1), y = round(rlnorm(10), 1))
#'x <- breakPoints(x = x, n = 4, type = "area", prop = TRUE)
#'x
#'
#'## Convert continuous data into categorical data
#'# Matrix
#'x <- y <- c(1:10)
#'x <- cbind(x, y)
#'breaks <- list(c(1, 2, 4, 8, 10), c(1, 5, 10))
#'y <- cont2cat(x, breaks)
#'y
#'# Data frame
#'x <- y <- c(1:10)
#'x <- data.frame(x, y)
#'breaks <- list(c(1, 2, 4, 8, 10), c(1, 5, 10))
#'y <- cont2cat(x, breaks, integer = TRUE)
#'y
#'# Vector
#'x <- c(1:10)
#'breaks <- c(1, 2, 4, 8, 10)
#'y <- cont2cat(x, breaks, integer = TRUE)
#'y
#'
#'## Stratification
#'x <- data.frame(x = round(rlnorm(10), 1), y = round(rnorm(10), 1))
#'x <- stratify(x = x, n = 4, type = "area", integer = TRUE)
#'x
# FUNCTION - CONVERT CONTINUOUS DATA INTO CATEGORICAL DATA #####################
#' @export
#' @rdname cont2cat 
cont2cat <-
  function (x, breaks, integer = FALSE) {
    
    # Check if suggested packages are installed
    if (!requireNamespace("SpatialTools", quietly = TRUE)) {
      stop(paste("Package 'SpatialTools' needed for this function to work. ",
                 "Please install it.", sep = ""), call. = FALSE)
    }
    
    # Process input
    x_cl <- class(x)
    x <- as.data.frame(x)
    n_col <- ncol(x)
    if (n_col == 1) breaks <- list(breaks)
    
    for (i in 1:n_col) {
      
      # Check break points (adapted from Hmisc::cut2)
      r <- range(x[, i], na.rm = TRUE)
      if (r[1] < min(breaks[[i]])) breaks[[i]] <- c(r[1], breaks[[i]])
      if (r[2] > max(breaks[[i]])) breaks[[i]] <- c(breaks[[i]], r[2])
      
      # Cut data
      x[, i] <- cut(x = x[, i], breaks = breaks[[i]], include.lowest = TRUE,
                    right = FALSE)
    }
    
    # Process output
    if (integer) x <- sapply(x, as.integer)
    if (n_col == 1) x <- as.vector(x)
    if (x_cl == "matrix") x <- as.matrix(x)
    
    # Output
    return (x)
  }
# FUNCTION - COMPUTE BREAK POINTS AND STRATA PROPORTIONS #######################
#' @export
#' @rdname cont2cat 
breakPoints <-
  function (x, n, type = "area", prop = FALSE) {
    
    vec <- is.vector(x)
    x <- as.data.frame(x)
    n_col <- ncol(x)
    
    # Equal-area strata
    if (type == "area") {
      
      # Compute the break points using discrete sample quantiles
      probs <- seq(0, 1, length.out = n + 1)
      breaks <- lapply(x, stats::quantile, probs, na.rm = TRUE, type = 3)
      
    } else { # Equal-range strata
      
      # Compute the break points
      breaks <- lapply(1:n_col, function(i)
        seq(min(x[, i]), max(x[, i]), length.out = n + 1))
      
      # Find and replace by the closest population value
      d <- lapply(1:n_col, function(i)
        SpatialTools::dist2(matrix(breaks[[i]]), matrix(x[, i])))
      d <- lapply(1:n_col, function(i) apply(d[[i]], 1, which.min))
      breaks <- lapply(1:n_col, function(i) breaks[[i]] <- x[d[[i]], i])
    }
    
    # Keep only the unique break points
    breaks <- lapply(breaks, unique)
    
    # Compute the proportion of points per marginal strata
    if (prop) {
      count <- lapply(1:n_col, function (i)
        graphics::hist(x[, i], breaks[[i]], plot = FALSE)$counts
      )
      prop <- lapply(1:n_col, function (i) count[[i]] / sum(count[[i]]))
      names(prop) <- colnames(x)
      res <- list(breaks = breaks, prop = prop)
      
    } else {
      if (vec) {
        res <- unlist(breaks)
      } else {
        res <- breaks
      }
    }
    
    # Output
    return (res)
  }
# FUNCTION - MARGINAL STRATIFICATION ###########################################
#' @export
#' @rdname cont2cat 
stratify <-
  function (x, n, type = "area", integer = FALSE) {
    
    # Compute break points
    breaks <- breakPoints(x = x, n = n, type = type)
    
    # Convert continuous data into categorical data
    x <- cont2cat(x = x, breaks = breaks, integer = integer)
    
    # Output
    return (x)
  }
