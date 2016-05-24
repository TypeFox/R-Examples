fixNA <- function(x, y, spline = TRUE, verbose = FALSE) {
  
  #        Number of missing values "nNA"
  nNA <- sum(is.na(y))
  #test if too many NA values are present
  if (nNA/length(y) > 0.3) 
    warning("NA values constitutes more than 0.3 of the input data.\nApproximation may not be correct.")
  
  if (sum(is.na(tail(y))) > 2) 
    warning("More than 2 missing values in last 6 elements.\nApproximation may not be correct.")
  
  if (sum(is.na(head(y))) > 2) 
    warning("More than 2 missing values in first 6 elements.\nApproximation may not be correct.")
  
  # Test if at least one sequence of four missing elements in a data sequence is
  # present. Give a warning in such cases.
  na.seq <- y
  na.seq[is.na(na.seq) == FALSE] <- 0
  na.seq[is.na(na.seq) == TRUE] <- 1
        
  res.na.out <- sapply(1L:(length(na.seq) - 4), function(i) {
	ifelse(sum(na.seq[i:(i + 4)]) == 4, TRUE, FALSE)
     }
  )
    
  na.fail <- FALSE
  if (sum(res.na.out) >= 2) {
      na.fail <- TRUE 
  }
  
  if (na.fail) {
    warning("Sequence of more than 4 missing values in data.\nApproximation may not be correct.")
  }
  # Test if y contains somethings else than NAs
  #if(complete.cases(y) != is.na(y))
  #  warnings("Use ordinate contain non-number elements (e.g., strings)")
  
  
  
  # Indicate if information about the number of missing
  # values is needed
  if (verbose) 
    message(paste(nNA, "missing value(s) imputed.", sep = " "))
  
  # If NAs are present in the data set, substitie them by
  # linear approximation or by splines
  if ((nNA > 0) && (class(try(approx(x, y, n = length(x)), 
                              silent = TRUE))!="try-error")) {
    if (!spline) y[which(is.na(y))] <- approx(x, y, 
                                              n = length(x))$y[which(is.na(y))]
    if (spline) y[which(is.na(y))] <- spline(x, y, 
                                             n = length(y))$y[which(is.na(y))]
  } else {
    # If imputation fails use 0 instead
    y[which(is.na(y))] <- 0
  }
  
  #   if (length(which(is.na(y) == TRUE)) == 0)
  #     y <- y
  y
}

setGeneric("fixNA")

# setMethod("fixNA", signature(x = "numeric", y = "numeric"), 
#           function(x, y, spline = TRUE, verbose = FALSE) {
#             # Test if x and y have identical lengths.
#             if (is.null(y)) 
#               stop("Enter ordinate value")
#             
#             if (length(x) != length(y)) 
#               stop("Use abscissa and ordinate data with same number of elements")
#             
#             
#             fixNA(x, y, spline, verbose)
#           })

setMethod("fixNA", signature(x = "data.frame", y="missing"), 
          function(x, y, spline = TRUE, verbose = FALSE) { 
            if (ncol(x) != 2) 
              stop("'x' must have two columns.")
            
            fixNA(x[, 1], x[, 2], spline, verbose)
          })

setMethod("fixNA", signature(x = "matrix", y = "missing"), 
          function(x, y, spline = TRUE, verbose = FALSE) { 
            if (ncol(x) != 2) 
              stop("'x' must have two columns.")
            
            fixNA(x[, 1], x[, 2], spline, verbose)
          })


