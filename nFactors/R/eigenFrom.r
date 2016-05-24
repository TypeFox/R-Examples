eigenFrom <-
function(x) {
 classType <- class(x)
 
 res <- switch (classType,
  data.frame  = "data",
  matrix      = "correlation",
  numeric     = "eigenvalues",
  stop("Not a data.frame, a matrix, or a numeric vector")
  )
  
 switch (res,
  data        = if (dim(x)[2] <= 2) stop("At least 3 variables must be supplied"),
  correlation = if (dim(x)[2] <= 2) stop("At least 3 variables must be supplied"),
  eigenvalues = if (length(x) <= 2) stop("A vector of 3 eigenvalues or more must be supplied")
  )
  
 if (res == "correlation") if (any(x[lower.tri(x)] != t(x)[lower.tri(t(x))])) {
  stop("A correlation/covariance matrix must be symetric, empirical data must
        come from a data.frame, or eigenvalues must directly come from a vector.
        Verify the documentation about the eigenFrom function.")
  }
  
 invisible(res)
 }
  



  