checkData <-
  function(data){
    if (data.class(data) != "matrix" && data.class(data) != "data.frame")
      stop("Data are not matrix or data.frame")
    matrix.data <- as.matrix(data)
    if (any(is.na(matrix.data))) stop("Missing values are not allowed")
    if (mode(matrix.data)!="numeric") stop("Data must be numeric")
    if (any(matrix.data < 0)) stop("All scores should be nonnegative")
    if (any(matrix.data %% 1 !=0)) stop("All scores must be integers")
    matrix.data <- matrix.data - min(matrix.data)
    return(matrix.data)
  }