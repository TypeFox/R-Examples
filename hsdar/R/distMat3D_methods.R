setMethod("show", signature(object = "DistMat3D"), 
          function(object)
{
  mat <- as.matrix(object, lyr = 1)
  colnames(mat) <- paste("[,",1:ncol(mat),"]", sep = "")
  row.names(mat) <- paste("[",1:ncol(mat),",]", sep = "")
  if (ncol(mat) > 10)
  {
    print(mat[1:5,c(1:5, c((ncol(mat)-5):ncol(mat)))])
    cat(paste("    ... (", ncol(mat)-10, " rows/cols omitted)\n\n", sep = ""))
    print(mat[c((ncol(mat)-5):ncol(mat)),c(1:5, c((ncol(mat)-5):ncol(mat)))])
  } else {
    print(mat)
  }
  if (object@nlyr > 1)
  {
    cat("\n")
    if (object@nlyr > 2)
      cat(paste("    ... (", object@nlyr-2, " layers omitted)\n\n", sep = ""))
    mat <- as.matrix(object, lyr = object@nlyr)
    if (ncol(mat) > 10)
    {
      print(mat[1:5,c(1:5, c((ncol(mat)-5):ncol(mat)))])
      cat(paste("    ... (", ncol(mat)-10, " rows/cols omitted)\n\n", sep = ""))
      print(mat[c((ncol(mat)-5):ncol(mat)),c(1:5, c((ncol(mat)-5):ncol(mat)))])
    } else {
      print(mat)
    }
  }
}         
)

setMethod("ncol", signature(x = "DistMat3D"),
          function(x)
{
  return(x@ncol)
}
)

setMethod("nrow", signature(x = "DistMat3D"),
          function(x)
{
  return(x@ncol)
}
)

setMethod("dim", signature(x = "DistMat3D"),
          function(x)
{
  return(c(x@ncol, x@ncol, x@nlyr))
}
)

setMethod("as.matrix", signature(x = "DistMat3D"), 
          function(x, lyr = 1)
{
  val_vec <- c(((lyr-1) * (sum(1:x@ncol)-x@ncol) + 1):
               ((lyr) * (sum(1:x@ncol)-x@ncol)))
  mat <- matrix(NA, ncol = x@ncol, nrow = x@ncol) 
  mat[lower.tri(mat)] <- x@values[val_vec]
  return(mat)
}         
)

setMethod("as.array", signature(x = "DistMat3D"), 
          function(x)
{
  arr <- array(NA, dim = c(x@ncol, x@ncol, x@nlyr))
  for (i in 1:x@nlyr)
    arr[,,i] <- as.matrix(x, lyr = i)
  return(arr)
}         
)

setMethod("apply", signature(X = "DistMat3D"),
          function(X, MARGIN, FUN, ...)
{
  if (length(MARGIN) > 3)
    stop("Got unvalid MARGINs")
  if (length(MARGIN) == 3)
    stop("Got unvalid MARGINs")
  if (length(MARGIN) == 2)
  {
    MARGIN <- sort(MARGIN) 
    if (MARGIN[1] == 1 & MARGIN[2] == 2)
    {
      MARGIN  <- 1
    } else {
      stop("Got unvalid MARGINs")
    }
  } else {
    if (MARGIN == 2)
      MARGIN <- 1
  } 
    
  if (MARGIN == 1)
  {
    dummy_matrix <- matrix(data = NaN, nrow = dim(X)[1],
                           ncol = dim(X)[2])
    indices <- matrix(c(rep.int(1:dim(X)[1],
                                dim(X)[1])[row(dummy_matrix) >
                                col(dummy_matrix)],
                        rep.int(1:(dim(X)[1]-1),(dim(X)[1]-1):1)),
                      ncol = 2
                    )
    return(apply(indices, MARGIN = 1,
                 FUN = function(i, FUN, data, ...)
                 {
                   data <- data[i[1],i[2],]
                   return(FUN(data, ...))
                 },
                 FUN, X, ...)
           )      
  } else {
    indices <- matrix(1:X@nlyr, ncol = 1)
    return(apply(indices, MARGIN = 1,
                 FUN = function(i, FUN, data, ...)
                 {
                   data <- data[,,i]
                   data <- data[lower.tri(data)]
                   return(FUN(data, ...))
                 },
                 FUN, X, ...)
           )
  }
}
)

setMethod("[", signature(x = "DistMat3D"), 
          function(x, i, j, n)
{
  if (all(c(missing(i), missing(j), missing(n))))
    return(as.array(x))
  if (missing(i))
  {
    index_i <- c(1:x@ncol)
  } else {
    index_i <- i
  }
  if (missing(j))
  {
    index_j <- c(1:x@ncol)
  } else {
    index_j <- j
  }
 
  if (missing(n))
  {
    index_n <- c(1:x@nlyr)
  } else {
    index_n <- n
  }
  return(.doExtract(x, index_i, index_j, index_n))
}         
)

setReplaceMethod("[", signature(x = "DistMat3D"), 
          function(x, i, j, n, value)
{
  if (missing(i))
  {
    index_i <- c(1:x@ncol)
  } else {
    index_i <- i
  }
  if (missing(j))
  {
    index_j <- c(1:x@ncol)
  } else {
    index_j <- j
  }
 
  if (missing(n))
  {
    index_n <- c(1:x@nlyr)
  } else {
    index_n <- n
  }
  return(.doReplace(x, index_i, index_j, index_n, value))
}         
)



.doReplace <- function(x, i, j, n, value)
{
  res <- as.array(x)
  res[i, j, n] <- value
  return(distMat3D(res))
}

.doExtract <- function(x, i, j, n)
{
  res <- try(as.array(x)[i,j,n], silent = TRUE)
  if (!inherits(res, "try-error"))
    return(res)
    
  res <- as.matrix(x, lyr = n[1])
  res <- res[i,j]  
  if (length(n) > 1)
  {
    res_arr <- array(0, dim = c(length(i), length(j), length(n)))
    res_arr[,,1] <- res
    for (lyr in 2:length(n))
    {
      mat <- as.matrix(x, lyr = n[lyr])
      res_arr[,,lyr] <- mat[i,j]
    }
    return(res_arr)
  } else {
    return(res)
  }
}

# if (!isGeneric("distMat3D")) {
#   setGeneric("distMat3D", function(x, ...)
#   standardGeneric("distMat3D"))
# }

setMethod("distMat3D", signature(x = "array"), 
          function(x, lower_tri = TRUE)
{
  tri_fun <- if (lower_tri) lower.tri else upper.tri
  vals <- numeric(dim(x)[3] * (sum(1:dim(x)[1])-dim(x)[2]))
  for (lyr in 1:dim(x)[3])
  {
    val_vec <- c(((lyr-1) * (sum(1:dim(x)[1])-dim(x)[1]) + 1):
                ((lyr) * (sum(1:dim(x)[1])-dim(x)[1])))
    mat <- as.matrix(x[,,lyr])
    vals[val_vec] <- mat[tri_fun(mat)]
  }
  return(new("DistMat3D", values = vals, nlyr = dim(x)[3], ncol = dim(x)[1]))
}
)

setMethod("distMat3D", signature(x = "matrix"), 
          function(x, lower_tri = TRUE)
{
  tri_fun <- if (lower_tri) lower.tri else upper.tri
  vals <- x[tri_fun(x)]
  return(new("DistMat3D", values = vals, nlyr = 1, ncol = dim(x)[1]))
}
)

setMethod("distMat3D", signature(x = "numeric"),
          function(x, ncol, nlyr)
{
  return(new("DistMat3D", values = x, nlyr = nlyr, ncol = ncol))
}
)


setMethod("<", signature(e1='DistMat3D'),
    function(e1, e2){ 
    e1@values < e2
  }
)

setMethod("<=", signature(e1='DistMat3D'),
    function(e1, e2){ 
    e1@values <= e2
  }
)

setMethod(">", signature(e1='DistMat3D'),
    function(e1, e2){ 
    e1@values > e2
  }
)

setMethod(">=", signature(e1='DistMat3D'),
    function(e1, e2){ 
    e1@values >= e2
  }
)

setMethod("==", signature(e1='DistMat3D'),
    function(e1, e2){ 
    e1@values == e2
  }
)


setMethod("<", signature(e1='DistMat3D', e2='DistMat3D'),
    function(e1, e2){ 
    e1@values < e2@value
  }
)

setMethod("<=", signature(e1='DistMat3D', e2='DistMat3D'),
    function(e1, e2){ 
    e1@values <= e2@value
  }
)

setMethod(">", signature(e1='DistMat3D', e2='DistMat3D'),
    function(e1, e2){ 
    e1@values > e2@value
  }
)

setMethod(">=", signature(e1='DistMat3D', e2='DistMat3D'),
    function(e1, e2){ 
    e1@values >= e2@value
  }
)

setMethod("==", signature(e1='DistMat3D', e2='DistMat3D'),
    function(e1, e2){ 
    e1@values == e2@value
  }
)


setMethod("<", signature(e2='DistMat3D'),
    function(e1, e2){ 
    e1 < e2@values
  }
)

setMethod("<=", signature(e2='DistMat3D'),
    function(e1, e2){ 
    e1 <= e2@values
  }
)

setMethod(">", signature(e2='DistMat3D'),
    function(e1, e2){ 
    e1 > e2@values
  }
)

setMethod(">=", signature(e2='DistMat3D'),
    function(e1, e2){ 
    e1 >= e2@values
  }
)

setMethod("==", signature(e2='DistMat3D'),
    function(e1, e2){ 
    e1 == e2@values
  }
)
