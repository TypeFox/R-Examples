"xinter"<-
  function(x, y, z, increasing = TRUE)
{
    
    if(increasing == FALSE) {
        x <- -1 * x
        x <- x[length(x):1]
        y <- y[length(y):1]
    }
    zz <- .Fortran("xinter",
                   as.double(x),
                   as.double(y),
                   length(x),
                   as.double(z),
                   result = double(1),
                   PACKAGE = "bootstrap")
    if(increasing == FALSE) {
        zz$result <- -1 * zz$result
    }
    return(zz$result)
}

