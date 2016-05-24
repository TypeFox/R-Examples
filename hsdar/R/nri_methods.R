setMethod("$", signature(x = "Nri"),
          function(x, name)
  {
    slot(x, name)
  }
)

setMethod("wavelength", signature(object = "Nri"),
          function(object)
{
  return(object@wavelength)
}
)

setMethod("show", signature(object = "Nri"),
          function(object)
{
  x <- object
  cat(paste("Data: nri, dimension: ", dim(x$nri)[1], ", ",
            dim(x$nri)[2], ", ", dim(x$nri)[3], "\n", sep=""))
  print(x$nri)
  cat(paste("      wavelength of length =",
            length(wavelength(x)),"\n"))
  cat(paste("      fwhm",
            if (length(x$fwhm)==1) "is constant for all wavelength"
              else "for each wavelength","\n"))
  if (length(x@multivariate) > 0)
  {
    .print.glmnri(x@multivariate)
  }
  .printUsagehistory(x)
  invisible(dim(x$nri))
}
)

setMethod("print", signature(x = "Nri"), 
          function(x)
{
  show(x)
})

setMethod("as.matrix", signature(x = "Nri"),
          function(x, ..., named_matrix = TRUE)
{
  mat <- matrix(x$nri, nrow = dim(x$nri)[3], byrow = TRUE, ...)
  if (named_matrix)
  {
    bnd_nam_ch <- eval(parse(text = as.character(dimnames(x$nri)[1])))
    bnd_nam <- as.vector(vapply(bnd_nam_ch, function(b1, b2) {
        paste(b2, "_", b1, sep = "")
      }, character(length = length(bnd_nam_ch)), bnd_nam_ch))
    colnames(mat) <- bnd_nam
    rownames(mat) <- eval(parse(text = as.character(dimnames(x$nri)[3])))
  }
  bnd_idx <- data.frame(b1 = rep.int(c(1:dim(x$nri)[1]), dim(x$nri)[1]),
                        b2 = c(sapply(c(1:dim(x$nri)[1]),
                                      function(x,n) rep.int(x,n), dim(x$nri)[1]))
                       )
  return(mat[, bnd_idx[,1] < bnd_idx[,2]])
}
)

setMethod("as.data.frame", signature(x = "Nri"),
          function(x, row.names = NULL, optional = FALSE, na.rm = FALSE, ...)
{
  .ConvertNri <- function(x, ...)
  {
    lyr <- as.matrix(x)
    lt <- lower.tri(lyr)
    data <- matrix(0, ncol = sum(lt), nrow = x@nlyr)
    data[1,] <- lyr[lt]
    if (x@nlyr > 1)
    {
      for (i in 2:x@nlyr)
      {
        lyr <- as.matrix(x, lyr = i)
        data[i,] <- lyr[lt]
      }
    }
    return(data)
  }
  bnd_nam_data <- x@dimnames
  bnd_nam_ch <- character()
  for (i in 1:(length(bnd_nam_data[[1]])-1))
    for (k in (i+1):length(bnd_nam_data[[2]]))
      bnd_nam_ch <- c(bnd_nam_ch, paste(bnd_nam_data[[2]][k], bnd_nam_data[[1]][i], sep = "_"))
  nri_data <- as.data.frame(.ConvertNri(x@nri, ...), row.names = NULL, optional = FALSE, ...)
  names(nri_data) <- bnd_nam_ch
  if (na.rm)
  {
    rem <- apply(as.matrix(nri_data), 2, function (x) all(is.finite(x)))
    nri_data <- nri_data[,rem]
  }
  return(nri_data)
})