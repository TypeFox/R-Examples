fdt_cat.matrix <- function (x,
                            sort=TRUE,
                            decreasing=TRUE, ...)
{
  stopifnot(is.matrix(x))

  res <- list()

  x <- na.omit(x)
 
  for (i in 1:ncol(x)) {
    m <- as.matrix(x[ ,i])

    fdt <- make.fdt_cat.simple(m,
                               sort,
                               decreasing)
    
    tmpres <- list(table=list(fdt)[[1]]) 

    res <- c(res,
             list(tmpres))
  }

  if (is.null(colnames(x)))
    nms <- paste('Column',
                 1:ncol(x),
                 sep=':')
  else
    nms <- colnames(x)

  names(res) <- nms

  class(res) <- c('fdt_cat.multiple',
                  'fdt_cat',
                  'list')

  invisible(res)
}
