fdt.matrix <- function (x,
                        k,
                        breaks=c('Sturges', 'Scott', 'FD'),
                        right=FALSE, ...)
{
  stopifnot(is.matrix(x))

  res <- list()

  for (i in 1:ncol(x)) {
    m <- as.matrix(x[ ,i])

    fdt <- make.fdt.multiple(m,
                             k,
                             breaks,
                             right)

    tmpres <- list(table=fdt[[1]],
                   breaks=fdt[[2]])

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

  class(res) <- c('fdt.multiple',
                  'fdt',
                  'list')

  invisible(res)
} 
