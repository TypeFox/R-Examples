fdt_cat.default <- function (x,
                             sort=TRUE,
                             decreasing=TRUE, ...)
{
  x <- na.omit(x)

  res <- make.fdt_cat.simple(x,
                             sort,
                             decreasing)

  class(res) <- c('fdt_cat.default',
                  'fdt_cat',
                  'data.frame')

  invisible(res)
}  
