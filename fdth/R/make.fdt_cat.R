make.fdt_cat <- function(f,
                         categories=NULL,
                         sort=TRUE,
                         decreasing=TRUE)
{
  if (is.null(categories))
    if (is.null(names(f)))
      categories  <- paste('V',
                           1:length(f),
                           sep='')
    else
      categories <- names(f)

    if (sort) {
      order.f <- order(f,
                       decreasing=decreasing)

      f <- f[order.f]

      categories <- categories[order.f]
    }

  f   <- as.vector(f)             # Frequency
  rf  <- f / sum(f)               # Relative freq
  rfp <- rf * 100                 # Relative freq, %
  cf  <- cumsum(f)                # Cumulative freq
  cfp <- cumsum(rfp)              # Cumulative freq, %

  res <- data.frame(categories,   # Make final table
                    f,
                    rf,
                    rfp,
                    cf,
                    cfp)

  names(res) <- c('Category',
                  'f',
                  'rf',
                  'rf(%)',
                  'cf',
                  'cf(%)')

  class(res) <- c('fdt_cat.default',
                  'data.frame')

  invisible(res)
}
