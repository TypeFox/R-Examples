make.fdt_cat.simple <- function(x,
                                sort,
                                decreasing)
{
  stopifnot(is.character(x) |
            is.factor(x))

  x <- as.factor(x)

  f <- as.vector(table(x))

  if (sort) {
    order.f <- order(f,
                     decreasing=decreasing)

    f <- f[order.f]

    category <- levels(x)[order.f]
  } else {
    category <- levels(x)
  }

  rf  <- f / sum(f)             # Relative freq
  rfp <- rf * 100               # Relative freq, %
  cf  <- cumsum(f)              # Cumulative freq
  cfp <- cumsum(rfp)            # Cumulative freq, %

  res <- data.frame(category,   # Make final table
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

  return(res)
}
