mean.fdt <- function(x, ...)
{
  brk <- with(x,
              seq(breaks['start'],
                  breaks['end'],
                  breaks['h']))

  mids <- 0.5 * (brk[-1] + 
                 brk[-length(brk)])

  y <- with(x,
            table[, 2])

  res <- sum(y * mids) / sum(y)

  return(res)
}
