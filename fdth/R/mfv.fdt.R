mfv.fdt <- function(x, ...)
{
  getMFV <- function(pos)
  {
    fdt <- with(x,
                table)

    liMFV <- brk[pos]

    ## It is important if 'pos' is inside of the first class
    if (pos - 1 <= 0)
      D1 <- fdt[pos, 2]
    else
      D1 <- fdt[pos, 2] - fdt[(pos - 1), 2]

    nrows <- dim(fdt)[1]

    ## It is important if 'pos' is inside of the last class
    if (pos + 1 > nrows)
      D2 <- fdt[pos, 2]
    else
      D2 <- fdt[pos, 2] - fdt[(pos + 1), 2]

    MFV <- liMFV + (D1 / (D1 + D2)) * h
  }

  y <- with(x,
            table[, 2])

  brk <- with(x,
              seq(breaks['start'],
                  breaks['end'],
                  breaks['h']))

  h <- as.vector(with(x,
                      breaks['h']))

  posMFV <- grep(max(y),
                 y)

  res <- sapply(posMFV,
                getMFV)

  return(res)
}                        
