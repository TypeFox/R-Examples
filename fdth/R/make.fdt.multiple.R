make.fdt.multiple <- function (x,
                               k,
                               breaks=c('Sturges', 'Scott', 'FD'),
                               right)
{
  x <- na.omit(x)

  # User defines only x and/or 'breaks'
  if (missing(k)) {
    brk <- match.arg(breaks)

    switch (brk,
            Sturges = (k <- nclass.Sturges(x)),
            Scott   = (k <- nclass.scott(x)),
            FD      = (k <- nclass.FD(x)))

    tmp   <- range(x)
    start <- tmp[1] - abs(tmp[1])/100
    end   <- tmp[2] + abs(tmp[2])/100
    R     <- end - start
    h     <- R/k
  }

  # User defines 'x' and 'k'
  else {
    tmp   <- range(x)
    start <- tmp[1] - abs(tmp[1])/100
    end   <- tmp[2] + abs(tmp[2])/100
    R     <- end - start
    h     <- R/abs(k)
  }

  fdt <- make.fdt.simple(x,
                         start,
                         end,
                         h,
                         right)

  breaks <- c(start,
              end,
              h,
              ifelse (right,
                      1,
                      0))

  names(breaks) <- c('start',
                     'end',
                     'h',
                     'right')

  res <- list(table=fdt,
              breaks=breaks)

  return (res)
}

