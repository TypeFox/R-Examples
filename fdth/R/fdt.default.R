fdt.default <- function (x,
                         k,
                         start,
                         end,
                         h,
                         breaks=c('Sturges', 'Scott', 'FD'),
                         right=FALSE, ...)
{
  # User defines nothing or not 'x' isn't numeric -> stop
  stopifnot(is.numeric(x))

  x <- na.omit(x)

  # User defines only 'x'
  if (missing(k) && 
      missing(start) && 
      missing(end) && 
      missing(h)) {

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
  else if (missing(start) && 
           missing(end) && 
           missing(h)) {
    stopifnot(length(k) >= 1)
    tmp   <- range(x)
    start <- tmp[1] - abs(tmp[1])/100
    end   <- tmp[2] + abs(tmp[2])/100
    R     <- end - start
    h     <- R/abs(k)
  }

  # User defines 'x', 'start' and 'end'
  else if (missing(k) && 
           missing(h)) {
    stopifnot(length(start) >= 1,
              length(end) >=1)
    tmp <- range(x)
    R   <- end - start
    k   <- sqrt(abs(R))
    if (k < 5) k = 5 # min value of k
    h   <- R/k
  }

  # User defines 'x', 'start', 'end' and 'h'
  else if (missing(k)) {
    stopifnot(length(start) >= 1,
              length(end) >= 1,
              length(h) >= 1)
  }

  else stop('Please, see the function sintaxe!')

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

  class(res) <- c('fdt.default',
                  'fdt',
                  'list')

  invisible(res)
} 
