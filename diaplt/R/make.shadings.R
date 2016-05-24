# make.shadings
# http://code.google.com/p/cowares-excel-hello/source/browse/trunk/util_r/
#
# Copyright (C) 2013 Tomizono
# Fortitudinous, Free, Fair, http://cowares.nobody.jp
#                            http://paidforeveryone.wordpress.com/

# generate shading patterns of specified number
# output a list of density vector and angle vector
#
#  density=NA and angle=NA : default : no shadings, generates two NULLs
#  density=T or angle=T : shadings with automatic range
#  density=10 : use the number as fixed density and genrates auto angles
#  angle=30 : use the number as fixed angle and generates auto densities
#  density=c(3, 20) : generates densities between 3 and 20 lines per inch
#  angle=c(30, 90) : generates angles between 30 and 90 degree
#
# densities are shuffled by fixed shuffling variable

make.shadings <- function(n, density=NA, angle=NA, verbose=FALSE) {
  shadings <- list(density=NULL, angle=NULL)

  if(is.na(density) && is.na(angle)) return(shadings)

  label <- c('density', 'angle')
  start <- c(12, 10)
  end <- c(36, 160)
  shuffling <- c(TRUE, FALSE)
  par <- list(density, angle)

  for(i in 1:2) {
    x <- par[[i]]

    if(length(x) >= n) {
      shadings[[i]] <- x[1L:n]
    } else {
      rx <- if(is.numeric(x)) rep(x, 2)[1:2] else rep(NA, 2)
      if(is.na(rx[1])) rx[1] <- start[i]
      if(is.na(rx[2])) rx[2] <- end[i]

      sx <- seq(from=rx[1], to=rx[2], length.out=n)

      if(shuffling[i]) {
        # shuffle by halves
        # 1, k, 2, k+1, 3, k+2,,, n
        shuffle <- seq(from=1, to=n, by=2)
        shuffle <- c(shuffle, shuffle + 1)[1L:n]
      } else {
        shuffle <- 1L:n
      }
      shadings[[i]] <- sx[order(shuffle)]
    }
  }

  if(verbose) {
    cat('# shading enabled\n')
    print(shadings)
  }

  shadings
}

