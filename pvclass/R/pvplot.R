pvplot <- function(Y = NA, pv, alpha = NA) {
  par(mai=c(0.1, 0.1, 0.1, 0.1))
  pv <- sqrt(pv)
  n <- dim(pv)[1]
  L <- dim(pv)[2]
  if(is.na(alpha)) {
    alpha <- 0
  }
  if(is.null(dim(pv)[1])) {
    # only one observation
    n <- 1
    L <- length(pv)
    if(!is.na(Y[1])) {
      plot(c(-L/5, L, L, -L/5, -L/5), c(0, 0, n + n / 25, n + n / 25, 0),
           type = 'l', axes = FALSE, xlab = '', ylab = '')
      lines(x = c(-L / 5, L), y = c(n, n))
      for(j in seq_len(L)) {
        rect(j - .5 - pv[j] / 2, n - 1 + .5 - pv[j] / 2, j - .5 + pv[j] / 2,
             n - 1 + .5 + pv[j] / 2, col = 4 - 2 * (pv[j]^2 <= alpha), lty = 'blank')
      }
      for(j in seq_len(L)) {
        text(j - .5, n + n / 50, bquote(.(attr(Y, 'levels')[j])))
        lines(x = c(j - 1, j - 1), y = c(0, n + n / 25))
      }
      text(-L / 10, .5, bquote(.(substr(attr(Y, 'levels'), 1, 7)[Y]))) 
    } else {
      # Y = NA
      # plot number of observation
      par(mai = rep(0.1, 4))
      plot(c(-L / 5, L, L, -L / 5, -L / 5), c(0, 0, 1.25, 1.25, 0), type = 'l',
           axes = FALSE, xlab = '', ylab = '')
      lines(x = c(-L / 5, L), y = c(1, 1))
      text(-L / 10, .5, paste(1))
      for(j in seq_len(L)) {
        lines(x = c(j - 1, j - 1), y = c(0, 1.25))
        rect(j - .5 - pv[j] / 2, .5 - pv[j] / 2, j - .5 + pv[j] / 2,
             .5 + pv[j] / 2, col = 4 - 2 * (pv[j]^2 <= alpha), lty = 'blank')
      }
      for(j in seq_len(L)){
        text(j - .5, 1.125, bquote(.(attr(pv, 'names')[j])))
        lines(x = c(j - 1, j - 1), y = c(0, 1.25))
      }
    }
  } else {
    # more than one observation
    if(!is.na(Y[1])){ 
      # order Y and pv
      pv <- pv[order(Y), ]
      Y <- Y[order(Y)]
      # changes
      ch <- which(Y[seq_len(n-1)] != Y[2:n])
      
      plot(c(-L / 5, L, L, -L / 5, -L / 5), c(0, 0, n + n / 25, n + n / 25, 0),
           type = 'l', axes = FALSE, xlab = '', ylab = '')
      lines(x = c(-L / 5, L), y = c(n, n))
      for(j in seq_len(L)) {
        rect(j - .5 - pv[1, j] / 2, n - 1 + .5 - pv[1, j] / 2,
             j - .5 + pv[1, j] / 2, n - 1 + .5 + pv[1, j] / 2,
             col = 4 - 2 * (pv[1, j]^2 <= alpha), lty = 'blank')
      }
      for(i in seq_len(n)) {
        if(n < 101) {
          # draw horizontal lines
          if(any(ch == (i - 1))) {
            lines(x = c(-L / 5, L), y = c(n - i + 1, n - i + 1), lwd = 3)
          } else {
            lines(x = c(0, L), y = c(n - i + 1, n - i + 1))
          }
        } else {
          # only draw lines between different classes
          if(any(ch == (i-1))) {
            lines(x = c(-L / 5, L), y = c(n - i + 1, n - i + 1))
          }
        }
        for(j in seq_len(L)) {
          rect(j - .5 - pv[i, j] / 2, n - i + .5 - pv[i, j] / 2,
               j - .5 + pv[i, j] / 2, n - i + .5 + pv[i, j] / 2,
               col = 4 - 2 * (pv[i, j]^2 <= alpha), lty = 'blank')
        }
      }
      for(j in seq_len(L)) {
        text(j - .5, n + n / 50, bquote(.(attr(Y, 'levels')[j])))
        lines(x = c(j - 1, j - 1), y = c(0, n + n / 25))
      }
      ch <- c(0, ch, n)
      for(i in seq_len(length(ch)-1)) {
        text(-L / 10, n - ch[i] - (ch[i + 1] - ch[i]) / 2,
             bquote(.(substr(attr(Y, 'levels'), 1, 7)[i])))
      }
  } else {
    # Y = NA
    if(n > 25) {
      # do not plot number of observation
      par(mai = rep(0.1, 4))
      plot(c(0, L, L, 0, 0), c(0, 0, n + n / 25, n + n / 25, 0),
           type = 'l', axes = FALSE, xlab = '', ylab = '')
      for(i in seq_len(n)) {
        if(n < 101 | i == 1) {
          # draw horizontal lines
          lines(x = c(0, L), y = c(n - i + 1, n - i + 1))
        }
        for(j in seq_len(L)) {
          lines(x = c(j, j), y = c(0, n + n / 25))
          rect(j - .5 - pv[i, j] / 2, n - i + .5 - pv[i, j] / 2,
               j - .5 + pv[i, j] / 2, n - i + .5 + pv[i, j] / 2,
               col = 4 - 2 * (pv[i, j]^2 <= alpha), lty = 'blank')
        }
      }
      for(j in seq_len(L)) {
        text(j - .5, n + n / 50, bquote(.(attr(pv,'dimnames')[[2]][j])))
      }
    } else {
      # plot number of observation
      par(mai = rep(0.1, 4))
      plot(c(-L / 5, L, L, -L / 5, -L / 5), c(0, 0, n + 1, n + 1, 0),
           type = 'l', axes = FALSE, xlab = '', ylab = '')
      for(i in seq_len(n)) {
        lines(x = c(-L / 5, L), y = c(n - i + 1, n - i + 1))
        text(-L / 10, n - i + .5, paste(i))
        for(j in seq_len(L)) {
          lines(x = c(j - 1, j - 1), y = c(0, n + 1))
          rect(j - .5 - pv[i, j] / 2, n - i + .5 - pv[i, j] / 2,
               j - .5 + pv[i,j] / 2, n - i + .5 + pv[i, j] / 2,
               col = 4 - 2 * (pv[i, j]^2 <= alpha), lty = 'blank')
        }
      }
      if(is.null(attr(pv, 'dimnames'))) {
        attr(pv, 'dimnames')[[2]] <- seq_len(L)
      }
      for(j in seq_len(L)) {
        text(j - .5, n + .5, bquote(.(attr(pv, 'dimnames')[[2]][j])))
      }
    }
  }
}
}
