PlotACC <- function(ACC, sdates, toptitle = "", sizetit = 1, ytitle = "", 
                    limits = NULL, legends = NULL, freq = 12, biglab = FALSE, 
                    fill = FALSE, linezero = FALSE, points = TRUE, vlines = NULL, 
                    fileout = "output_PlotACC.eps") {
  if (length(dim(ACC)) != 5 | dim(ACC)[5] != 4) {
    stop("5 dim needed : c(nexp, nobs, nsdates, nltime, 4)")
  }
  nexp <- dim(ACC)[1]
  nobs <- dim(ACC)[2]
  nleadtime <- dim(ACC)[4]
  nsdates <- dim(ACC)[3]
  if (is.null(limits) == TRUE) {
    ll <- min(ACC, na.rm = TRUE)
    ul <- max(ACC, na.rm = TRUE)
    if (biglab) {
      ul <- ul + 0.3 * (ul - ll)
    } else {
      ul <- ul + 0.2 * (ul - ll)
    }
  } else {
    ll <- limits[1]
    ul <- limits[2]
  }
  yearinit <- as.integer(substr(sdates[1], 1, 4))
  moninit <- as.integer(substr(sdates[1], 5, 6))
  lastyear <- as.integer(substr(sdates[nsdates], 1, 4)) + (moninit + (
                         nleadtime - 1) * 12 / freq - 1) %/% 12
  lastmonth <- (moninit + (nleadtime - 1) * (12 / freq) - 1) %% 12 + 1
  empty_ts <- ts(start = c(yearinit, (moninit - 1) %/% (12 / freq) + 1), 
                 end = c(lastyear, (lastmonth - 1) %/% (12 / freq) + 1), 
                 frequency = freq)
  color <- c("red4", "dodgerblue4", "lightgoldenrod4", "deeppink4", 
             "mediumpurple4", "green4", "orange4", "lightblue4", "mediumorchid4", 
             "olivedrab4")
  colorblock <- c("red1", "dodgerblue1", "lightgoldenrod1", "deeppink1", 
                  "mediumpurple1", "green1", "orange1", "lightblue1", 
                  "mediumorchid1", "olivedrab1")
  postscript(fileout, width = 550, height = 300)
  if (biglab) {
    par(mai = c(1, 1.1, 0.5, 0), mgp = c(2.8, 0.9, 0))
    par(cex = 1.3, cex.lab = 2, cex.axis = 1.8)
    cexmain <- 2.2
    legsize <- 1.5
  } else {
    par(mai = c(0.8, 0.8, 0.5, 0.1), mgp = c(2, 0.5, 0))
    par(cex = 1.3, cex.lab = 1.5, cex.axis = 1.1)
    cexmain <- 1.5
    legsize <- 1
  }
  plot(empty_ts, ylim = c(ll, ul), xlab = "Time (years)", ylab = ytitle, 
       main = toptitle, cex.main = cexmain * sizetit)
  for (jexp in 1:nexp) {
    for (jobs in 1:nobs) {
      numcol <- jobs + (jexp - 1) * nobs
      for (jdate in 1:nsdates) {
        year0 <- as.integer(substr(sdates[jdate], 1, 4))
        mon0 <- as.integer(substr(sdates[jdate], 5, 6))
        start <- (year0 - yearinit) * freq + 1
        end <- start + nleadtime - 1
        var <- array(dim = c(3, length(empty_ts)))
        var[, start:end] <- t(ACC[jexp, jobs, jdate, , 1:3])
        if (fill) {
          par(new = TRUE)
          bordup <- ACC[jexp, jobs, jdate, , 3]
          borddown <- ACC[jexp, jobs, jdate, , 1]
          tmp <- c(start:end)
          xout <- is.na(bordup + borddown)
          tmp <- tmp[which(xout == FALSE)]
          xx <- c(tmp, rev(tmp))
          bordup <- bordup[which(xout == FALSE)]
          borddown <- borddown[which(xout == FALSE)]
          yy <- c(bordup, rev(borddown))
          if (jdate == 1) {
            matplot(t(var), type = "l", lty = 1, lwd = 1, ylim = c(ll, ul), 
                    col = color[numcol], xlab = "", ylab = "", axes = FALSE)
          }
          polygon(xx, yy, col = colorblock[numcol], border = NA)
        }
        if (points) {
          par(new = TRUE)
          plot(var[2, ], type = "p", lty = 1, lwd = 6, ylim = c(ll, ul), 
               col = color[numcol], xlab = "", ylab = "", axes = FALSE,
               cex = 0.6)
          par(new = TRUE)
          plot(var[1, ], type = "p", pch = 6, lwd = 3, ylim = c(ll, ul), 
               col = color[numcol], xlab = "", ylab = "", axes = FALSE,
               cex = 0.6)
          par(new = TRUE)
          plot(var[3, ], type = "p", pch = 2, lwd = 3, ylim = c(ll, ul), 
               col = color[numcol], xlab = "", ylab = "", axes = FALSE,
               cex = 0.6)
          par(new = TRUE)
          for (jind in start:end) {
            lines(c(jind, jind), var[c(1, 3), jind], lwd = 1, 
                  ylim = c(ll, ul), col = color[numcol], xlab = "", 
                  ylab = "", axes = FALSE)
          }
        } else {
          par(new = TRUE)
          plot(var[2, ], type = "l", lty = 1, lwd = 4, ylim = c(ll, ul), 
               col = color[numcol], xlab = "", ylab = "", axes = FALSE)
          par(new = TRUE)
          plot(var[1, ], type = "l", lty = 1, lwd = 1, ylim = c(ll, ul), 
               col = color[numcol], xlab = "", ylab = "", axes = FALSE)
          par(new = TRUE)
          plot(var[3, ], type = "l", lty = 1, lwd = 1, ylim = c(ll, ul), 
               col = color[numcol], xlab = "", ylab = "", axes = FALSE)
        }
      }
    }
  }
  if (linezero) {
    abline(h = 0, col = "black")
  }
  if (is.null(vlines) == FALSE) {
    for (x in vlines) {
      abline(v = x, col = "black")
    }
  }
  if (is.null(legends) == FALSE) {
    if (points) {
      legend(0, ul, legends[1:(nobs * nexp)], lty = 3, lwd = 10, 
             col = color[1:(nobs * nexp)], cex = legsize)
    } else {
      legend(0, ul, legends[1:(nobs * nexp)], lty = 1, lwd = 4, 
             col = color[1:(nobs * nexp)], cex = legsize)
    }
  }
  dev.off()
}
