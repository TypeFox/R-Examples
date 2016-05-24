PlotVsLTime <- function(var, toptitle = '', ytitle = '', monini = 1, freq = 12, 
               nticks = NULL, limits = NULL, 
               listexp = c('exp1', 'exp2', 'exp3'), 
               listobs = c('obs1', 'obs2', 'obs3'), biglab = FALSE, hlines = NULL, 
               leg = TRUE, siglev = FALSE, fileout = 'output_plotvsltime.eps', 
               sizetit = 1, show_conf = TRUE) {
  #
  #  Get some arguments
  # ~~~~~~~~~~~~~~~~~~~~
  #
  if (length(dim(var)) == 3) {
    var <- InsertDim(var, posdim = 2, lendim = 1)
  }
  nleadtime <- dim(var)[4]
  nexp <- dim(var)[1]
  nobs <- dim(var)[2]
  if (is.null(limits) == TRUE) {
    ll <- min(var, na.rm = TRUE)
    ul <- max(var, na.rm = TRUE)
    if (biglab) {
      ul <- ul + 0.4 * (ul - ll)
    } else {
      ul <- ul + 0.3 * (ul - ll)
    }
  } else {
    ll <- limits[1]
    ul <- limits[2]
  }
  lastyear <- (monini + (nleadtime - 1) * 12 / freq - 1) %/% 12
  lastmonth <- (monini + (nleadtime - 1) * 12 / freq - 1) %% 12 + 1
  empty_ts <- ts(start = c(0000, (monini - 1) %/% (12 / freq) + 1), 
                 end = c(lastyear, (lastmonth - 1) %/% (12 / freq) + 1),
                 frequency = freq)
  empty <- array(dim = length(empty_ts))
  #
  #  Define some plot parameters
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  if (is.null(nticks)) {
    if (biglab) {
      nticks <- 5
    } else {
      nticks <- 10
    }
  }
  labind <- seq(1, nleadtime, max(nleadtime %/% nticks, 1))
  months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep",
              "Oct", "Nov", "Dec")
  labyear <- ((labind - 1) * 12 / freq + monini - 1) %/% 12
  labmonth <- months[((labind - 1) * 12 / freq + monini -1 ) %% 12 + 1]
  for (jx in 1:length(labmonth)) {
    y2o3dig <- paste("0", as.character(labyear[jx]), sep = "")
    labmonth[jx] <- paste(labmonth[jx], "\nYr ", substr(y2o3dig, nchar(y2o3dig)
                    - 1, nchar(y2o3dig)), sep = "")
  }
  color <- c("red1", "dodgerblue1", "green1", "orange1", "lightblue1", 
             "deeppink1", "mediumpurple1", "lightgoldenrod1", "olivedrab1", 
             "mediumorchid1")
  type <- c(1, 3, 2, 4)
  thickness <- array(dim = c(4, 4))
  thickness[, 1] <- c(1, 2, 1, 1.5)
  thickness[, 2] <- c(8, 12, 8, 10)
  thickness[, 3] <- thickness[, 1]
  thickness[, 4] <- c(4, 6, 4, 5)
  if (siglev == TRUE) {
    lines <- c("n", "l", "n", "l")
  } else {
    lines <- c("l", "l", "l", "n")
  }
  #
  #  Define plot layout
  # ~~~~~~~~~~~~~~~~~~~~
  #
  postscript(fileout, width = 550, height = 300)
  if (biglab) {
    par(mai = c(1.25, 1.4, 0.5, 1), mgp = c(4, 2.5, 0))
    par(cex = 1.3, cex.lab = 2, cex.axis = 1.8)
    cexmain <- 2.2
    legsize <- 1.5
  } else {
    par(mai = c(1, 1.1, 0.5, 0), mgp = c(3, 1.8, 0))
    par(cex = 1.3, cex.lab = 1.5, cex.axis = 1.1)
    cexmain <- 1.5
    legsize <- 1
  }
  plot(empty, ylim = c(ll, ul), xlab = "Time (months)", ylab = ytitle, 
       main = toptitle, cex.main = cexmain*sizetit, axes = FALSE)
  axis(1, at = labind, labels = labmonth)
  axis(2)
  box()
  if (is.null(hlines) != TRUE) { 
    for (jy in 1:length(hlines)) {
      par(new = TRUE) 
      abline(h = hlines[jy])
    }
  }
  #
  #  Loop on experimental & observational data
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #                
  legendnames <- array(dim = nobs * nexp)
  legendthick <- array(dim = nobs * nexp)
  legendsty <- array(dim = nobs * nexp)
  legendcol <- array(dim = nobs * nexp)
  ind <- 1
  if (show_conf == TRUE) {
    start_line <- dim(var)[3]
    end_line <- 1
  } else {
    start_line <- 2
    end_line <- 2
  }
  for (jt in seq(start_line, end_line, -1)) {
    ind <- 1
    for (jexp in 1:nexp) {
      for (jobs in 1:nobs) {
        par(new = TRUE)
        plot(var[jexp, jobs, jt, ], type = lines[jt], ylim = c(ll, ul), 
             col = color[jexp], lty = type[jobs], lwd = thickness[jobs, jt], 
             ylab = "", xlab = "", axes = FALSE)
        legendnames[ind] <- paste(listexp[jexp], 'vs', listobs[jobs])
        legendthick[ind] <- thickness[jobs, 1] * 3
        legendsty[ind] <- type[jobs]
        legendcol[ind] <- color[jexp]
        ind <- ind + 1
      }
    }
  }
  if (leg) {
    if (nobs == 1) {
      legendnames <- listexp[1:nexp]
    }  
    legend(1, ul, legendnames, lty = legendsty, lwd = legendthick,
           col = legendcol, cex = legsize)
  }
  dev.off()
}
