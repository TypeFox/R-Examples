PlotClim <- function(exp_clim, obs_clim = NULL, toptitle = '', ytitle = '', 
                     monini = 1, freq = 12, limits = NULL, 
                     listexp = c('exp1', 'exp2', 'exp3'), 
                     listobs = c('obs1', 'obs2', 'obs3'), biglab = FALSE, leg = TRUE, 
                     fileout = 'output_plotclim.eps', sizetit = 1) {
  #
  #  Get some arguments
  # ~~~~~~~~~~~~~~~~~~~~
  #
  if (length(dim(exp_clim)) != 2 & length(dim(exp_clim)) != 3 ) {
    stop("2 or 3 dim needed : c(nexp, nltime) or c(nexp, nmemb, nltime)")
  }
  if (length(dim(exp_clim)) < 3) {
    exp_clim <- InsertDim(exp_clim, 2, 1)
  }
  nleadtime <- dim(exp_clim)[3]
  nexp <- dim(exp_clim)[1]
  if (is.null(obs_clim)) { 
    nobs <- 0
  } else { 
    nobs <- dim(obs_clim)[1]
    if (length(dim(obs_clim)) != 2 & length(dim(obs_clim)) != 3 ) {
      stop("2 or 3 dim needed : c(nobs, nltime) or c(nobs, nmemb, nltime)")
    } 
    if (length(dim(obs_clim)) < 3) {
      obs_clim <- InsertDim(obs_clim, 2, 1)
    }
    if (dim(obs_clim)[3] != nleadtime) {
      stop("obs and exp must have same number of ltimes")
    }
  }
  if (is.null(limits) == TRUE) {
    ll <- min(min(exp_clim, na.rm = TRUE), min(obs_clim, na.rm = TRUE), na.rm = TRUE)
    ul <- max(max(exp_clim, na.rm = TRUE), max(obs_clim, na.rm = TRUE), na.rm = TRUE)
    if (biglab) {
      ul <- ul + 0.3 * (ul - ll)
    } else {
      ul <- ul + 0.2 * (ul - ll)
    }
  } else {
    ll <- limits[1]
    ul <- limits[2]
  }
  lastyear <- (monini + (nleadtime - 1) * 12 / freq - 1) %/% 12
  lastmonth <- (monini + (nleadtime - 1) * 12 / freq - 1) %% 12 + 1
  #
  #  Define some plot parameters
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  if (biglab) {
    labind <- seq(1, nleadtime, max(nleadtime %/% 5, 1))
  } else {
    labind <- seq(1, nleadtime, max(nleadtime %/% 10, 1))
  }
  months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep",
              "Oct", "Nov", "Dec")
  labyear <- ((labind - 1) * 12 / freq + monini - 1) %/% 12
  labmonth <- months[((labind - 1) * 12 / freq + monini - 1) %% 12 + 1]
  for (jx in 1:length(labmonth)) {
    y2o3dig <- paste("0", as.character(labyear[jx]), sep = "") 
    labmonth[jx] <- paste(labmonth[jx], "\nYr ", substr(y2o3dig,
                          nchar(y2o3dig) - 1, nchar(y2o3dig)), sep = "")
  }
  empty_ts <- ts(start = c(0000, (monini - 1) %/% (12 / freq) + 1), 
                 end = c(lastyear, (lastmonth - 1) %/% (12 / freq) + 1), 
                 frequency = freq)
  empty <- array(dim = length(empty_ts))
  color <- c("red1", "dodgerblue1", "green1", "orange1", "lightblue1",
             "deeppink1", "mediumpurple1", "lightgoldenrod1", "olivedrab1", 
             "mediumorchid1")
  type <- c(1, 3, 2, 4)
  thickness <- c(1, 3, 1, 2)
  #
  #  Define plot layout
  # ~~~~~~~~~~~~~~~~~~~~
  # 
  postscript(fileout, width = 550, height = 300)
  if (biglab) {
    par(mai = c(1.25, 1.4, 0.5, 0), mgp = c(4, 2.5, 0))
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
       main = toptitle, cex.main = cexmain * sizetit, axes = FALSE)
  axis(1, at = labind, labels = labmonth)
  axis(2)
  box()
  #
  #  Loops on experimental and observational data
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  for (jexp in 1:nexp) {
    for (jmemb in 1:dim(exp_clim)[2]) {
      par(new = TRUE)
      plot(exp_clim[jexp, jmemb, ], type = "l", lty = 1, lwd = 2, 
           ylim = c(ll, ul), col = color[jexp], ylab = "", xlab = "", 
           axes = FALSE)
    }
  }
  if (nobs > 0) {
    for (jobs in 1:nobs) {
      for (jmemb in 1:dim(obs_clim)[2]) {
        par(new = TRUE)
        plot(obs_clim[jobs, jmemb, ], lty = type[jobs], lwd = thickness[jobs],
             type = "l", ylim = c(ll, ul), col = 1, ylab = "", xlab = "", 
             axes = FALSE)
      }
    }
    if (leg) {
      legend(1, ul, c(listexp[1:nexp], listobs[1:nobs]), 
             lty = c(array(1, dim = nexp), type[1:nobs]),
             lwd = c(array(2, dim = nexp), thickness[1:nobs]),
             col = c(color[1:nexp], array(1, dim = nobs)), cex = legsize)
    }
  } else {
    if (leg) {
      legend(1, ul, listexp[1:nexp], lty = 1, lwd = 2, col = color[1:nexp],
             cex = legsize)
    }
  }
  dev.off()
}
