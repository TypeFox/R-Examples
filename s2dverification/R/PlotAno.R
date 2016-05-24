PlotAno <- function(exp_ano, obs_ano = NULL, sdates, 
                    toptitle = c('', '', '', '', '', '', '', '', '', '', '', '',
                                 '', '', ''),
                    ytitle = c('', '', '', '', '', '', '', '', '', '', '', '', 
                               '', '', ''), 
                    limits = NULL, legends = NULL, freq = 12, biglab = FALSE, 
                    fill = TRUE, memb = TRUE, ensmean = TRUE, linezero = FALSE, 
                    points = FALSE, vlines = NULL, 
                    fileout = c('output1_plotano.eps', 
                    'output2_plotano.eps', 'output3_plotano.eps', 
                    'output4_plotano.eps', 'output5_plotano.eps'), 
                    sizetit = 1) {           
  #
  #  Get some arguments
  # ~~~~~~~~~~~~~~~~~~~~
  #
  if (length(dim(exp_ano)) != 4 ) {
    stop("4 dim needed : c(nexp/nobs, nmemb, nsdates, nltime)")
  }
  nexp <- dim(exp_ano)[1]
  nmemb <- dim(exp_ano)[2]
  nleadtime <- dim(exp_ano)[4]
  nsdates <- dim(exp_ano)[3]
  if (is.null(obs_ano) == FALSE) { 
    nobs <- dim(obs_ano)[1] 
    if (length(dim(obs_ano)) != 4 ) {
      stop("4 dim needed : c(nexp/nobs, nmemb, nsdates, nltime)") 
    }
    if (dim(obs_ano)[3] != nsdates | dim(obs_ano)[4] != nleadtime ) {
      stop("obs and exp must have same number of sdates & ltimes") 
    }
  } else {
    nobs <- 0
  }
  if (is.null(limits) == TRUE) {
    if (memb) {
      ll <- min(min(exp_ano, na.rm = TRUE), min(obs_ano, na.rm = TRUE), na.rm = TRUE)
      ul <- max(max(exp_ano, na.rm = TRUE), max(obs_ano, na.rm = TRUE), na.rm = TRUE)
    }
    else{
      ll <- min(min(Mean1Dim(exp_ano, 2), na.rm = TRUE), min(obs_ano, na.rm = TRUE), 
                na.rm = TRUE)
      ul <- max(max(Mean1Dim(exp_ano, 2), na.rm = TRUE), max(obs_ano, na.rm = TRUE),
                na.rm = TRUE)
    }
    if (nobs > 0) {
      if (biglab) {
        ul <- ul + 0.3 * (ul - ll)
      } else {
        ul <- ul + 0.2 * (ul - ll)
      }
    }
  } else {
    ll <- limits[1]
    ul <- limits[2]
  }
  #
  #  Define some plot parameters
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #
  yearinit <- as.integer(substr(sdates[1], 1, 4))
  moninit <- as.integer(substr(sdates[1], 5, 6))
  lastyear <- as.integer(substr(sdates[nsdates], 1, 4)) + (moninit + (
                         nleadtime - 1) * 12 / freq - 1) %/% 12
  lastmonth <- (moninit + (nleadtime - 1) * (12 / freq) - 1) %% 12 + 1
  empty_ts <- ts(start = c(yearinit, (moninit - 1) %/% (12 / freq) + 1), 
                 end = c(lastyear, (lastmonth - 1) %/% (12 / freq) + 1), 
                 frequency = freq)
  color <- c("red4", "orange4", "lightgoldenrod4", "olivedrab4", "green4",
             "lightblue4", "dodgerblue4", "mediumpurple4", "mediumorchid4",
             "deeppink4")
  color <- c(color, color, color, color, color, color, color, color, color, 
             color, color)
  colorblock <- c("red1", "orange1", "lightgoldenrod1", "olivedrab1", "green1",
                  "lightblue1", "dodgerblue1", "mediumpurple1", "mediumorchid1",
                   "deeppink1")
  colorblock <- c(colorblock, colorblock, colorblock, colorblock, colorblock,
                  colorblock, colorblock, colorblock, colorblock, colorblock)
  type <- c(1, 3, 2, 4)
  thickness <- c(1, 3, 2, 2)
  #
  #  Loop on the experiments : one plot for each
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # 
  for (jexp in 1:nexp) {
    #
    #  Define plot layout
    # ~~~~~~~~~~~~~~~~~~~~
    #
    postscript(fileout[jexp], width = 550, height = 300)
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
    plot(empty_ts, ylim = c(ll, ul), xlab = "Time (years)", ylab = ytitle[jexp],
          main = toptitle[jexp], cex.main = cexmain * sizetit)
    # 
    #  Plot experimental data + all observational datasets sdate by sdate
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #
    for (jdate in 1:nsdates) {
      year0 <- as.integer(substr(sdates[jdate], 1, 4))
      mon0 <- as.integer(substr(sdates[jdate], 5, 6))
      start <- (year0 - yearinit) * freq + 1
      end <- start + nleadtime - 1
      var <- array(dim = c(nmemb, length(empty_ts)))
      var[, start:end] <- exp_ano[jexp, , jdate, ]
      #
      #  Compute parameters for filling max-min over members 
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
      #
      if (fill) { 
        par(new = TRUE)
        bordup <- array(dim = nleadtime)
        borddown <- array(dim = nleadtime)
        for (jt in 1:nleadtime) {
          bordup[jt] <- max(exp_ano[jexp, , jdate, jt], na.rm = TRUE)
          borddown[jt] <- min(exp_ano[jexp, , jdate, jt], na.rm = TRUE)
        }
        tmp <- c(start:end)
        xout <- is.na(bordup + borddown) 
        tmp <- tmp[which(xout == FALSE)]
        xx <- c(tmp, rev(tmp))
        bordup <- bordup[which(xout == FALSE)]
        borddown <- borddown[which(xout == FALSE)]
        yy <- c(bordup, rev(borddown))
        #
        #  Plotting 
        # ~~~~~~~~~~
        # 
        if (jdate == 1) {
          matplot(t(var), type = "l", lty = 1, lwd = 1, ylim = c(ll, ul),
                  col = color[jdate], xlab = "", ylab = "", axes = FALSE)
        } 
        # Max-min member range
        polygon(xx, yy, col = colorblock[jdate], border = NA) 
      }
      if (ensmean) {  # Ensemble-mean
        par(new = TRUE)
        if (points) {
          plot(Mean1Dim(t(var), 2), type = "p", lty = 1, lwd = 4, 
               ylim = c(ll, ul), col = color[jdate], xlab = "", ylab = "",
               axes = FALSE)
        } else {
          plot(Mean1Dim(t(var), 2), type = "l", lty = 1, lwd = 4, 
               ylim = c(ll, ul), col = color[jdate], xlab = "", ylab = "",
               axes = FALSE)
        }
      }
      if (memb) {
        par(new = TRUE)  # All members
        if (points) { 
          matpoints(t(var), type = "p", lty = 1, lwd = 1, pch = 20, 
                    ylim = c(ll, ul), col = color[jdate], xlab = "", ylab = "",
                    axes = FALSE)
        } else {
          matplot(t(var), type = "l", lty = 1, lwd = 1, ylim = c(ll, ul),
                  col = color[jdate], xlab = "", ylab = "", axes = FALSE)
        }
      }
      if (nobs > 0) {   
        for (jobs in 1:nobs) { 
          for (jmemb in 1:dim(obs_ano)[2]) {
            var <- array(dim = length(empty_ts))
            var[start:end] <- obs_ano[jobs, jmemb, jdate, ]
            par(new = TRUE)  # Observational datasets
            if (points) {
              plot(var, type = "p", lty = 1, lwd = 4, pch = 20, 
                   ylim = c(ll, ul), col = 1, ylab = "", xlab = "", 
                   axes = FALSE)
            } else {
              plot(var, lty = type[jobs], lwd = thickness[jobs], type = "l",
                   ylim = c(ll, ul), col = 1, ylab = "", xlab = "", 
                   axes = FALSE)
            }
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
        legend(0, ul, legends[1:nobs], lty = 3, lwd = 10, col = 1, 
               cex = legsize)
      } else {
        legend(0, ul, legends[1:nobs], lty = type[1:nobs], 
               lwd = thickness[1:nobs], col = 1, cex = legsize)
      }
    }
    dev.off()
  }
}
