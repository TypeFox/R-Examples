setGeneric("summary")
setGeneric("plot")

#just for writing comfort, self explanatory
setClassUnion("numericOrNULL",c("numeric","NULL"))
setOldClass("summary.lm")
setOldClass("htest")

#amptest class, amptester function -------------------------------
setClass("amptest", contains = "numeric", representation(.Data = "numeric", 
                                                         decisions = "logical",
                                                         noiselevel = "numeric",
                                                         background = "numericOrNULL",
                                                         polygon = "numeric",
                                                         slope.ratio = "numeric"))

setMethod("show", signature(object = "amptest"), function(object) {
  print(slot(object, ".Data"))
})

setMethod("summary", signature(object = "amptest"), function(object, print = TRUE) {
  #print(slot(object, ".Data")) I think it's too much and repeats show method
  if(print) {
    cat(paste0("\nAmplification significance (threshold test): ", 
               slot(object, "decisions")[["tht.dec"]]))
    cat(paste0("\nAmplification significance (signal level test): ", 
               slot(object, "decisions")[["slt.dec"]]))
    cat(paste0("\nAmplification significance (resids growth test): ", 
               slot(object, "decisions")[["rgt.dec"]]))
    cat(paste0("\nNoise detected: ", slot(object, "decisions")[["shap.noisy"]]))
    cat(paste0("\nNoise level: ", slot(object, "noiselevel")))
    cat(paste0("\nLinearity: ", slot(object, "decisions")[["lrt.test"]]))
    cat(paste0("\nPolygon: ", slot(object, "polygon")))
    bcg <- slot(object, "background")
    if (is.null(bcg)) {
      cat(paste0("\nBackground: not defined")) 
    } else {
      bcg <- paste0(bcg, collapse = ", ")
      cat(paste0("\nBackground: (", bcg, ")"))
    }
  }
  invisible(c(tht.dec = slot(object, "decisions")[["tht.dec"]],
              slt.dec = slot(object, "decisions")[["slt.dec"]],
              rgt.dec = slot(object, "decisions")[["rgt.dec"]],
              shap.noisy = slot(object, "decisions")[["shap.noisy"]],
              noiselevel = slot(object, "noiselevel"),
              lrt.test = slot(object, "decisions")[["lrt.test"]], 
              polygon = slot(object, "polygon")))
})


setMethod("plot", signature(x = "amptest"), 
          function(x, abscissa = 1L:length(slot(x, ".Data")), ...) {
            y.pos <- slot(x, ".Data")
            nh <- trunc(length(abscissa) * 0.20)
            nt <- trunc(length(abscissa) * 0.15)
            
            y.pos.head <- head(y.pos, n = nh)
            y.pos.tail <- tail(y.pos, n = nt)
            
            lb.pos <- median(y.pos.head) + 2 * mad(y.pos.head)
            ub.pos <- median(y.pos.tail) - 2 * mad(y.pos.tail)
            
            res.shapiro.pos <- shapiro.test(y.pos)
            
            res.wt.pos <- wilcox.test(head(y.pos, n = nh), tail(y.pos, n = nt), 
                                      alternative = "less")
            
            layout(matrix(c(1,1,1,2,3,4), ncol = 3, byrow = TRUE), respect = TRUE)
            
            y.lims <- range(y.pos)
            y.lims[1] <- min(y.lims[1], lb.pos)
            y.lims[2] <- max(y.lims[1], ub.pos)
            
            y.lims[1] <- y.lims[1] - 0.4*(lb.pos + abs(ub.pos))
            y.lims[2] <- y.lims[2] + 0.4*(lb.pos + abs(ub.pos))
            plot(abscissa, y.pos, xlim = c(abscissa[1], 
                                           abscissa[length(abscissa)] * 1.3), 
                 ylim = y.lims, xlab = "Cycle", ylab = "RFU", 
                 main = "Input data", type = "b", pch = 19)
            
            abline(v = c(nh, length(abscissa) - nt), lty = 3)
            
            abline(h = lb.pos, lty = 2, col = "red")
            text(abscissa[9], lb.pos * 1.2, "Noise\nmedian + 2 * mad", col = "red", 
                 cex = 1.3)
            
            abline(h = ub.pos, lty = 2, col = "green")
            text(abscissa[length(abscissa)] * 1.1, ub.pos * 0.90, 
                 "Signal\nmedian - 2 * mad", col = "green", cex = 1.3)
            
            arrows(4.5, 12.5, 42.5, 12.5, length = 0.1, angle = 90, code = 3)
            text(25, 14.5, paste("W = ", format(res.wt.pos[["statistic"]], digits = 4), 
                                 "\np-value = ", 
                                 format(res.wt.pos[["p.value"]], digits = 6)))
            text(5,5, paste("Fold change: \n", round(ub.pos/lb.pos, 2)))
            
            res.pos <- rtg.test(y.pos)
            plot(1L:ncol(res.pos), res.pos[nrow(res.pos), ], xaxt='n', 
                 xlab = "Cycle interval", ylab = "", ylim = c(0, 1), main = "RGt", 
                 pch = 19)
            mtext("Correlation coeffcient", 2, 3)
            mtext("(studentized residuals and RFU)", 2, 2)
            nice.labs <- sapply(1L:ncol(res.pos), function(i) 
              paste0(res.pos[c(1, nrow(res.pos) - 1), i], collapse = "-"))
            axis(1, 1L:ncol(res.pos), labels = nice.labs)
            abline(h = 0.8, lty = "66")
            
            qqnorm(y.pos, pch = 19, main = paste("W = ", 
                                                 format(res.shapiro.pos[["statistic"]], 
                                                        digits = 6), "\np-value = ", 
                                                 format(res.shapiro.pos[["p.value"]], 
                                                        digits = 6)))
            qqline(y.pos, col = "orange", lwd = 2)
            
            plot(RGt(y.pos), xlab = "Cycle", ylab = expression(R^2), main = "LRt", 
                 pch = 19, type = "b")
            abline(h = 0.8, col = "black", lty = 2)
          })

###


RGt <- function(y) {
  ws <- ceiling((15 * length(y)) / 100)
  if (ws < 5) 
    ws <- 5
  if (ws > 15) 
    ws <- 15
  y.tmp <- na.omit(y[-c(1:5)])
  x <- 1:length(y.tmp)
  suppressWarnings(
    res.reg <- sapply(1L:(length(y.tmp)), function (i)  {
      round(summary(lm(y.tmp[i:c(i + ws)] ~ x[i:c(i + ws)]))[["r.squared"]], 4)
    }
    )
  )
  
  # Binarize R^2 values. Everything larger than 0.8 is positve
  res.LRt <- res.reg
  # Define the limits for the R^2 test
  res.LRt[res.LRt < 0.8] <- 0
  res.LRt[res.LRt >= 0.8] <- 1
  # Seek for a sequence of at least six positve values (R^2 >= 0.8)
  # The first five measure points of the amplification curve are skipped
  # because most technologies and probe technologies tend to overshot
  # in the start (background) region.
  res.out <- sapply(5L:(length(res.LRt) - 6), function(i) {
    ifelse(sum(res.LRt[i:(i + 4)]) == 5, TRUE, FALSE)
  }
  )
  cbind(1L:(length(y.tmp)), res.reg)
}


rtg.test <- function(y) {
  nh <- trunc(length(y) * 0.2)
  if (nh < 5) nh <- 5
  rgts <-sapply(0L:round(length(y)/8, 0), function (j) {
    cyc <- 1:nh + j
    reg <- lm(y[cyc] ~ cyc)
    c(cyc, cor(rstudent(reg), y[cyc]))
  })
  rgts
}


#der class ----------------------
setClass("der", contains = "matrix", representation(.Data = "matrix", 
                                                    method = "character"))


setMethod("summary", signature(object = "der"), function(object, 
                                                         digits = getOption("digits") - 3, 
                                                         print = TRUE) {
  data <- slot(object, ".Data")
  FDM <- data[data[, "d1y"] == max(data[, "d1y"]), "x"] 
  SDM <- data[data[, "d2y"] == max(data[, "d2y"]), "x"]
  SDm <- data[data[, "d2y"] == min(data[, "d2y"]), "x"]
  SDC <- sqrt(SDM * SDm)
  if (print) {
    cat(paste0("Smoothing method: ", slot(object, "method")))
    cat(paste0("\nFirst derivative maximum: ", format(FDM, digits = digits)))
    cat(paste0("\nSecond derivative maximum: ", format(SDM, digits = digits)))
    cat(paste0("\nSecond derivative minimum: ", format(SDm, digits = digits)))
    cat(paste0("\nSecond derivative center: ", format(SDC, digits = digits)))
  }
  res <- c(FDM, SDM, SDm, SDC)
  names(res) <- c("FDM", "SDM", "SDm", "SDC")
  invisible(res)
})

setMethod("show", signature(object = "der"), function(object) {
  print(slot(object, ".Data"))
})

setMethod("plot", signature(x = "der"), function(x, what = 1:3, add = FALSE,
                                                 legend = TRUE,
                                                 plot.colors = c("black", "red", "blue"), 
                                                 ...) {
  if (!all(what %in% 1:3)) 
    stop("'what' must contain values from set: 1, 2, 3.")
  
  if (length(unique(what)) != length(what)) 
    stop("'what' must contain unique values.")
  
  if (length(plot.colors) != 3) 
    stop("'plot.colors' must contain three colors.")
  #smallest and biggest fluorescence values
  ylims <- range(sapply(what + 1, function(i) {
    x[, i]
  }))
  
  if(!add) {
    plot(x = range(x[, 1]),  y = ylims, cex = 0, ...)
  }
  


  for (i in what)
    points(x[, c(1, i + 1)], col = plot.colors[i], pch = 20, type = "b")
  
  
  if (legend)
    legend("topleft", c("Raw data", "First derivative", "Second derivative")[what], 
           pch = rep(20,3), col = plot.colors)
})


#bg class, bg.max function -------------------------------
setClass("bg", contains = "matrix", representation(.Data = "matrix", 
                                                   bg.start = "numeric",
                                                   bg.stop = "numeric",
                                                   bg.corr = "numeric",
                                                   fluo = "numeric",
                                                   amp.stop = "numeric"))


setMethod("summary", signature(object = "bg"), function(object, print = TRUE) {
  if (print) {
    cat(paste0("Background start: ", slot(object, "bg.start")))
    cat(paste0("\nBackground stop: ", slot(object, "bg.stop")))
    cat(paste0("\nBackground correction: ", slot(object, "bg.corr")))
    cat(paste0("\nEnd of the amplification reaction: ", slot(object, "amp.stop")))
    cat(paste0("\nFluorescence at the end of the amplification reaction: ", 
               format(slot(object, "fluo"), options("digits")[["digits"]])))
  }
  invisible(c(bg.start = slot(object, "bg.start"), 
              bg.stop = slot(object, "bg.stop"),
              bg.corr = slot(object, "bg.corr"),
              amp.stop = slot(object, "amp.stop"),
              fluo = slot(object, "fluo")))
})

setMethod("show", signature(object = "bg"), function(object) {
  print(slot(object, ".Data"))
})

setMethod("plot", signature(x = "bg"), function(x, what = 1:3, add = FALSE, 
                                                indicators = TRUE, 
                                                legend = TRUE, stan.labs = TRUE, 
                                                plot.colors = c("black", "red", "blue"), 
                                                ...) {
  if (stan.labs) {
    plot(new("der", .Data = x, method = "supsmu"), what = what, 
         add = add, legend = FALSE, plot.colors = plot.colors, xlab = "Cycle", 
         ylab = "Fluorescence", ...)
  } else {
    plot(new("der", .Data = x, method = "supsmu"), what = what, 
         add = add, legend = FALSE, plot.colors = plot.colors, ...)
  }
  
  if (indicators) {
    abline(v = slot(x, "bg.start"))
    text(slot(x, "bg.start"), 0.2, "Background start", pos = 4)
    abline(v = slot(x, "bg.stop"), col = "blue")
    text(slot(x, "bg.stop"), 0.25, "Background stop", pos = 4, col = "blue")
    abline(v = slot(x, "amp.stop"), col = "green")
    text(slot(x, "amp.stop"), 0.3, "Plateau transition", pos = 4, col = "green")
  }
  if (legend)
    legend("topleft", c("Raw data", "First derivative", "Second derivative")[what], 
           pch = rep(20,3), col = plot.colors)
})



#refMFI class, MFIaggr function -------------------------------
setClass("refMFI", contains = "matrix", representation(.Data = "matrix", 
                                                       density = "density", 
                                                       qqnorm.data = "data.frame",
                                                       stats = "numeric"))

setMethod("show", signature(object = "refMFI"), function(object) {
  print(slot(object, ".Data"))
})

setMethod("qqnorm", signature(y = "refMFI"), function(y, main = "Normal Q-Q Plot",
                                                      xlab = "Theoretical Quantiles",
                                                      ylab = "Sample Quantiles",
                                                      plot.it = TRUE, datax = FALSE, ...) {
  qqnorm(unlist(slot(y, "qqnorm.data")), main = main, xlab = xlab,
         ylab = ylab, plot.it = plot.it, datax = datax, ...)
})

setMethod("qqline", signature(y = "refMFI"), function(y, datax = FALSE, 
                                                      distribution = qnorm,
                                                      probs = c(0.25, 0.75), qtype = 7, ...) {
  qqline(unlist(slot(y, "qqnorm.data")), datax = datax, distribution = distribution,
         probs = probs, qtype = qtype, ...)
})


setMethod("summary", signature(object = "refMFI"), 
          function(object, digits = getOption("digits") - 3, print = TRUE) {
            stats <- slot(object, "stats")
            nice.names <- c("Mean", "Median", "Standard deviation", 
                            "Median Absolute Deviation", "Interquartile Range", 
                            "Medcouple", "Skewness", 
                            "SNR", "VRM", "Number of NAs", "Intercept", "Slope", 
                            "R squared", "Breusch-Pagan Test p-value"
            )
            if (print) {
              for(i in 1L:length(stats))
                cat(paste0(nice.names[i], ": ", 
                           format(stats[i], digits = digits), "\n"))
            }
            invisible(stats)
          })


setMethod("plot", signature(x = "refMFI"), function(x, CV = FALSE, type = "p", 
                                                    pch = 19, length = 0.05, 
                                                    col = "black") {
  res <- slot(x, ".Data")
  res.dens <- slot(x, "density")
  llul <- rownames(slot(x, "qqnorm.data"))
  stats <- slot(x, "stats")
  ncol.y <- ncol(slot(x, "qqnorm.data"))
  
  default.par <- par(no.readonly = TRUE)
  
  #Plot the Coefficient of Variance
  layout(matrix(c(1,2,1,3), 2, 2, byrow = TRUE))
  
  error.plot.text <- paste0("ROI samples: ", format(ncol.y, nsmall = 3), "\n",
                       "ROI mean : SD\n", format(stats[1], nsmall = 3), " \u00b1 ", 
                       format(stats[3], nsmall = 2), "\n",
                       "ROI median : MAD\n", format(stats[2], nsmall = 3), " \u00b1 ", 
                       format(stats[4], nsmall = 2))
  
  density.plot.text <- paste("ROI cycles: ", 
                             llul[1], " to ", llul[2], 
                             "\n", "bw ", 
                             round(res.dens[["bw"]], 3), 
                             "; ", "N = ", res.dens[["n"]])
  
  
  par(mar=c(7.8, 4.1, 4.1, 2.1))
  if (CV) {
    plot(res[, 1], res[, 4], xlab = "", ylab = "CV", 
         type = type, pch = pch, col = col,
         main = "Error plot")
    
  } else {
    plot(res[, 1], res[, 2], ylim = c(min(res[, 2] - res[, 3]), 
                                      max(res[, 2] + res[, 3])), 
         xlab = "", ylab = "MFI", 
         type = type, pch = pch, col = col,
         main = "Error plot")
    
    #Plot the location with error bars.
    arrows(res[, 1], res[, 2] + res[, 3], res[, 1], 
           res[, 2] - res[, 3], angle = 90, code = 3, 
           length = length, col = col)
    deviation.measure <- strsplit(colnames(res)[3], "(", fixed = TRUE)[[1]][2]
    deviation.measure <- substr(deviation.measure, 1, nchar(deviation.measure) - 1)
    mtext(paste0("Deviation: ", deviation.measure), 4, cex = 0.75)
  }
  
  # Add a range for the ROI
  abline(v = llul, col = "lightgrey", lwd = 1.25)
  
  
  mtext(error.plot.text, side = 1, line = 6.8, cex = 0.75)
  mtext("Cycle", side = 1, line = 2, cex = 0.8)

  par(mar=c(4.1, 4.1, 4.1, 2.1))
  
  # "Calculate" the Quantile-Quantile plots and density plots
  # and plot the results
  plot(res.dens, xlab = "", main = "Density", col = col)
  mtext("RFU", side = 1, line = 2, cex = 0.8)
  mtext(density.plot.text, side = 1, line = 4.85, cex = 0.75)
  
  # Analysis of the quantiles
  qqnorm(x, xlab = "", pch = pch, col = col)
  mtext("Theoretical Quantiles", side = 1, line = 2, cex = 0.8)
  mtext(paste0("\nBreusch-Pagan Test p-value: ", format(stats["heter.p"], digits = 4)),
        side = 1, line = 3, cex = 0.75)
  qqline(x, col = col)
  
  par(default.par)
})

#combo plot
setMethod("plot", signature(x = "refMFI", y = "refMFI"), function(x, y, CV = FALSE, type = "p", 
                                                    pch = 19, length = 0.05, 
                                                    col = "black") {
  res <- list(slot(x, ".Data"), slot(y, ".Data"))
  res.dens <- list(slot(x, "density"), slot(y, "density"))
  llul <- list(rownames(slot(x, "qqnorm.data")), rownames(slot(y, "qqnorm.data")))
  stats <- list(slot(x, "stats"), slot(y, "stats"))
  ncol.y <- list(ncol(slot(x, "qqnorm.data")), ncol(slot(y, "qqnorm.data")))
  
  default.par <- par(no.readonly = TRUE)
  
  #Plot the Coefficient of Variance
  layout(matrix(c(1,2,1,3), 2, 2, byrow = TRUE))
  
  error.plot.text <- paste0(
    #"samples: ", format(ncol.y[[1]], nsmall = 3), "; ",
    #format(ncol.y[[2]], nsmall = 3), "\n",
    "ROI mean : SD:\n A ", format(stats[[1]][1], nsmall = 3), " \u00b1 ", 
    format(stats[[1]][3], nsmall = 2), "\n B ",
    format(stats[[2]][1], nsmall = 3), " \u00b1 ", 
    format(stats[[2]][3], nsmall = 2), "\n",
    "ROI median : MAD:\n A ", format(stats[[1]][2], nsmall = 3), " \u00b1 ", 
    format(stats[[1]][4], nsmall = 2), "\n B ",
    format(stats[[2]][2], nsmall = 3), " \u00b1 ", 
    format(stats[[2]][4], nsmall = 2))
  
  density.plot.text <- paste("ROI cycle ", 
                             llul[[1]][1], " to ", llul[[1]][2], "; ",
                             llul[[2]][1], " to ", llul[[2]][2],
                             "\n", "bw ", 
                             round(res.dens[[1]][["bw"]], 3), "; ", 
                             round(res.dens[[2]][["bw"]], 3),
                             "\n", "N = ", res.dens[[1]][["n"]], "; ", 
                             res.dens[[2]][["n"]])
  
  
  par(mar=c(7.9, 4.1, 4.1, 2.1))
  if (CV) {
    plot(res[[1]][, 1], res[[1]][, 4], xlab = "", ylab = "CV", 
         type = type, pch = pch, col = col,
         main = "Error plot", 
         xlim = range(res[[1]][, 1], res[[2]][, 1]),
         ylim = range(res[[1]][, 4], res[[2]][, 4])) 
    points(res[[2]][, 1], res[[2]][, 2], pch = pch, col = adjustcolor(col, alpha.f = 0.25))
    
  } else {
    plot(res[[1]][, 1], res[[1]][, 2], 
         ylim = c(min(c(res[[1]][, 2] - res[[1]][, 3], res[[2]][, 2] - res[[2]][, 3])), 
                  max(c(res[[1]][, 2] + res[[1]][, 3], res[[2]][, 2] + res[[2]][, 3]))),
         xlim = range(res[[1]][, 1], res[[2]][, 1]),
         xlab = "", ylab = "MFI", 
         type = type, pch = pch, col = col,
         main = "Error plot")
    points(res[[2]][, 1], res[[2]][, 2], pch = pch, col = adjustcolor(col, alpha.f = 0.25))
    #Plot the location with error bars.
    arrows(res[[1]][, 1], res[[1]][, 2] + res[[1]][, 3], res[[1]][, 1], 
           res[[1]][, 2] - res[[1]][, 3], angle = 90, code = 3, 
           length = length, col = col)
    
    arrows(res[[2]][, 1], res[[2]][, 2] + res[[2]][, 3], res[[2]][, 1], 
           res[[2]][, 2] - res[[2]][, 3], angle = 90, code = 3, 
           length = length, col = adjustcolor(col, alpha.f = 0.25))
    
    deviation.measure <- strsplit(colnames(res[[1]])[3], "(", fixed = TRUE)[[1]][2]
    deviation.measure <- substr(deviation.measure, 1, nchar(deviation.measure) - 1)
    mtext(paste0("Deviation: ", deviation.measure), 4, cex = 0.75)
  }
  
  
  # Add a range for the ROI
  abline(v = llul[[1]], col = "lightgrey", lwd = 1.25)
  abline(v = llul[[2]], col = "lightgrey", lwd = 1.25, lty = 6)
  mtext(error.plot.text, side = 1, line = 7.2, cex = 0.75)
  mtext("Cycle", side = 1, line = 2, cex = 0.8)
  
  #add legend
  legend("topleft", c("A", "B"), pch = c(pch, pch), lty = c(1, 1), 
         col = c(col, adjustcolor(col, alpha.f = 0.25)),
         bg = "white")
  
  par(mar=c(4.1, 4.1, 4.1, 2.1))
  
  # "Calculate" the Quantile-Quantile plots and density plots
  # and plot the results
  plot(res.dens[[1]], xlab = "", main = "Density", col = col,
       xlim = range(c(res.dens[[1]][["x"]], res.dens[[2]][["x"]])),
       ylim = range(c(res.dens[[1]][["y"]], res.dens[[2]][["y"]])))
  lines(res.dens[[2]], col = adjustcolor(col, alpha.f = 0.25))
  mtext("RFU", side = 1, line = 2, cex = 0.8)
  mtext(density.plot.text, side = 1, line = 4.85, cex = 0.75)
  
  # Analysis of the quantiles
  qqnorm1 <- qqnorm(x, plot.it = FALSE)
  qqnorm2 <- qqnorm(y, plot.it = FALSE)
  plot(qqnorm1[["x"]], qqnorm1[["y"]], xlim = range(c(qqnorm1[["x"]], qqnorm2[["x"]])),
         ylim = range(c(qqnorm1[["y"]], qqnorm2[["y"]])), col = col, xlab = "",
         pch = pch, ylab = "Sample Quantiles")
  points(qqnorm2[["x"]], qqnorm2[["y"]], col = adjustcolor(col, alpha.f = 0.25),
         pch = pch)
  mtext("Theoretical Quantiles", side = 1.5, line = 2, cex = 0.75)
  mtext(paste0("\nBreusch-Pagan Test p-value: A ", 
               format(stats[[1]]["heter.p"], digits = 4), "; B ",
               format(stats[[2]]["heter.p"], digits = 4)),
        side = 1, line = 3, cex = 0.75)
  qqline(x, col = col)
  qqline(y, col = adjustcolor(col, alpha.f = 0.25))
  par(default.par)
})



#th class, th.cyc function -------------------------------
setClass("th", contains = "matrix", representation(.Data = "matrix", 
                                                   stats = "summary.lm", 
                                                   input = "matrix"))

setMethod("show", signature(object = "th"), function(object) {
  print(slot(object, ".Data"))
})

setMethod("summary", signature(object = "th"), function(object) {
  cat("Cycle threshold: ", slot(object, ".Data")[, "cyc.th"], "\n")
  cat("Fluorescence threshold: ", slot(object, ".Data")[, "atFluo"], "\n")
  print(slot(object, "stats"))
})

#eff class, effcalc function -------------------------------

setClass("eff", contains = "matrix", representation(.Data = "matrix", 
                                                    amplification.efficiency = "numeric", 
                                                    regression = "lm",
                                                    correlation.test = "htest"))

setMethod("show", signature(object = "eff"), function(object) {
  print(slot(object, ".Data"))
})

setMethod("summary", signature(object = "eff"), function(object) {
  cat("Amplification efficiency: ", slot(object, "amplification.efficiency"), "\n")
  summary(slot(object, "regression"))
})

setMethod("plot", signature(x = "eff"), function(x, xlab = "log10(Concentration)", 
                                                 ylab = "Cq", 
                                                 main = "Efficiency Plot", 
                                                 trend = TRUE, res.fit = "topright", CI = FALSE, 
                                                 level = 0.95, type = "p", 
                                                 pch = 19, er.length = 0.05, 
                                                 col = "black") {
  res <- slot(x, ".Data")
  lm.res <- slot(x, "regression")
  cortest <- slot(x, "correlation.test")
  AE <- slot(x, "amplification.efficiency")
  
  #get significance level
  if (cortest[["p.value"]] < 0.001) {
    sign.out <- "; p < 0.001"
  }
  if (0.001 <= cortest[["p.value"]] && cortest[["p.value"]] < 0.01) {
    sign.out <- "; p < 0.01"
  }
  if (0.01 <= cortest[["p.value"]] && cortest[["p.value"]] < 0.05) {
    sign.out <- "; p < 0.05"
  }
  if (cortest[["p.value"]] >= 0.05) {
    sign.out <- "; p > 0.05"
  }
  
  # Extract coordinates from data matrix
  coords <- c(range(res[, 1]), 
              range(res[, 2]))
  # Store default graphic parameters
  #   default.par <- par()
  
  # Perform a linear regression based on the values of the 
  # calculated mean/median
  # Calculate goodness of fit
  Rsquared <- round(summary(lm.res)[["r.squared"]], 3)
  
  # Plot the Coefficient of Variance
  plot(res[, 1], res[, 2], ylim = c(min(res[, 2] - res[, 3]), 
                                    max(res[, 2] + res[, 3])), xlab = xlab, 
       ylab = ylab, type = type, pch = pch, col = col,
       main = main)
  
  if (CI) {
    # add area and border lines of confidence interval
    # fix, does not yet work properly
    #polygon(c(rev(x.ci), x.ci),
    #  c(predict.ci[, 3], rev(predict.ci[, 2])),
    #  col = "lightgrey", border = NA)
    # prediction and parameters for confidence interval
    x.ci <- seq(coords[1], coords[2], length.out= nrow(res))
    predict.ci <- predict.lm(lm.res, newdata = data.frame(x.ci), 
                             interval = "confidence", level = level)
    polygon(c(x.ci, rev(x.ci)), c(rev(predict.ci[, 3]), predict.ci[, 2]),
            border = NA, col = adjustcolor("lightblue", alpha.f = 0.4))
    #     lines(rev(x.ci), predict.ci[ ,2], col = "lightblue")
    #     lines(rev(x.ci), predict.ci[ ,3], col = "lightblue")
  }
  # Add error bar to the location parameters
  arrow.data <-res[res[, 3] != 0, ] 
  arrows(arrow.data[, 1], arrow.data[, 2] + arrow.data[, 3], 
         arrow.data[, 1], 
         arrow.data[, 2] - arrow.data[, 3], 
         angle = 90, code = 3, length = er.length, 
         col = col)
  
  # Add trend line of linear regression to plot
  if (trend) {
    abline(lm.res)
  }
  
  # Add "legend" with amplification efficiency, goodness of fit to plot
  # ToDo: expression(italic(R)^2 == Rsquared) for ... "R^2 = ", Rsquared ...?
  if (!is.null(res.fit))
    legend(res.fit, paste0("Efficiency = ", AE, " %", "\n",
                              "R^2 = ", Rsquared, "\n",
                              "r = ", round(cortest[["estimate"]], 3), sign.out, "\n"),
           bty = "n")
  
  #   # Restore default graphic parameters
  #   par(default.par)
})