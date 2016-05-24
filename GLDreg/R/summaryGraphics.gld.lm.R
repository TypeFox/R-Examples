summaryGraphics.gld.lm<-
  function (overall.fit.obj, alpha = 0.05, label = NULL, ColourVersion = TRUE, 
            diagnostics = TRUE, range = c(0.01, 0.99)) 
  {
    fit.obj <- overall.fit.obj[[1]]
    simu.obj <- overall.fit.obj[[3]]
    parEstimate <- fit.obj[[3]][match(dimnames(simu.obj)[[2]], 
                                      names(fit.obj[[3]]))]
    num.par <- length(parEstimate)
    if (is.null(label)) {
      parNames <- names(parEstimate)
    }
    else {
      parNames <- label[match(dimnames(simu.obj)[[2]], names(fit.obj[[3]]))]
    }
    empty.panel <- function() {
      plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), 
           xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "o")
      invisible()
    }
    columnCols <- c("orange", "red", "navy", "green4")
    columnHeadings <- c("Parameter", "Estimate", "Density", paste("Summary \n", 
                                                                  paste((1 - alpha) * 100, "%", sep = ""), "interval:"))
    if (!ColourVersion) {
      columnCols <- rep("black", num.par)
    }
    headInds <- 1:4
    panels.per.par <- 4
    par(mfrow = c((num.par + 1), panels.per.par))
    for (i in 1:panels.per.par) {
      par(ann = F, mar = rep(0, 4), xaxt = "n", yaxt = "n", 
          xpd = TRUE)
      empty.panel()
      text(0.5, 0.5, columnHeadings[headInds[i]], cex = 3, 
           col = columnCols[headInds[i]])
    }
    for (j in 1:num.par) {
      par(ann = F, mar = c(0, 0, 0, 0), xaxt = "n", yaxt = "n")
      empty.panel()
      text(0.5, 0.5, parNames[[j]][1], cex = 2.5, col = columnCols[1])
      empty.panel()
      text(0.5, 0.5, signif(parEstimate[[j]][1], 3), cex = 2.5, 
           col = columnCols[2])
      par(ann = T, xaxt = "n", mar = c(2, 2, 1, 1), yaxt = "n", 
          xpd = TRUE)
      if (sum(columnCols == "black") == num.par) {
        hist(as.numeric(simu.obj[, j]), prob = T, nclass = "scott", 
             main = "", xlab = "", xaxs = "i", ylab = "", 
             col = "white", border = "black", col.axis = "black")
      }
      if (sum(columnCols == "black") != num.par) {
        hist(as.numeric(simu.obj[, j]), prob = T, nclass = "scott", 
             main = "", xlab = "", xaxs = "i", ylab = "", 
             col = columnCols[3], border = "pink", col.axis = columnCols[3])
      }
      box("figure")
      par(ann = F, mar = c(0, 0, 0, 0), xaxt = "n", yaxt = "n")
      empty.panel()
      interval <- signif(quantile(simu.obj[, j], c(0.5 - (1 - 
                                                            alpha)/2, 0.5 + (1 - alpha)/2), type = 8), 3)
      text(0.5, 0.5, paste(interval, collapse = ",  "), cex = 2, 
           col = columnCols[4])
    }
    if (diagnostics) {
      dev.new()
      resid <- fit.obj[[5]]
      param <- fit.obj[[7]]
      gld.values <- fit.obj[[3]][(length(fit.obj[[3]]) - 3):length(fit.obj[[3]])]
      par(mfrow = c(1, 2))
      qqplot.gld(resid, gld.values, param = param, range = range, 
                 main = "QQ plot for residuals Version 1")
      legend("bottomright", c(paste(toupper(param), "GLD"), 
                              paste("(", paste(signif(gld.values, 3), collapse = ","), 
                                    ")", sep = "")), bty = "n")
      qqplot.gld(resid, gld.values, param = param, range = range, 
                 type = "str.qqplot", main = "QQ plot for residuals Version 2")
      pval <- ks.gof(resid, "pgl", lambda1 = gld.values, param = param)$p.value
      legend("topleft", paste("KS test\np-value=", format.pval(pval)), 
             bty = "n")
    }
  }
