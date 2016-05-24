plot.basta <-
    function(x, plot.trace = TRUE, trace.name = "theta", fancy = FALSE,
        ...){
  # Define number of variables and colors:
  Palette <- c('#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', 
      '#FFFF33', '#A65628', '#F781BF', '#999999')
  nv <- ifelse(fancy, length(x$survQuant), 
      ifelse(!plot.trace, length(x$survQuant), x$settings['nsim']))
  if ("col" %in% names(list(...))) {
    if (length(list(...)$col) < nv) {
      op <- options()
      options(warn = 1)
      ncwarn <- ifelse(fancy, "covariates", ifelse(!plot.trace, "covariates",
              "simulations"))
      warning(sprintf("Insufficient number of colors. Not all %s will be displayed.",
              ncwarn), call. = FALSE)
      options(op)
      Bord <- list(...)$col[1:nv]
    } else {
      Bord <- list(...)$col[1:nv]
    }
  } else {
    if (nv <= 9) {
      Bord <- Palette[1:nv]
    } else {
      Bord <- rainbow(nv)
    }
  }
  Cols <- adjustcolor(Bord, alpha.f = 0.5)
  
  if (fancy) {
    catNames <- names(x$survQuant)      
  } else {
    if (plot.trace) {
      catNames <- NA
    } else {
      catNames <- names(x$survQuant)
    }
  }
  if ("names.legend" %in% names(list(...))) {
    if (length(list(...)$names.legend) != nv) {
      stop(sprintf("Wrong length in 'names.legend'. Correct length is %s elements.",
              nv), call. = FALSE)
    } else {
      legNames <- list(...)$names.legend
    }
    
  } else {
    legNames <- substr(catNames, 2, nchar(catNames))
  }
  if (fancy) {
    allThetaNames <- c("a0", "a1", "c", "b0", "b1", "b2")
    allThetaExpr  <- expression(italic(a[0]), italic(a[1]), italic(c), 
                                italic(b[0]), italic(b[1]), italic(b[2]))
    #catNames <- names(x$survQuant)
    lenCat <- length(catNames)
    varNames <- substr(colnames(x$params), 1, 2)
    idTh <- is.element(substr(varNames, 1, 1), c("a", "b", "c"))
    thetaMat <- x$params[,idTh]
    model <- as.character(x$modelSpecs['model'])
    shape <- as.character(x$modelSpecs['shape'])
    
    if (model == "EX" | model == 1) {
      nTheta          <- 1
      idAllTheta      <- 4
    } else if (model %in% c("GO", "WE") | model %in% c(2, 3)) {
      nTheta          <- 2
      idAllTheta      <- c(4, 5)
    }  else {
      nTheta          <- 3
      idAllTheta      <- c(4, 5, 6)
    }
    if (shape == "Makeham" | shape == 2) {
      nTheta          <- nTheta + 1
      idAllTheta      <- c(3, idAllTheta)
    } else if(shape == "bathtub" | shape == 3) {
      nTheta          <- nTheta + 3
      idAllTheta      <- c(1:3, idAllTheta)
    }
    lows <- c(-Inf, 0, 0, -Inf, 0, 0)
    names(lows) <- allThetaNames
    
    # Build plots:  
    op <- par(no.readonly = TRUE)
    layout(mat = matrix(data = c(rep(1:nTheta, each = 2), 
                rep(rep(c(nTheta + 1, nTheta + 2), 
                        each = nTheta), 2)), 
            nrow  = nTheta * 2, 
            ncol  = 3), 
        widths  = rep(2, 3), 
        heights = rep(1, nTheta))
    par(mar=c(3, 3, 0.5, 0.5))
    for(i in 1:nTheta){
      dez <- list()
      ylz <- rep(NA, lenCat)
      xlz <- matrix(0, lenCat, 2)
      for(j in 1:lenCat){
        idth <- paste(allThetaNames[idAllTheta[i]],
                      catNames[j], sep = "")
        if (lenCat == 1) {
          idth <- allThetaNames[idAllTheta[i]]
        }
        dez[[catNames[j]]] <- density(thetaMat[, idth])
        ylz[j] <- max(dez[[catNames[j]]]$y)
        xlz[j,] <- range(dez[[j]]$x)
      }
      xr <- range(xlz)
      xl <- c(floor(max(c(lows[i], xr[1])) * 10) / 10, ceiling(xr[2] * 10) / 10)
      xd <- ceiling(diff(xl) * 10) / 10
      plot(x = dez[[1]], xlab = "", ylab = "", xlim = xl, ylim = c(0, max(ylz)), 
          lwd = 3, axes = FALSE, main = "", col = NA)
      for(j in 1:lenCat) {
        polygon(x = c(dez[[j]]$x, dez[[j]]$x[1]), 
            y = c(dez[[j]]$y, dez[[j]]$y[1]), 
            col = Cols[j], border = Bord[j], lwd = 1.5)
      }
      axis(side = 1, at = seq(xl[1], xl[2], length = 5), 
          line = 0.5, labels = NA, tcl = 0.4)
      axis(side = 1, at = seq(xl[1], xl[2], length = 3), lwd = NA)
      mtext(text = allThetaExpr[idAllTheta[i]], side = 2, line = 0, 
          at = max(ylz) * 0.8, las = 2, cex = 1.25)
    }
    
    # Plot survival probability:
    par(mar = c(4, 7, 0.5, 0.5))
    xv <- lapply(1:lenCat, 
        function(idcovs) as.numeric(colnames(x$mortQuant[[idcovs]])))
    mxv <- ceiling(max(unlist(xv)) / 5) * 5
    plot(x = c(0, mxv), y = range(0, 1), col = NA, axes = FALSE, xlab = "", 
        ylab = "")
    for(i in 1:lenCat){
      polygon(x = c(xv[[i]], rev(xv[[i]])), 
          y = c(x$survQuant[[i]][2, ], rev(x$survQuant[[i]][3, ])), 
          col = Cols[i], border = Bord[i])
      lines(x = xv[[i]], y = x$survQuant[[i]][1, ], col = Bord[i], lty = 3)
    }
    if (lenCat > 1) {
      legend(x = 'topright', legend = legNames, 
          pch = 15, pt.cex = 3, cex = 1.5, col = Cols, bty = 'n')
    }
    
    axis(side = 1, at = seq(0, mxv, 5), labels = NA, tcl = 0.4, line = 0.5)
    axis(side = 2, at = seq(0, 1, 0.2), tcl = 0.4, las = 2, cex.axis = 1.2)
    mtext(text = expression(paste("Survival, ", italic(S(x)))), 
        side = 2, line = 3.5, cex = 1.25)
    
    # Plot mortality rates:
    #ylmx <- c(0, max(c(0.5, round(max(unlist(x$mortQuant))))))
    ylmx <- c(0, max(c(0.5, ceiling(max(unlist(x$mortQuant)/ 10) * 10)) ))
    plot(x = c(0, mxv), y = ylmx, col = NA, ylim = ylmx, xlab = "", ylab = "",
        axes = FALSE)
    for(i in 1:lenCat){
      polygon(x = c(xv[[i]], rev(xv[[i]])), 
          y = c(x$mortQuant[[i]][2, ], rev(x$mortQuant[[i]][3, ])), 
          col = Cols[i], border = Bord[i])
      lines(x = xv[[i]], y = x$mortQuant[[i]][1, ], col = Bord[i], lty = 3)
    }
    axis(side = 1, at = seq(0, mxv, 5), labels = NA, tcl = 0.4, line = 0.5)
    axis(side = 1, at = seq(0, mxv, 5), lwd = NA, cex.axis = 1.2)
    axis(side = 2, at = seq(0, ylmx[2], length = 5), 
        tcl = 0.4, las = 2, cex.axis = 1.2)
    mtext(text = expression(paste("Mortality, ", italic(mu(x)))), 
        side = 2, line = 3.5, cex = 1.25)
    mtext(text = expression(paste("Age ", italic(x), " (years)")), 
        side = 1, cex = 1.25, line = 3)
    par(op)
  } else {
    # NOT FANCY...
    # 1. Trace plots for parameters:
    if(plot.trace) {
      traces <- c("theta","gamma","pi")
      parNames <- substr(colnames(x$params), 1, 2)
      if (trace.name == "theta") {
        idPars <- which(!(parNames %in% c("ga", "pi")))
      } else if (trace.name == "gamma") {
        if (trace.name == "gamma" & x$modelSpecs[3] %in% c("inMort", "noCov")) {
          stop("\nTrace plots cannot be drawn for 'gamma' parameters.",
              "\nNo proportional hazards arguments were evaluated.", 
              call. = FALSE)
        } else {
          idPars <- which(parNames == "ga")
        }
      } else if (trace.name == "pi") {
        idPars <- which(parNames == "pi")
      } else {
        stop("Wrong 'trace.name' argument. Valid arguments are:\n", 
            "'theta', 'gamma' or 'pi'.\n", call. = FALSE)
      }
      nRows <- ifelse(trace.name == "theta", length(unique(parNames[idPars])),
          ceiling(length(idPars) / 2))
      nCols <- ifelse(trace.name == "theta", table(parNames[idPars])[1], 
          ceiling(length(idPars) / nRows))
      yLim <- sapply(colnames(x$params)[idPars], function(par) 
            range(sapply(1:x$settings["nsim"], function(sim)
                      range(x$parsForPlot[[sim]][, par]))))
      op <- par(mfrow = c(nRows, nCols), mar = c(3, 3, 2.5, 1))
      for (i in 1:length(idPars)) {
        plot(c(1, x$settings['niter']), col = NA, yLim[, i], xlab = "", 
            ylab = "", frame.plot = FALSE, main = colnames(x$params)[idPars[i]])
        for (j in 1:x$settings['nsim']) {
          lines(seq(1, x$settings['niter'], x$settings['thinning']), 
              x$parsForPlot[[j]][, idPars[i]], col = Bord[j])
        }
      }
      par(op)
      
      # 2. Plot survival and mortality:
    } else {
      nv <- length(x$survQuant)
      op <- par(mfrow = c(2, 1), mar = c(4, 4, 3, 2))
      xv <- lapply(1:nv, 
          function(idcovs) as.numeric(colnames(x$mortQuant[[idcovs]])))
      if (!("xlim" %in% names(list(...)))) {
        xlim <- c(0, max(unlist(xv)))
        idRange <- lapply(1:nv, function(idcovs) length(xv[[idcovs]]))
      } else {
        xlim <- list(...)$xlim
        idRange <- lapply(1:nv, 
            function(idcovs) which(abs(xv[[idcovs]] - xlim[2]) == 
                      min(abs(xv[[idcovs]] - xlim[2])))[1])
        xv <- lapply(1:nv, 
            function(idcovs) xv[[idcovs]][1:idRange[[idcovs]]])
      }
      if ("noCI" %in% names(list(...))) {
        noCI <- list(...)$noCI
        if (class(noCI) != "logical") {
          noCI <- FALSE
        }
      } else {
        noCI <- FALSE
      }
      for (dem in c("survQuant", "mortQuant")) {
        if (noCI) {
          idQuants <- 1
        } else {
          idQuants <- 1:3
        }
        ylim <- c(0, max(sapply(1:length(x[[dem]]), function(idem) 
                      max(x[[dem]][[idem]][idQuants, 2:idRange[[idem]]]))))
        xlab <- ifelse(dem == "mortQuant", "Age", "")
        ylab <- ifelse(dem == "mortQuant", expression(mu(x)), expression(S(x)))
        main <- ifelse(dem == "mortQuant", "Mortality", "Survival")
        plot(xlim, ylim, col = NA, xlab = xlab, ylab = ylab, 
            frame.plot = FALSE, main = main, ylim = c(0, ylim[2]))
        if (!noCI) {
          for (cov in 1:length(x$mortQuant)) {
            polygon(c(xv[[cov]][1:idRange[[cov]]], 
                    rev(xv[[cov]][1:idRange[[cov]]])), 
                c(x[[dem]][[cov]][2, 1:idRange[[cov]]], 
                    rev(x[[dem]][[cov]][3, 1:idRange[[cov]]])),
                col = Cols[cov], border = Bord[cov], lwd = 0.5, lty = 1)
          }
        }
        for (cov in 1:length(x$mortQuant)) {
          lines(xv[[cov]][1:idRange[[cov]]], 
              x[[dem]][[cov]][1, 1:idRange[[cov]]], col = Bord[cov], lwd = 3)
        }
        if (dem == "survQuant" & length(x$mortQuant) > 1) {
          legend('topright', legNames, 
              col = Bord, lwd = 3, bty = 'n')
        }
      }
      par(op)
    }
  }
}