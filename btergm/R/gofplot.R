# This file contains print and plot methods for various gof objects.

# print method for boxplot gof objects
print.boxplot <- function(x, ...) {
  message(x$label, "\n")
  printCoefmat(x$stats, has.Pvalue = TRUE, cs.ind=numeric(0), 
      signif.legend = FALSE)
  message(paste("\nNote: Small p-values indicate a significant difference",
      "\n      between simulations and observed network(s)."))
}

# print method for roc gof objects
print.roc <- function(x, ...) {
  message(x$label, "\n")
  auc <- cbind(x$auc.roc, x$auc.roc.rgraph)
  colnames(auc) <- c("ROC model", "ROC random")
  rownames(auc) <- 1:nrow(auc)
  print(auc)
}

# print method for roc gof objects
print.pr <- function(x, ...) {
  message(x$label, "\n")
  auc <- cbind(x$auc.pr, x$auc.pr.rgraph)
  colnames(auc) <- c("PR model", "PR random")
  rownames(auc) <- 1:nrow(auc)
  print(auc)
}

# print method for roc gof objects
print.rocpr <- function(x, ...) {
  message(x$label, "\n")
  auc <- cbind(x$auc.roc, x$auc.roc.rgraph, x$auc.pr, x$auc.pr.rgraph)
  colnames(auc) <- c("ROC model", "ROC random", "PR model", "PR random")
  rownames(auc) <- 1:nrow(auc)
  print(auc)
}

# print method for univariate gof objects
print.univariate <- function(x, ...) {
  message(x$label, "\n")
  message("Observed:")
  print(summary(x$obs))
  message("\nSimulated:")
  print(summary(x$sim))
}

# print method for gof objects
print.gof <- function(x, ...) {
  for (i in 1:length(x)) {
    cat("\n")
    print(x[[i]])
  }
}

# plot function for gof objects (= lists of boxplot etc. objects)
plot.gof <- function(x, mfrow = TRUE, ...) {
  if (mfrow == TRUE) {
    l <- length(x)
    if (l == 0) {
      stop("'x' does not contain any gof objects.")
    } else if (l == 1) {
      # do nothing
    } else if (l == 2) {
      par(mfrow = c(1, 2))
    } else if (l == 3) {
      par(mfrow = c(1, 3))
    } else if (l == 4) {
      par(mfrow = c(2, 2))
    } else if (l == 5 || l == 6) {
      par(mfrow = c(2, 3))
    } else if (l > 6 && l < 10) {
      par(mfrow = c(3, 3))
    } else if (l > 9 && l < 13) {
      par(mfrow = c(4, 3))
    } else if (l > 12 && l < 17) {
      par(mfrow = c(4, 4))
    } else if (l > 16 && l < 21) {
      par(mfrow = c(5, 4))
    } else if (l > 20 && l < 26) {
      par(mfrow = c(5, 5))
    } else {
      warning("Too many gof objects for tiling. Switching off 'mfrow'.")
    }
  }
  for (i in 1:length(x)) {
    plot(x[[i]], ...)
  }
}

# plot function for boxplot objects
plot.boxplot <- function(x, relative = TRUE, transform = function(x) x, 
    xlim = NULL, main = x$label, xlab = x$label, ylab = "Frequency", 
    border = "darkgray", boxplot.lwd = 0.8, outline = FALSE, median = TRUE, 
    median.col = "black", median.lty = "solid", median.lwd = 2, mean = TRUE, 
    mean.col = "black", mean.lty = "dashed", mean.lwd = 1, ...) {
  
  # transform data
  mat <- t(x$raw)
  mat <- apply(mat, 1:2, transform)
  if (ncol(x$stats == 9)) {
    x$stats[, 2] <- transform(x$stats[, 2])
  }
  x$stats[, 1] <- transform(x$stats[, 1])
  
  # convert to relative frequencies
  if (relative == TRUE) {
    rs <- rowSums(mat)
    for (i in 1:nrow(mat)) {
      if (rs[i] == 0) {
        mat[i, ] <- 0
      } else {
        mat[i, ] <- mat[i, ] / rs[i]
      }
    }
    cs <- colSums(x$stats)[1:2]
    for (i in 1:2) {
      if (cs[i] == 0) {
        x$stats[, i] <- 0
      } else {
        x$stats[, i] <- x$stats[, i] / cs[i]
      }
    }
  }
  
  # find minimum and maximum values for plotting
  if (ncol(x$stats) == 9) { # several observed networks
    obs.mean <- x$stats[, 1]
    obs <- x$stats[, 2]
  } else { # only one observed network
    obs <- x$stats[, 1]
  }
  if (any(is.infinite(c(mat, obs)))) {
    stop(paste("Simulated or observed values contain infinite values.", 
        "Check the 'transform' argument."))
  }
  
  # impose xlim
  if (!is.null(xlim)) {
    if (colnames(mat)[ncol(mat)] == "Inf") {
      if (xlim > nrow(x$stats) - 2) {
        xlim <- nrow(x$stats) - 1
        warning("'xlim' was out of bounds. Replaced by maximum possible.")
      }
      mat <- mat[, c(1:(xlim), ncol(mat))]
      obs <- obs[c(1:(xlim), length(obs))]
      if (exists("obs.mean")) {
        obs.mean <- obs.mean[c(1:(xlim), length(obs.mean))]
      }
    } else {
      if (xlim > ncol(mat) - 1) {
        warning("'xlim' was out of bounds. Replaced by maximum possible.")
      } else {
        mat <- mat[, c(1:(xlim + 1))]
        obs <- obs[c(1:(xlim + 1))]
        if (exists("obs.mean")) {
          obs.mean <- obs.mean[c(1:(xlim + 1))]
        }
      }
    }
  }
  
  # set empty main title if NULL
  if (is.null(main)) {
    main <- ""
  }
  
  # plot boxplots and curves
  boxplot(mat, border = border, xlab = xlab, ylab = ylab, main = main, 
      outline = outline, lwd = boxplot.lwd, ...)
  if (ncol(x$stats) == 9 && mean == TRUE) {
    lines(obs.mean, lwd = mean.lwd, type = "l", lty = mean.lty, col = mean.col)
  }
  if (median == TRUE) {
    lines(obs, lwd = median.lwd, type = "l", lty = median.lty, col = median.col)
  }
}

# plot function for roc objects
plot.roc <- function(x, add = FALSE, main = x$label, avg = c("none", 
    "horizontal", "vertical", "threshold"), spread.estimate = c("boxplot", 
    "stderror", "stddev"), lwd = 3, rgraph = TRUE, col = "#bd0017", 
    random.col = "#bd001744", ...) {
  
  plot(x$roc, avg = avg[1], spread.estimate = spread.estimate[1], 
      add = add, col = col, main = main, lwd = lwd, 
      boxplot.boxcol = col, ylim = c(0, 1), ...)
  if (rgraph == TRUE) {
    plot(x$roc.rgraph, avg = avg[1], spread.estimate = "none", 
        col = random.col, add = TRUE, lwd = lwd, ylim = c(0, 1), ...)
  }
}

# plot function for pr objects
plot.pr <- function(x, add = FALSE, main = x$label, avg = c("none", 
    "horizontal", "vertical", "threshold"), spread.estimate = c("boxplot", 
    "stderror", "stddev"), lwd = 3, rgraph = TRUE, col = "#5886be", 
    random.col = "#5886be44", pr.poly = 0, ...) {
  
  plot(x$pr, avg = avg[1], spread.estimate = spread.estimate[1], 
      add = add, col = col, main = main, lwd = lwd, 
      boxplot.boxcol = col, ylim = c(0, 1), ...)
  if (rgraph == TRUE) {
    plot(x$pr.rgraph, avg = avg[1], spread.estimate = "none", 
        col = random.col, add = TRUE, lwd = lwd, ylim = c(0, 1), ...)
  }
  
  # fit polynomial curve through PR curve
  if (pr.poly > 0) {
    for (i in 1:length(x$pr@y.values)) {
      prec <- x$pr@y.values[[i]]
      prec[1] <- NaN
      rec <- x$pr@x.values[[i]]
      p <- data.frame(poly(rec, pr.poly, raw = TRUE))
      fit <- lm(prec ~ ., data = p)
      yhat <- predict(fit, newdata = p)
      for (j in 1:length(yhat)) {
        if (yhat[j] > 1) {
          yhat[j] <- 1
        }
        if (yhat[j] < 0) {
          yhat[j] <- 0
        }
      }
      lines(rec, yhat, type = "l", lty = "dashed", col = "red", 
          lwd = lwd)
      points(rec[1], yhat[1], pch = 1, col = "red", cex = 2)
      points(x$pr@x.values[[i]][1], x$pr@y.values[[i]][1], pch = 1, 
          col = col, cex = 2)
    }
  }
}

# plot function for rocpr objects
plot.rocpr <- function(x, main = x$label, roc.avg = c("none", "horizontal", 
    "vertical", "threshold"), roc.spread.estimate = c("boxplot", "stderror", 
    "stddev"), roc.lwd = 3, roc.rgraph = TRUE, roc.col = "#bd0017", 
    roc.random.col = "#bd001744", pr.avg = c("none", "horizontal", 
    "vertical", "threshold"), pr.spread.estimate = c("boxplot", "stderror", 
    "stddev"), pr.lwd = 3, pr.rgraph = TRUE, pr.col = "#5886be", 
    pr.random.col = "#5886be44", pr.poly = 0, ...) {
  
  plot(x$roc, avg = roc.avg[1], spread.estimate = roc.spread.estimate[1], 
      col = roc.col, main = main, lwd = roc.lwd, boxplot.boxcol = roc.col, 
      ylim = c(0, 1), ylab = "TPR / PPV", xlab = "FPR / TPR", ...)
  if (roc.rgraph == TRUE) {
    plot(x$roc.rgraph, avg = roc.avg[1], spread.estimate = "none", 
        col = roc.random.col, add = TRUE, lwd = roc.lwd, ylim = c(0, 1), ...)
  }
  
  plot(x$pr, avg = pr.avg[1], spread.estimate = pr.spread.estimate[1], 
      col = pr.col, main = main, lwd = pr.lwd, boxplot.boxcol = pr.col, 
      ylim = c(0, 1), add = TRUE, ...)
  if (pr.rgraph == TRUE) {
    plot(x$pr.rgraph, avg = pr.avg[1], spread.estimate = "none", 
        col = pr.random.col, add = TRUE, lwd = pr.lwd, ylim = c(0, 1), ...)
  }
  
  # fit polynomial curve through PR curve
  if (pr.poly > 0) {
    for (i in 1:length(x$pr@y.values)) {
      prec <- x$pr@y.values[[i]]
      prec[1] <- NaN
      rec <- x$pr@x.values[[i]]
      p <- data.frame(poly(rec, pr.poly, raw = TRUE))
      fit <- lm(prec ~ ., data = p)
      yhat <- predict(fit, newdata = p)
      for (j in 1:length(yhat)) {
        if (yhat[j] > 1) {
          yhat[j] <- 1
        }
        if (yhat[j] < 0) {
          yhat[j] <- 0
        }
      }
      lines(rec, yhat, type = "l", lty = "dashed", col = "red", 
          lwd = pr.lwd)
      points(rec[1], yhat[1], pch = 1, col = "red", cex = 2)
      points(x$pr@x.values[[i]][1], x$pr@y.values[[i]][1], pch = 1, 
          col = pr.col, cex = 2)
    }
  }
}

# plot function for univariate objects
plot.univariate <- function(x, main = x$label, sim.hist = TRUE, sim.bar = TRUE, 
    sim.density = TRUE, obs.hist = FALSE, obs.bar = TRUE, obs.density = TRUE, 
    sim.adjust = 1, obs.adjust = 1, sim.lwd = 2, obs.lwd = 2, 
    sim.col = "black", obs.col = "red", ...) {
  d.sim <- density(x$sim, adjust = sim.adjust)
  d.sim.max.y <- max(d.sim$y)
  d.sim.max.x <- max(d.sim$x)
  d.sim.min.y <- min(d.sim$y)
  d.sim.min.x <- min(d.sim$x)
  d.sim.median <- median(x$sim)
  if (length(x$obs) > 1) {
    d.obs <- density(x$obs, adjust = obs.adjust)
    d.obs.max.y <- max(d.obs$y)
    d.obs.max.x <- max(d.obs$x)
    d.obs.min.y <- min(d.obs$y)
    d.obs.min.x <- min(d.obs$x)
    d.obs.median <- median(x$obs)
  } else {
    d.obs.max.y <- 0
    d.obs.max.x <- x$obs
    d.obs.min.y <- 0
    d.obs.min.x <- x$obs
  }
  h <- hist(x$sim, plot = FALSE)
  y.max <- max(c(d.sim.max.y, d.obs.max.y, h$density))
  x.max <- max(c(d.sim.max.x, d.obs.max.x, h$breaks))
  y.min <- min(c(d.sim.min.y, d.obs.min.y, h$density))
  x.min <- min(c(d.sim.min.x, d.obs.min.x, h$breaks))
  plot(0, xlim = c(x.min, x.max), ylim = c(y.min, y.max), 
      xlab = x$label, ylab = "Density", main = main, type = "n")
  if (sim.hist == TRUE) {
    hist(x$sim, add = TRUE, lwd = sim.lwd, freq = FALSE)
  }
  if (sim.bar == TRUE) {
    lines(rep(d.sim.median, 2), c(y.min, y.max), col = sim.col, lwd = sim.lwd)
  }
  if (sim.density == TRUE) {
    lines(d.sim, lwd = sim.lwd)
  }
  if (obs.hist == TRUE) {
    hist(x$obs, add = TRUE, lwd = sim.lwd, freq = FALSE, border = obs.col)
  }
  if (obs.bar == TRUE) {
    if (length(x$obs) > 1) {
      lines(rep(d.obs.median, 2), c(y.min, y.max), col = obs.col, lwd = obs.lwd)
    } else {
      lines(rep(x$obs, 2), c(y.min, y.max), col = obs.col, lwd = obs.lwd)
    }
  }
  if (obs.density == TRUE) {
    if (length(x$obs) > 1) {
      lines(d.obs, col = obs.col, lwd = obs.lwd)
    }
  }
}
