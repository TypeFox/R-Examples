plot.simfrail <- function(x, type=c("residuals","hazard"), ...) {
  sim <- x
  if (type == "residuals") {
    plot.simfrail.residuals(sim, ...)
  } else if (type == "hazard") {
    plot.simfrail.hazard(sim, ...)
  }
}

plot.simfrail.residuals <- function(sim, n.Lambda=3, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || 
      !requireNamespace("reshape2", quietly = TRUE)) {
    stop("Plotting requires the ggplot2, and reshape2 packages")
  }
  Lambda.cols <- names(sim)[grepl("^Lambda", names(sim))]
  
  # All BUT n.Lambda
  if (n.Lambda < 0) {
    n.Lambda <- max(length(Lambda.cols) + n.Lambda, 0)
  } 
  
  if (n.Lambda == 0) {
    # No Lambda
    Lambda.cols <- NULL
  } else if (length(Lambda.cols) > n.Lambda) {
    # Evenly spaced n.Lambda
    idx <- round(seq(0, length(Lambda.cols), 
                     length.out=(n.Lambda+2))[2:(n.Lambda+1)])
    
    Lambda.cols <- Lambda.cols[idx]
  }
  hat.cols <- c(names(sim)[grepl("^hat.beta|^hat.theta", names(sim))], 
                paste("hat.", Lambda.cols, sep=""))
  value.cols <- c(names(sim)[grepl("^beta|^theta", names(sim))], Lambda.cols)
  
  residuals <- data.frame(cbind(N=sim$N,
                     vapply(value.cols,
    function(col) sim[[paste("hat.", col, sep="")]] - sim[[col]], rep(0, nrow(sim)))))
  
  # Select N and residuals columns, start with res
  n.vars <- length(value.cols)
  
  res.melt <- melt(residuals, id = c("N"))
  cases <- c(t(unique(res.melt["N"])))
  res.melt$x <- factor(0.5*floor(res.melt$N))
  
  p <- ggplot(res.melt, aes_string(x='x', y='value', fill='variable')) + 
    geom_boxplot(notch=TRUE) +
    facet_grid(.~variable) +
    labs(x="N",y="Residual") + 
    theme(legend.position="none",axis.text.x=element_text(angle=-90, vjust=0.4,hjust=1))
  
  p
}

plot.simfrail.hazard <- function(sim, CI=0.95, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || 
      !requireNamespace("reshape2", quietly = TRUE)) {
    stop("Plotting requires the ggplot2, reshape2 packages")
  }
  
  hats <- sim[,grepl("^hat.Lambda", names(sim))]
  se <- sim[,grepl("^se.Lambda", names(sim))]
  values <- colMeans(sim[,grepl("^Lambda", names(sim))])
  
  Lambda.times <- sapply(names(values), function(w) gsub(".+\\.", "", w), USE.NAMES=FALSE)
  Lambda.times <- as.numeric(Lambda.times)
  names(hats) <- Lambda.times
  melthats <- melt(t(hats))
  names(melthats) <- c("Time","instance", "value")
  melthats$type <- sprintf("Empirical (%.2f CI)", CI)
  
  values <- data.frame(x=Lambda.times, y=values)
  values$type <- "Actual"
  
  Z.score <- qnorm((1-CI)/2)
  mean.se <- colMeans(se)
  se <- data.frame(x=Lambda.times, lower=values$y-Z.score*mean.se, upper=values$y+Z.score*mean.se)
  se$type <- sprintf("Estimated %.2f CI", CI)
  
  p <- ggplot(melthats, aes_string(x='Time',y='value',color='type')) +
    stat_summary(fun.data=mean_cl_boot, geom="smooth") +
    geom_line(aes_string(x='x', y='y', color='type'), values) +
    theme(legend.position=c(0,1),
          legend.justification=c(0,1)) +
    ylab("Cumulative baseline hazard")
  
  if (all(is.na(mean.se))) {
    p <- p + scale_colour_manual("Legend", values=c("black","blue"))
  } else {
    p <- p + 
          geom_line(aes_string(x='x', y='upper', color='type'), se) +
          geom_line(aes_string(x='x', y='lower', color='type'), se) +
          scale_colour_manual("Legend", values=c("black","blue","brown"))
  }
  
  p
}