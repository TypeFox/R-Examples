plot.fitfrail <- function(x, type=c("hazard", "trace"), ...) {
  fit <- x
  if (type == "hazard") {
    plot.fitfrail.hazard(fit, ...)
  } else if (type == "trace") {
    plot.fitfrail.trace(fit, ...)
  }
}

plot.fitfrail.trace <- function(fit, show.loglik=TRUE, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || 
      !requireNamespace("reshape2", quietly = TRUE)) {
    stop("Plotting requires the ggplot2, reshape2 packages")
  }
  
  trace <- fit$VARS$trace
  loglik <- trace[,c("Iteration", "loglik")]
  trace$loglik <- NULL
  
  melttrace <- melt(trace, "Iteration", variable.name="Parameter")
  breaks <- seq(1, nrow(trace), by=max(1, round(nrow(trace)/10)))
  
  p <- ggplot(melttrace, aes_string(x='Iteration', y='value', color='Parameter')) +
    geom_line() +
    xlab("Iteration") + 
    ylab("Estimate") + 
    scale_x_continuous(breaks = breaks) +
    theme(legend.justification = c(1, 0.5), legend.position = c(1, 0.5)) +
    ggtitle("Parameter estimate trace")
  
  if (show.loglik) {
    if (!requireNamespace("gridExtra", quietly = TRUE)) {
      stop("Plotting requires the gridExtra package")
    }
    p2 <- ggplot(loglik, aes_string(x='Iteration', y='loglik')) +
      geom_line() + 
      xlab("Iteration") + 
      ylab("Log-liklihood") + 
      scale_x_continuous(breaks = breaks) +
      ggtitle("Log-likelihood trace")
    
    gridExtra::grid.arrange(p, p2, ncol=2)
  } else {
    p
  }
}

plot.fitfrail.hazard <- function(fit, CI=0, end=NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Plotting requires the ggplot2 package")
  }
  
  Lambda <- fit$Lambda
  
  if (is.null(end)) {
    end <- Lambda$time[nrow(Lambda)] + mean(diff(fit$Lambda$time))
  }
  
  rownames(Lambda) <- 1:nrow(Lambda)
  ymax <- max(Lambda$Lambda[Lambda$time <= end])
  Lambda <- Lambda[Lambda$Lambda<=ymax,]
  
  p <- qplot(time, Lambda, data = Lambda, geom="step") +
    scale_x_continuous(breaks=seq(0, round(end), by=max(10, round(end/10)))) +
    geom_rug(sides="b", size=0.5) +
    xlab("Time") + 
    ylab("Cumulative baseline hazard") + 
    theme(legend.position="none")
  
  if ((CI > 0)&&(CI < 1)) {
    COV <- vcov(fit, boot=TRUE, ...)
    se.Lambda <- sqrt(diag(COV))[(fit$VARS$n.gamma+1):nrow(COV)]
    se.Lambda <- se.Lambda[1:nrow(Lambda)]
    
    Z.score <- qnorm(1-(1-CI)/2)
    UB <- Lambda$Lambda + Z.score*se.Lambda
    LB <- Lambda$Lambda - Z.score*se.Lambda
    LB <- vapply(LB, function(x) max(0, x), 0)
    
    time <- rep(Lambda$time, rep(2, nrow(Lambda)))
    time <- c(time[1:(length(time)-1)], rev(time[1:(length(time)-1)]))
    
    LB <- rep(LB, rep(2, nrow(Lambda)))[2:(2*nrow(Lambda))]
    UB <- rep(UB, rep(2, nrow(Lambda)))[2:(2*nrow(Lambda))]
    
    df <- data.frame(px=time,
                     py=c(UB, rev(LB)))
    
    p <- p + geom_polygon(aes_string(x='px', y='py'), data=df, fill="black", alpha=0.1)
  }
  
  p
}