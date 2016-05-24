plot.graf <-
function(x, vars = NULL, resolution = 50, CI = 0.95, prior = FALSE, data = TRUE, jitter = 1,
			peak = FALSE, ...) {
  # find factors
  facs <- x$facs
  # if not plotting at peak, find column menas, accounting for factors
  if (!peak) {
	if(length(facs) > 0) {
		x$peak[(1:ncol(x$obsx))[-facs]] <- colMeans(x$obsx[, -facs, drop = FALSE])
		x$peak[, facs] <- sapply(data.frame(x$obsx[, facs, drop = FALSE]),
				function(x) names(sort(table(x), decreasing = TRUE))[1])
		for(i in facs) x$peak[, i] <- factor(x$peak[, i], levels = levels(x$obsx[, i]))
	} else {
		x$peak <- colMeans(x$obsx)
	}
  }
  pars <- x$peak
  k <- ncol(x$obsx)
  if (is.null(vars)) vars <- 1:k
  X <- as.data.frame(lapply(pars, rep, resolution))

  # calculate scaled data for mean function
  scaledX <- X
  if (!is.null(x$scaling)) {
    notfacs <- (1:k)
	if(length(facs) > 0) notfacs <- notfacs[-facs]
    for(i in 1:length(notfacs)) {
#      if (i %in% notfacs) {
        scaledX[, notfacs[i]] <- (scaledX[, notfacs[i]] - x$scaling[1, i]) / x$scaling[2, i]
#      }
    }
  }
  for (i in vars) {
    predX <- X
#    predscaledX <- scaledX
    if(i %in% facs) {
      # for factors
      seqx <- levels(X[, i])

      m <- length(seqx)
      predX[1:m, i] <- seqx
      predX <- predX[1:m, ]
      
#      predscaledX[1:m, i] <- seqx
#      predscaledX <- predscaledX[1:m, ]

      p <- predict(x, newdata = predX, CI = CI)
      mn <- p[, 1]
      up <- p[, 3]
      low <- p[, 2]
      width = 0.2
      
      x.plot <- as.numeric(predX[, i])

      nam <- ifelse(is.null (colnames(x$x)[i]), "covariate", colnames(x$x)[i]) 
      # set up plot
      plot(mn ~ predX[, i], ylim = c(0, 1), type = 'p', xlab = nam, border ='white',
           ylab = 'probability of presence', ...)

      # plot CIs
      rect(xleft = x.plot - width, ybottom = low, xright = x.plot + width, ytop = up, col = 'light grey',
           border = NA)
      # plot mean
      for(j in 1:m) {
        lines(x = c(x.plot[j] - width, x.plot[j] + width), y = rep(mn[j], 2), col = rgb(0.4, 0.4, 0.4), lwd = 2)
      }

      if (prior) {
        xs <- sort(x.plot)
        xs[c(1, m)] <- c(-1, m + 1)        
        lines(x$mnfun(predX) ~ xs, lty = 2, col = rgb(0.4, 0.4, 0.4))
      }

      if (data) {
        rug(jitter(x$x[x$obsy == 1, i], jitter), side = 3, col = 'dark grey')
        rug(jitter(x$x[x$obsy == 0, i], jitter), side = 1, col = 'dark grey')
        box()
      }
      
    } else {
      # for continuous variables
      seqx <- seq(min(x$obsx[, i]), max(x$obsx[, i]), length.out = resolution)
  
      predX[, i] <- seqx
      
 #     predscaledX[, i] <- (seqx - x$scaling[1, i]) / x$scaling[2, i]

      
      p <- predict(x, newdata = predX, CI = CI)
      mn <- p[, 1]
      up <- p[, 3]
      low <- p[, 2]
      nam <- ifelse(is.null (colnames(x$x)[i]), "covariate", colnames(x$x)[i]) 
      plot(mn ~ seqx, ylim = c(0, 1), type = 'n', xlab = nam,
           ylab = 'probability of presence', ...)
      polygon(x = c(seqx, rev(seqx)), y = c(up, rev(low)), col = 'light grey',
              border = NA)
      lines(mn ~ seqx, col = rgb(0.4, 0.4, 0.4), lwd = 2)
      if (prior) {
        lines(x$mnfun(predX) ~ seqx, lty = 2, col = 'dark grey')
      }
      if (data) {
        rug(jitter(x$obsx[x$obsy == 1, i], jitter), side = 3, col = 'dark grey')
        rug(jitter(x$obsx[x$obsy == 0, i], jitter), side = 1, col = 'dark grey')
        box()
      }
    }
  }
}