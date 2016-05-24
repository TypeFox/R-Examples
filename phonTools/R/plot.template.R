# Copyright (c) 2015 Santiago Barreda
# All rights reserved.

plot.template = function (x, xsampa = FALSE, territorial = FALSE, ...){
  if (!territorial | is.null(x$territory)){
    sigma = x$covariance[1:2, 1:2]
    plot(x$means[, 1:2], type = "n", xlab = "F1", ylab = "F2", xlim = x$ranges[1,], ylim = x$ranges[2,])
    if (!xsampa) text(x$means[, 1:2], labels = rownames(x$means))
    if (xsampa) points(x$means[, 1:2], pch = xsampatoIPA(rownames(x$means)))
    for (i in 1:nrow(x$means)) {
      t = seq(0, 6.3, 0.05)
      xs = rbind(cos(t), sin(t))
      A = eigen(sigma)$vectors %*% (diag(sqrt(eigen(sigma)$values)) * 1.96)
      points = t(x$means[i, 1:2] + A %*% xs)
      lines(points)
    }
  }
  if (territorial & !is.null(x$territory)){
    colors = colors()[c(24, 506, 118, 610, 30, 124, 556, 258, 290, 151, 84, 657, 404)]
    plot (x$means[,1:2], type = 'n', xlim = x$ranges[1,], ylim = x$ranges[2,],xaxs='i',yaxs='i')
    for (i in 1:nrow(x$means)) polygon (x$territory[[i]][,1], x$territory[[i]][,2], col = colors[i+1], border = 1, lwd = 2)
    if (!xsampa) text(x$means[, 1:2], labels = rownames(x$means), col = 'white', cex = 2)
    if (xsampa) points(x$means[, 1:2], pch = xsampatoIPA(rownames(x$means)), col = 'white', cex = 2)
  }
  if (territorial & is.null(x$territory))
    stop ('Sorry, no territorial map exists, please use the territorialmap() function to create one.\n')
}


