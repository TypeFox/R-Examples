# Copyright (c) 2015 Santiago Barreda
# All rights reserved.

sdellipse = function (points, stdev = 1.96, density = .01, add = TRUE, 
show = TRUE, means = NULL, se = FALSE, ...){
  if (ncol (points) != 2) stop ('Points input must have exactly two columns.')
  if (!is.null(means) & nrow(points) > 2) stop ('Covariance matrix must be 2 by 2.')
  if (!is.null(means) & length(means) > 2) stop ('Exactly two means must be specified.')

  t = seq (0,2*pi+density,density)  
  x = rbind (cos(t), sin(t))
  if (is.null(means)){
    sigma = var (points)
    if (se) sigma = sigma / nrow(points)
  }
  if (!is.null(means)){
    sigma = points
    if (is.numeric(se)) sigma = sigma / se
  }
  A = eigen(sigma)$vectors %*% (diag(sqrt(eigen(sigma)$values)) * stdev)
  points = t(colMeans(points) + A%*%x)
  if (is.null(means)) points = t(colMeans(points) + A%*%x)
  if (!is.null(means)) points = t(means + A%*%x)
  
  if (add == TRUE & show == TRUE) lines (points, ...)
  if (add == FALSE & show == TRUE) plot (points, type = 'l', ...)
  invisible (points)
}

