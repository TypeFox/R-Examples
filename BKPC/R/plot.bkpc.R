plot.bkpc <-
function(x, type = "default", n.burnin = 0, ...){
  
  n.samples <- dim(x$beta)[1]
  if (n.burnin >= n.samples) stop("error: too many burn-in iterations specified")
  if (n.burnin < 0) n.burnin <- 0
  
  if (type == "tracePlot"){
    matplot(x$beta[(n.burnin + 1) : n.samples, ], main = expression(beta), ylab = "", xlab = "Iterations", type='l', ...)
    matplot(x$tau[(n.burnin + 1) : n.samples, ],  main = expression(tau),  ylab = "", xlab = "Iterations", type='l', ...)
    matplot(x$z[(n.burnin + 1) : n.samples, ],  main = "z", xlab = "Iterations", ylab = "", type='l', ...)
    plot(x$sigmasq[(n.burnin + 1) : n.samples,1,  drop = FALSE],  main = expression(sigma^2), ylab = "", xlab = "Iterations", type='l', ...)
  
      }
  
  else if (type == "default"){
  plot(apply(x$beta[(n.burnin + 1) : n.samples, ], 2, median), pch  = 20, xlab = "", ylab = expression(beta),
       ylim = c(min(x$beta[(n.burnin + 1) : n.samples, ]), max(x$beta[(n.burnin + 1) : n.samples, ])), col = 2)
  
  u <- apply(x$beta[(n.burnin + 1) : n.samples, ], 2, quantile, probs = c(0.9))
  l <- apply(x$beta[(n.burnin + 1) : n.samples, ], 2, quantile, probs = c(0.1))
  for(i in 1 : dim(x$beta)[2])lines(c(i, i), c(l[i], u[i]), ...)
  abline(h = 0)
  
  u <- apply(x$tau[(n.burnin + 1) : n.samples, ], 2, quantile, probs = c(0.9))
  l <- apply(x$tau[(n.burnin + 1) : n.samples, ], 2, quantile, probs = c(0.1))
  plot(apply(x$tau[(n.burnin + 1) : n.samples, ], 2, median), pch  = 20, xlab = "", ylab = expression(tau),
       ylim = c(min(x$tau[(n.burnin + 1) : n.samples, ]), max(x$tau[(n.burnin + 1) : n.samples, ])), col = 2)
  for(i in 1 : dim(x$tau)[2])lines(c(i, i),  c(l[i], u[i]), ...)
  
  
  u <- apply(x$z[(n.burnin + 1) : n.samples, ], 2, quantile, probs = c(0.9))
  l <- apply(x$z[(n.burnin + 1) : n.samples, ], 2, quantile, probs = c(0.1))
  plot(apply(x$z[(n.burnin + 1) : n.samples, ], 2, median), pch  = 20, xlab = "", ylab = expression(z),
       ylim = c(min(x$z[(n.burnin + 1) : n.samples, ]), max(x$z[(n.burnin + 1) : n.samples, ])), col = 2)
  for(i in 1 : dim(x$z)[2])lines(c(i, i), c(l[i], u[i]), ...)
  
  plot(apply(x$sigmasq[(n.burnin + 1) : n.samples,1,  drop = FALSE], 2, median), pch  = 20, xlab = "", ylab = expression(sigma^2),
       ylim = c(min(x$sigmasq[(n.burnin + 1) : n.samples,1,  drop = FALSE]), max(x$sigmasq[(n.burnin + 1) : n.samples, 1,  drop = FALSE])), col = 2)
  lines(c(1, 1), c(apply(x$sigmasq[(n.burnin + 1) : n.samples, 1,  drop = FALSE], 2, quantile, probs = c(0.1)),
                   apply(x$sigmasq[(n.burnin + 1) : n.samples, 1,  drop = FALSE], 2, quantile, probs = c(0.9))), ...)
  
  }
  else if (type == "boxPlot"){
    boxplot(x$beta[(n.burnin + 1) : n.samples, ], ylab = expression(beta), ...)
    boxplot(x$tau[(n.burnin + 1) : n.samples, ], main = expression(tau), ...)
    boxplot(x$z[(n.burnin + 1) : n.samples, ], ylab = "z", ...)
    boxplot(x$sigmasq[(n.burnin + 1) : n.samples, 1,  drop = FALSE], ylab = expression(sigma^2), ...)
    }
    else  stop("error: Plot type not supported for a bkpc object")
}
