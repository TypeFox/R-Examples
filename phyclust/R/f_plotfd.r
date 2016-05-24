### This file contains the function for fitness distirbution plot. 
### Plot the posterior distribution by each mixture component with
### sorted cluster.

plotfd <- function(X, pcobj, merge = TRUE,
    main = NULL, xlab = NULL, ylab = NULL, ...){
  N <- pcobj$N
  K <- pcobj$K

  post.class <- sort(pcobj$class.id)
  post <- pcobj$Z.normalized[order(pcobj$class.id),]
  for(k in 1:K){
    id <- which(post.class == k)
    id.new <- order(post[id, k], decreasing = TRUE)
    post[id,] <- post[id[id.new],]
  }

  xlim <- c(1, N)
  ylim <- c(0, 1)

  if(is.null(main)) main <- "Distribution Plot"
  if(is.null(xlab)) xlab <- "Sequence (sorted)"
  if(is.null(ylab)) ylab <- "Posterior"

  color <- .Color[1:K %% length(.Color) + 1]

  x <- 1:N
  if(merge){
    if(K == 1){
      plot(NULL, NULL, xlim = xlim, ylim = ylim,
           main = main, xlab = xlab, ylab = "K = 1")
      y <- post[, 1]
      lines(list(x = x, y = y), col = color[1])
    } else{
      plot(NULL, NULL, xlim = xlim, ylim = ylim,
           main = main, xlab = xlab, ylab = ylab)
      for(k in 1:K){
        lines(list(x = x, y = post[, k]), col = color[k], lwd = 3, lty = k)
      }
    }
  } else{
    v <- cumsum(pcobj$n.class)[1:(K-1)] + 0.5
    if(K == 1){
      plot(NULL, NULL, xlim = xlim, ylim = ylim,
           main = main, xlab = xlab, ylab = "K = 1")
      y <- post[, 1]
      lines(list(x = x, y = y), col = color[1])
      abline(v = v, lty = 3)
    } else{
      par(mfrow = c(K, 1))
      plot(NULL, NULL, xlim = xlim, ylim = ylim,
           main = main, xlab = "", ylab = "K = 1", axes = FALSE)
      lines(list(x = x, y = post[, 1]), col = color[1], lwd = 2)
      axis(2)
      abline(v = v, lty = 3)
      for(k in 2:(K - 1)){
        plot(NULL, NULL, xlim = xlim, ylim = ylim,
            xlab = "", ylab = paste("K =", k), axes = FALSE)
        lines(list(x = x, y = post[, k]), col = color[k], lwd = 2)
        axis(2)
        abline(v = v, lty = 3)
      }
      plot(NULL, NULL, xlim = xlim, ylim = ylim,
           main = "", xlab = xlab, ylab = paste("K =", K), axes = FALSE)
      lines(list(x = x, y = post[, K]), col = color[K], lwd = 2)
      axis(1)
      axis(2)
      abline(v = v, lty = 3)
    }
  }
} # End of plotfd().
