### This files contains functions for ploting some 2D clustering graphs
### Written: Wei-Chen Chen on 2008/10/14.

plotem <- function(emobj, x, main = NULL, xlab = NULL, ylab = NULL, ...){
  if(is.null(main)){
    main <- paste("n=", nrow(x), " K=", emobj$nclass, sep = "")
  }
  if(ncol(x) == 2){
    if(is.null(xlab)) xlab <- "x"
    if(is.null(ylab)) ylab <- "y"
    if(!is.null(emobj$class)){
      plot2d(x, emobj, main = main, xlab = xlab, ylab = ylab, ...)
    } else{
      plot2d(x, emobj, k = 1, main = main, xlab = xlab, ylab = ylab, ...)
    }
  } else{
    if(is.null(xlab)) xlab <- "Dimension"
    if(is.null(ylab)) ylab <- "Value"
    if(!is.null(emobj$class)){
      plotmd(x, class = emobj$class, main = main, xlab = xlab, ylab = ylab, ...)
    } else{
      plotmd(x, main = main, xlab = xlab, ylab = ylab, ...)
    }
  }
}

plot2d <- function(x, emobj = NULL, k = NULL, color.pch = 1,
    append.BN = TRUE, ...){
  p <- ncol(x)
  if(p != 2){
    stop("two-dimensional case only.")
  }

  plot(x, type = "n", ...)
  if(!is.null(emobj) && !is.null(k)){
    plotpt(x, color.pch = k, ...)
    if(append.BN){
      append.BN(emobj, lty = 2, col = 1, lwd = 1, cex = 1)
    }
  } else if(!is.null(emobj) && !is.null(emobj$class)){
    plotpt(x, color.pch = emobj$class, ...)
    if(append.BN){
      append.BN(emobj, lty = 2, col = 1, lwd = 1, cex = 1)
    }
  } else{
    plotpt(x, color.pch = color.pch, ...)
  }
}


plotpt <- function(x, color.pch = 1, ...){
#  color.my <- c("#00EEEEFF", "#A0EEEEFF", "#AAC6EEFF", "#EEC6BBFF",
#                "#EE66EEFF", "#DEAA00FF", "#EECF5EFF", "grey65", "grey80",
#                "#60FF00FF", "#C9FF00FF")
  # color.my <- c("#00EEEEFF", "#AAC6EEFF", "#EE66EEFF", "#EEC6BBFF",
  #               "grey65", "#60FF00FF", "#A0EEEEFF", "#EECF5EFF", "grey80",
  #               "#DEAA00FF", "#C9FF00FF")
  color.my <- color.class

  pch.my <- 0:14
  n.pch <- length(pch.my)
  n.color <- length(color.my)
  tl.symbol <- n.pch * n.color
  color.pch <- color.pch %% tl.symbol
#  pch <- pch.my[(color.pch %/% n.color + 1) %% n.pch + 1]
  pch <- pch.my[color.pch %% n.pch + 1]
  color <- color.my[color.pch %% n.color + 1]

  points(x, col = color, pch = pch, ...)
}


append.BN <- function(emobj, lty = 2, col = 1, lwd = 1, cex = 1){
  mu <- t(emobj$Mu)
  sigma <- LTSigma2variance(emobj$LTSigma)

  for(k in 1:emobj$nclass){
    plotBN(mu[,k], sigma[,,k], lty = lty, col = col, lwd = lwd, cex = cex)
  }
}


### This function are modified from "mvn2dplot()" in "mclust" package.
plotBN <- function (mu, sigma, k = 15, col = 1, lty = 1, lwd = 1, cex = 1){
    p <- length(mu)
    if (p != 2) 
        stop("two-dimensional case only")
    if (any(unique(dim(sigma)) != p)) 
        stop("mu and sigma are incompatible")
    ev <- eigen(sigma, symmetric = TRUE)
    s <- sqrt(rev(sort(ev$values)))
    V <- t(ev$vectors[, rev(order(ev$values))])
    theta <- (0:k) * (pi/(2 * k))
    x <- s[1] * cos(theta)
    y <- s[2] * sin(theta)
    xy <- cbind(c(x, -x, -x, x), c(y, y, -y, -y))
    xy <- xy %*% V
    xy <- sweep(xy, MARGIN = 2, STATS = mu, FUN = "+")

    l <- length(x)
    i <- 1:l
    for (k in 1:4) {
        lines(xy[i, ], col = col, lty = lty)
        i <- i + l
    }
    x <- s[1]
    y <- s[2]
    xy <- cbind(c(x, -x, 0, 0), c(0, 0, y, -y))
    xy <- xy %*% V
    xy <- sweep(xy, MARGIN = 2, STATS = mu, FUN = "+")
    lines(xy[1:2, ], lty = lty, col = col, lwd = lwd)
    lines(xy[3:4, ], lty = lty, col = col, lwd = lwd)
    points(mu[1], mu[2], pch = 8, col = col, cex = cex)
}


draw.colors <- function(color = NULL){
  if(is.null(color)){
    color <- colors()
  }
  n.colors <- length(color)
  tl.x <- ceiling(sqrt(n.colors))

  plot(NULL, NULL, xlim = c(0, tl.x), ylim = c(0, tl.x), axes = FALSE,
       xlab = "index", ylab = "index", main = "colors()")
  box()

  axis(1, 0:(tl.x - 1) + 0.5, 0:(tl.x - 1) * tl.x + 1)
  axis(2, 0:(tl.x - 1) + 0.5, 1:tl.x)

  rect.my <- function(i){
    xleft <- (i - 1) %/% tl.x
    ybottom <- (i - 1) %% tl.x
    rect(xleft, ybottom, xleft + 1, ybottom + 1,
         col = color[i], border = TRUE)
  }

  for(i in 1:n.colors){ rect.my(i) }
}

