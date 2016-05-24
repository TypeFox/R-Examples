### This file contains the functions for histograms.

plothist <- function(X, X.class = NULL, Mu = NULL, fill.color = .Color,
    draw.all = TRUE, main = "Mutation counts",
    xlab = "Difference", ylab = "Counts", append = FALSE){
  X.tmp <- t(X)

  if(is.null(Mu)){
    if(is.null(X.class)){
      common <- X.tmp[, 1]
    } else{
      common <- X.tmp[, X.class == 1][, 1]
    }
  } else{
    common <- Mu[1,]
  }

  if(is.null(X.class)){
    color <- .Color[1]

    X.table <- table(colSums(X.tmp != common))
    X.ht <- cbind(as.numeric(names(X.table)), as.numeric(X.table))
    plothist.my(X.ht, main = main, xlab = xlab, ylab = ylab, fill.color = color)
  } else{
    K <- length(unique(X.class))
    color <- .Color[0:(K - 1) %% length(.Color) + 1]

    X.ht <- NULL
    xlim <- NULL
    ylim <- NULL
    for(k in 1:K){
      X.k.tmp <- X.tmp[, X.class == k]
      X.table <- table(colSums(X.k.tmp != common))
      X.ht[[k]] <- cbind(as.numeric(names(X.table)), as.numeric(X.table))
      xlim <- c(xlim, X.ht[[k]][, 1])
      ylim <- c(ylim, X.ht[[k]][, 2])
    }
    xlim <- range(xlim) + c(-0.6, 0.6)
    ylim <- c(0, max(ylim) + 0.1)

    if(draw.all){
      if(! append){
        par(mfrow = c(K + 1, 1), mar = c(2, 4, 0.5, 2))
      }
      X.table <- table(colSums(X.tmp != common))
      X.ht.all <- cbind(as.numeric(names(X.table)), as.numeric(X.table))
      ylim <- c(0, max(X.ht.all[, 2]) + 0.1)
      plothist.my(X.ht.all, xlim = xlim, ylim = ylim,
                 main = "", xlab = "", ylab = "")
    } else{
      if(! append){
        par(mfrow = c(K, 1), mar = c(2, 4, 0.5, 2))
      }
    }

    for(k in 1:K){
      plothist.my(X.ht[[k]], xlim = xlim, ylim = ylim, fill.color = color[k],
                 main = "", xlab = "", ylab = "")
    }
  }
} # End of plothist().


plothist.my <- function(X.ht, xlim = NULL, ylim = NULL, fill.color = NULL,
  main = NULL, xlab = NULL, ylab = NULL){
  if(is.null(xlim)) xlim <- range(X.ht[, 1]) + c(-0.6, 0.6)
  if(is.null(ylim)) ylim <- c(0, max(X.ht[, 2]) + 0.1)

  plot(NULL, NULL, type = "n", xlim = xlim, ylim = ylim, axes = FALSE,
       main = main, ylab = ylab, xlab = xlab)

  for(i in 1:nrow(X.ht)){
    x.left <- X.ht[i, 1] - 0.5
    x.right <- X.ht[i, 1] + 0.5
    y.top <- X.ht[i, 2]
    y.bottom <- 0
    rect(x.left, y.bottom, x.right, y.top, col = fill.color, border = 1)
  }

  axis(1)
  axis(2)
} # End of plothist.my().

