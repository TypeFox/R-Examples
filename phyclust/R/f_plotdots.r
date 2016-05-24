### This file contains the function for the dots plot.
### Plot the difference for all sequences based on the first sequence.

plotdots <- function(X, X.class = NULL, Mu = NULL, code.type = .code.type[1],
    diff.only = TRUE, fill = FALSE, label = TRUE,
    with.gap = FALSE, xlim = NULL, ylim = NULL,
    main = "Dots Plot", xlab = "Sites", ylab = "Sequences",
    missing.col = "gray95", ...){
  if(sum(code.type %in% .code.type) != 1){
    stop("The code.type is not found.")
  }

  if(! is.null(X.class)){
    if(length(X.class) != nrow(X)){
      stop("length(X.class) != N.")
    }

    id.X <- order(X.class)
    X <- X[id.X,]
    K <- max(X.class)
    X.nc.cum <- rep(0, K)
    for(k in 1:K) X.nc.cum[k] <- sum(X.class %in% 1:k)
  }

  if(is.null(Mu)){
    Mu <- find.consensus(X, code.type = code.type, with.gap = with.gap)
  } else{
    if(length(Mu) != ncol(X)){
      stop("length(Mu) != L.")
    }
  }

  if(diff.only){
    tl.diff <- apply(X, 2, function(x) length(unique(x[x != .missing.code$mid])))
    X <- X[, tl.diff > 1]
    Mu <- Mu[tl.diff > 1]
  }

  N <- nrow(X)
  L <- ncol(X)

  if(is.null(xlim)){
    xlim <- c(1, L + 1)
#  } else{
#    if(xlim[1] < 1){
#      xlim[1] <- 1
#    }
  }
  if(is.null(ylim)){
    if(label){
      ylim <- c(N + 3, -3)
    } else{
      ylim <- c(N + 1, -1)
    }
  }

  if(code.type == .code.type[1]){
    my.col <- c("green3", "blue2", "#CC00CC", "red2", "gray", missing.col)
  } else if(code.type == .code.type[2]){
    my.col <- c("green3", "blue2", "gray", missing.col)
  } else{
    stop("code.type is not implemented.")
  }


  X <- X + 1
  Mu <- Mu + 1
  missing.code.mid <- .missing.code$mid + 1
  plot(NULL, NULL, type = "n", xlim = xlim, ylim = ylim,
       main = main, xlab = xlab, ylab = ylab)

  ### Consensus sequence.
  y.top <- -3
  for(j in 1:L){
    tmp.Mu <- ifelse(Mu[j] != missing.code.mid, Mu[j], length(my.col))
    x.left <- j
    rect(x.left, y.top + 4, x.left + 1, y.top,
         col = my.col[tmp.Mu],
         border = NA)
  }
  abline(h = y.top + 4, lty = 3, lwd = 0.5)

  ### Different sites.
  plot.column.fill <- function(j){
    x.left <- j
    x.right <- j + 1
    tmp.X <- X[, j]
    tmp.X[tmp.X == missing.code.mid] <- length(my.col)
    for(i in 1:N){
      rect(x.left, i + 1, x.right, i, col = my.col[tmp.X[i]], border = NA)
    }
  }
  plot.column <- function(j){
    x.left <- j
    x.right <- j + 1
    tmp.X <- X[, j]
    tmp.X[tmp.X == missing.code.mid] <- length(my.col)
    for(i in which(tmp.X != Mu[j])){
      rect(x.left, i + 1, x.right, i, col = my.col[tmp.X[i]], border = NA)
    }
  }

  flag.diff <- rowSums(t(X) != Mu) > 0
  if(xlim[2] > L){
    xlim[2] <- L
  } else if(xlim[1] < 1){
    xlim[1] <- 1
  }
  J <- xlim[1]:xlim[2]
  if(diff.only){
    J <- J[flag.diff[J]]
  }
  if(fill){
    lapply(J, plot.column.fill)
  } else{
    lapply(J, plot.column)
  }

  ### Segregating sites.
  if(label){
    y.top <- N + 1
    for(j in 1:L){
      x.left <- j
      if(flag.diff[j]){
        rect(x.left, y.top + 4, x.left + 1, y.top, col = "orange", border = NA)
      }
    }
    abline(h = y.top, lty = 3, lwd = 0.5)
  }

  if(! is.null(X.class)) abline(h = X.nc.cum + 1, lty = 3, lwd = 0.5)
} # End of plotdots().

