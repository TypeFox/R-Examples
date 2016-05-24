gridOne <- function(Mu.1, S.1, Pi.1, grid.x, grid.y){
  p <- length(Mu.1)
  S.1 <- matrix(S.1, nrow = p, ncol = p)

  n.grid.x <- length(grid.x)
  n.grid.y <- length(grid.y)
  z <- array(rep(0, n.grid.x * n.grid.y), c(n.grid.x, n.grid.y))

  S.1.det <- det(S.1)
  S.1.inv <- solve(S.1)
  C <- Pi.1 / (2 * pi)^(p/2) / S.1.det

  for(i.x in 1:n.grid.x){
    for(i.y in 1:n.grid.y){
      X.c <- matrix(c(grid.x[i.x], grid.y[i.y]) - Mu.1, ncol = 1)
      z[i.x, i.y] <- C * exp(-0.5 * t(X.c) %*% S.1.inv %*% X.c)
    }
  }

  z
} # End of gridOne().

fill.ppcontour <- function(z,
    x = seq(0, 1, length.out = nrow(z[,, 1])),
    y = seq(0, 1, length.out = ncol(z[,, 1])),
    power = 0.4, xlab = "", ylab = "", main = ""){

  x.lim <- range(x)
  y.lim <- range(y)
  x.lim <- x.lim + c(0.02, -0.02) * (x.lim[2] - x.lim[1])
  y.lim <- y.lim + c(0.02, -0.02) * (y.lim[2] - y.lim[1])

  plot(NULL, NULL, xlab = "", ylab = "", main = "",
       xlim = x.lim, ylim = y.lim, axes = FALSE)

  for(i.k in 1:dim(z)[3]){
    image(x, y, z[,, i.k]^power,
          col = col.ppcontour$fill[[(i.k - 1) %% 7 + 1]](),
          add = TRUE)
  }

  for(i.k in 1:dim(z)[3]){
    contour(x, y, z[,, i.k]^power, nlevels = 10,
            drawlabels = FALSE,
            # col = col.ppcontour$line[[(i.k - 1) %% 7 + 1]],
            col = rgb(255, 255, 255, 50, maxColorValue = 255),
            lwd = 0.2, add = TRUE)
  }

  ### Merged contours.
  # zz <- z[,, 1]
  # for(i.k in 2:dim(z)[3]){
  #   zz <- zz + z[,, i.k]
  # }
  # contour(x, y, zz^power, nlevels = 10,
  #         drawlabels = FALSE,
  #         # col = col.ppcontour$line[[(i.k - 1) %% 7 + 1]],
  #         col = rgb(255, 255, 255, 100, maxColorValue = 255),
  #         lwd = 0.5, add = TRUE)

  box()
  # axis(1)
  # axis(2)

  invisible()
} # End of fill.ppcontour().

fill.pppoints <- function(da, class, class.true){
  for(i.k in 1:length(unique(class))){
    da.tmp <- da[class == i.k,]
    class.true.tmp <- class.true[class == i.k]
    points(da.tmp[, 1], da.tmp[, 2],
           col = col.ppcontour$point[[(i.k - 1) %% 7 + 1]],
           pch = col.ppcontour$symbol[(class.true.tmp - 1) %% 7 + 1],
           cex = 2.0)
  }
} # End of fill.pppoints().

fill.ppmu <- function(Mu){
  for(i.k in 1:nrow(Mu)){
    points(Mu[i.k, 1], Mu[i.k, 2],
           pch = "x",
           col = rgb(255, 255, 255, maxColorValue = 255),
           cex = 2.0)
  }
} # End of fill.ppmu().

plotppcontour <- function(da, Pi, Mu, S, class, class.true = NULL,
    n.grid = 128, angle = 0, xlab = "", ylab = "", main = ""){
  if(is.null(class.true)){
    class.true <- class
  }

  # Rotate plot if needed.
  rotation <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)),
                     nrow = 2, byrow = TRUE)
  K <- length(Pi)
  Mu <- Mu %*% rotation
  tmp.S <- S
  for(i.k in 1:K){
    tmp <- eigen(tmp.S[,, i.k])
    tmp$vectors <- t(t(tmp$vectors) %*% rotation)
    S[,, i.k] <- tmp$vectors %*% diag(tmp$values) %*% t(tmp$vectors)
  }
  da <- da %*% rotation

  # Add extra space for plot.
  # max.x <- max.y <- max(da)
  # min.x <- min.y <- min(da)
  min.x <- min(da[, 1])
  max.x <- max(da[, 1])
  range.x <- max.x - min.x
  min.y <- min(da[, 2])
  max.y <- max(da[, 2])
  range.y <- max.y - min.y
  bounds <- c(min.x - range.x * 0.1,
              max.x + range.x * 0.1,
              min.y - range.y * 0.1,
              max.y + range.y * 0.1)

  # Build grid points.
  grid.x <- seq(bounds[1], bounds[2], length.out = n.grid)
  grid.y <- seq(bounds[3], bounds[4], length.out = n.grid)

  K <- dim(S)[3]
  z <- array(rep(NA, n.grid * n.grid * K), c(n.grid, n.grid, K))
  ret <- NULL
  for(i.k in 1:K){
    z[,, i.k] <- gridOne(Mu[i.k,], S[,, i.k], Pi[i.k], grid.x, grid.y)
  }

  fill.ppcontour(z, x = grid.x, y = grid.y,
                 xlab = xlab, ylab = ylab, main = main)
  fill.pppoints(da, class, class.true)
  fill.ppmu(Mu)

  invisible()
} # End of plotppcontour().

