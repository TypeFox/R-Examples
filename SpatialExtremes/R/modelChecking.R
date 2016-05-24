qqextcoeff <- function(fitted, estim = "ST", marge = "emp",
                       xlab = "Semi-Empirical", ylab = "Model", ...){

  if (!any("maxstab" %in% class(fitted)))
    stop("This functin is only available for 'maxstab' objects")

  data <- fitted$data
  coord <- fitted$coord
  ext.coeff <- fitted$ext.coeff

  if (fitted$iso){
    dist <- distance(coord)
    exco.mod <- ext.coeff(dist)
  }

  else {
    dist <- distance(coord, vec = TRUE)
    exco.mod <- apply(dist, 1, ext.coeff)
  }

  exco.emp <- fitextcoeff(data, coord, plot = FALSE, estim = estim,
                          marge = marge)$ext.coeff[,"ext.coeff"]

  plot(exco.emp, exco.mod, xlab = xlab, ylab = ylab, ...)
  abline(0, 1)
}

qqgev <- function(fitted, xlab, ylab, ...){

  data <- fitted$data
  n.site <- ncol(data)

  gev.param <- t(apply(data, 2, gevmle))
  pred <- predict(fitted)

  if (missing(xlab))
    xlab <- c(expression(mu[MLE]), expression(sigma[MLE]), expression(xi[MLE]))

  if (missing(ylab))
    ylab <- c(expression(mu[Model]), expression(sigma[Model]), expression(xi[Model]))

  op <- par(mfrow=c(1,3))
  on.exit(par(op))

  if (length(unique(pred[,"loc"])) != 1){
    xlim <- ylim <- range(c(gev.param[,"loc"], pred[,"loc"]))
    plot(gev.param[,"loc"], pred[,"loc"], xlab = xlab[1], ylab = ylab[1], ...,
         xlim = xlim, ylim = ylim)
    abline(0, 1)
  }

  else{
    hist(gev.param[,"loc"], xlab = xlab[1], main = "")
    axis(3, at = pred[1, "loc"], labels = ylab[1])
  }

  if (length(unique(pred[,"scale"])) !=1){
    xlim <- ylim <- range(c(gev.param[,"scale"], pred[,"scale"]))
    plot(gev.param[,"scale"], pred[,"scale"], xlab = xlab[2], ylab = ylab[2],
         ..., xlim = xlim, ylim = ylim)
    abline(0, 1)
  }

  else{
    hist(gev.param[,"scale"], xlab = xlab[2], main = "")
    axis(3, at = pred[1, "scale"], labels = ylab[2])
  }

  if (length(unique(pred[,"shape"])) != 1){
    xlim <- ylim <- range(c(gev.param[,"shape"], pred[,"shape"]))
    plot(gev.param[,"shape"], pred[,"shape"], xlab = xlab[3], ylab = ylab[3],
         ..., xlim = xlim, ylim = ylim)
    abline(0, 1)
  }

  else{
    hist(gev.param[,"shape"], xlab = xlab[3], main = "")
    axis(3, at = pred[1, "shape"], labels = ylab[3])
  }
}

plot.copula <- function(x, ..., sites){
  n.site <- ncol(x$data)
  n.obs <- nrow(x$data)

  ##The graph is as follows :
  ## The diagonal are return level plots
  ## The upper part is the extremal coeff. function
  ## The lower one is qq plots on max for pairs
  ## The two remaining plots are hist. of pairwise distance and location map
  if (missing(sites))
    sites <- sample(1:n.site, 4)

  else if (length(sites) != 4)
    stop("'sites' must have length 4")

  op <- par(no.readonly = TRUE)
  layout(matrix(c(1,6,7,9,13,2,8,10,5,5,3,11,5,5,12,4), 4))
  par(mar = c(4,4,1,0.5))
  on.exit(par(op))

  ## Return level plots
  ##covariates <- cbind(x$coord[sites,], x$marg.cov[sites,,drop=FALSE])
  ##gev.param <- predict(x, covariates, std.err = FALSE)[,c("loc", "scale", "shape")]
  gev.param <- predict(x, std.err = FALSE)[sites, c("loc", "scale", "shape")]

  for (i in 1:4){
    boot <- matrix(NA, nrow = 1000, ncol = n.obs)
    loc <- gev.param[i,1]
    scale <- gev.param[i,2]
    shape <- gev.param[i,3]
    probs <- 1:n.obs / (n.obs + 1)

    for (j in 1:1000)
      boot[j,] <- sort(rgev(n.obs, loc, scale, shape))

    ci <- apply(boot, 2, quantile, c(0.025, 0.975))
    matplot(1 / (1 - probs), t(ci), pch ="-", col = 1,
            xlab = "Return Period", ylab = "Return level", log = "x")
    fun <- function(T) qgev(1 - 1/T, loc, scale, shape)
    curve(fun, from = 1.001, to = 100, add = TRUE)
    points(1 / (1 - probs), sort(x$data[,sites[i]]))

  }

  ##F-madogram
  fmadogram(x$data, x$coord, which = "ext", col = "lightgrey")
  fmadogram(fitted = x, which = "ext", add = TRUE, n.bins = n.site)

  ##Pairwise maxima on the Gumbel scale
  model <- x$model

  DoF <- x$par["DoF"]
  nugget <- x$par["nugget"]
  range <- x$par["range"]
  smooth <- x$par["smooth"]
  sim.copula <- rcopula(n.obs * 1000, x$coord[sites,], x$copula, x$cov.mod,
                        nugget = nugget, range = range, smooth = smooth, DoF = DoF)

  sim.copula <- array(log(sim.copula), c(n.obs, 1000, 4))

  gumb <- log(apply(x$data[,sites], 2, gev2frech, emp = TRUE))
  ##Plot of the pairwise maxima
  for (i in 1:3){
    for (j in (i+1):4){
      pair.max <- sort(apply(gumb[,c(i, j)], 1, max))
      sim.pair.max <- apply(pmax(sim.copula[,,i], sim.copula[,,j]), 2, sort)
      dummy <- rowMeans(sim.pair.max)
      ci <- apply(sim.pair.max, 1, quantile, c(0.025, 0.975))
      matplot(dummy, t(ci), pch = "-", col = 1, , xlab = "Model",
              ylab = "Observed")
      points(dummy, pair.max)
      abline(0, 1)
      h <- distance(x$coord[sites[c(i,j)],])
      legend("bottomright", paste("h =", round(h, 2)), bty = "n")
    }
  }

  ##Plot of the blockwise maxima
  block.max <- sort(apply(gumb, 1, max))
  sim.block.max <- sim.copula[,,1]
  for (i in 2:4)
    sim.block.max <- pmax(sim.block.max, sim.copula[,,i])

  sim.block.max <- apply(sim.block.max, 2, sort)
  dummy <- rowMeans(sim.block.max)
  ci <- apply(sim.block.max, 1, quantile, c(0.025, 0.975))
  matplot(dummy, t(ci), pch = "-", col = 1, xlab = "Model", ylab = "Observed")
  points(dummy, block.max)
  abline(0, 1)

  ##Plot of the locations - identifying the one selected
  plot(x$coord, type = "n")
  points(x$coord[-sites,])
  points(x$coord[sites,], pch = c("1", "2", "3", "4"), col = "blue")

}

plot.maxstab <- function(x, ..., sites){
  n.site <- ncol(x$data)

  ##The graph is as follows :
  ## The diagonal are return level plots
  ## The upper part is the extremal coeff. function
  ## The lower one is qq plots on max for pairs
  ## The two remaining plots are hist. of pairwise distance and location map
  if (missing(sites))
    sites <- sample(1:n.site, 4)

  else if (length(sites) != 4)
    stop("'sites' must have length 4")

  op <- par(no.readonly = TRUE)
  layout(matrix(c(1,6,7,9,13,2,8,10,5,5,3,11,5,5,12,4), 4))
  par(mar = c(4,4,1,0.5))
  on.exit(par(op))

  ## Return level plots
  ##covariates <- cbind(x$coord[sites,], x$marg.cov[sites,,drop=FALSE])
  ##gev.param <- predict(x, covariates, std.err = FALSE)[,c("loc", "scale", "shape")]
  gev.param <- predict(x, std.err = FALSE)[sites, c("loc", "scale", "shape")]

  for (i in 1:4){
    n.obs <- length(na.omit(x$data[,sites[i]]))## usefull since
                                               ## missing values are
                                               ## now allowed
    boot <- matrix(NA, nrow = 1000, ncol = n.obs)
    loc <- gev.param[i,1]
    scale <- gev.param[i,2]
    shape <- gev.param[i,3]
    probs <- 1:n.obs / (n.obs + 1)

    for (j in 1:1000)
      boot[j,] <- sort(rgev(n.obs, loc, scale, shape))

    ci <- apply(boot, 2, quantile, c(0.025, 0.975))
    matplot(1 / (1 - probs), t(ci), pch ="-", col = 1,
            xlab = "Return Period", ylab = "Return level", log = "x")
    fun <- function(T) qgev(1 - 1/T, loc, scale, shape)
    curve(fun, from = 1.001, to = 100, add = TRUE)
    points(1 / (1 - probs), sort(x$data[,sites[i]]))

  }

  ##F-madogram
  fmadogram(x$data, x$coord, which = "ext", col = "lightgrey")
  fmadogram(fitted = x, which = "ext", add = TRUE, n.bins = n.site)

  ##Pairwise maxima on the Gumbel scale
  n.obs <- nrow(x$data[,sites])
  model <- x$model
  notimplemented <- FALSE
  if (model == "Smith"){
    cov11 <- x$par["cov11"]
    cov12 <- x$par["cov12"]
    cov22 <- x$par["cov22"]
    sim.maxstab <- rmaxstab(n.obs * 1000, x$coord[sites,], "gauss",
                            cov11 = cov11, cov12 = cov12, cov22 = cov22)
  }

  else if (model == "Schlather"){
    nugget <- x$par["nugget"]
    range <- x$par["range"]
    smooth <- x$par["smooth"]
    sim.maxstab <- rmaxstab(n.obs * 1000, x$coord[sites,], x$cov.mod,
                            nugget = nugget, range = range, smooth = smooth)
  }

  else if (model == "Geometric"){
    sigma2 <- x$par["sigma2"]
    nugget <- x$par["nugget"]
    range <- x$par["range"]
    smooth <- x$par["smooth"]
    cov.mod <- paste("g", x$cov.mod, sep = "")
    sim.maxstab <- rmaxstab(n.obs * 1000, x$coord[sites,], cov.mod,
                            sigma2 = sigma2, nugget = nugget, range = range,
                            smooth = smooth)
  }

  else if (model == "Brown-Resnick"){
    range <- x$par["range"]
    smooth <- x$par["smooth"]
    sim.maxstab <- rmaxstab(n.obs * 1000, x$coord[sites,], "brown",
                            range = range, smooth = smooth)
  }

  else if (model == "Extremal-t"){
    DoF <- x$par["DoF"]
    nugget <- x$par["nugget"]
    range <- x$par["range"]
    smooth <- x$par["smooth"]
    cov.mod <- paste("t", x$cov.mod, sep = "")
    sim.maxstab <- rmaxstab(n.obs * 1000, x$coord[sites,], cov.mod,
                            DoF = DoF, nugget = nugget, range = range,
                            smooth = smooth)
  }

  if (notimplemented){
    for (i in 1:7){
      plot(0, 0, type = "n", bty = "n", axes = FALSE, xlab = "", ylab = "")
      text(0, 0, "Not implemented")
    }
  }

  else {
    sim.maxstab <- array(log(sim.maxstab), c(n.obs, 1000, 4))

    for (i in 1:4){
        idx.na <- which(is.na(x$data[,sites[i]]))
        sim.maxstab[idx.na,,i] <- NA
    }

    gumb <- log(apply(x$data[,sites], 2, gev2frech, emp = TRUE))
    ##Plot of the pairwise maxima
    for (i in 1:3){
      for (j in (i+1):4){
        pair.max <- sort(apply(gumb[,c(i, j)], 1, max))
        sim.pair.max <- apply(pmax(sim.maxstab[,,i], sim.maxstab[,,j]), 2, sort)
        dummy <- rowMeans(sim.pair.max)
        ci <- apply(sim.pair.max, 1, quantile, c(0.025, 0.975))
        matplot(dummy, t(ci), pch = "-", col = 1, , xlab = "Model",
                ylab = "Observed")
        points(dummy, pair.max)
        abline(0, 1)
        h <- distance(x$coord[sites[c(i,j)],])
        legend("bottomright", paste("h =", round(h, 2)), bty = "n")
      }
    }

    ##Plot of the blockwise maxima
    block.max <- sort(apply(gumb, 1, max, na.rm = TRUE))
    ##block.max[!is.finite(block.max)] <- NA
    sim.block.max <- sim.maxstab[,,1]
    for (i in 2:4)
      sim.block.max <- pmax(sim.block.max, sim.maxstab[,,i], na.rm = TRUE)

    sim.block.max <- apply(sim.block.max, 2, sort, na.last = NA)
    dummy <- rowMeans(sim.block.max)
    ci <- apply(sim.block.max, 1, quantile, c(0.025, 0.975))
    matplot(dummy, t(ci), pch = "-", col = 1, xlab = "Model", ylab = "Observed")
    points(dummy, block.max)
    abline(0, 1)
  }

  ##Plot of the locations - identifying the one selected
  plot(x$coord, type = "n")
  points(x$coord[-sites,])
  points(x$coord[sites,], pch = c("1", "2", "3", "4"), col = "blue")

}



