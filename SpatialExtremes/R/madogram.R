madogram <- function(data, coord, fitted, n.bins, gev.param = c(0, 1, 0),
                     which = c("mado", "ext"), xlab, ylab, col = c(1, 2),
                     angles = NULL, marge = "emp", add = FALSE, xlim = c(0, max(dist)), ...){

  if (missing(fitted) && missing(data) && missing(coord))
    stop("You must either specify a fitted model OR 'data' and 'coord'")

  fit.curves <- FALSE

  if (missing(fitted) & (marge == "model"))
    stop("'marge' can be set to 'model' only if you supplied 'fitted'")

  if (!missing(fitted)){
    data <- fitted$data
    coord <- fitted$coord

    if ((fitted$model == "Smith") && !fitted$iso)
      warning("The madogram is valid only for isotropic models. The fitted curves won't be plotted.")

    else{
      fit.curves <- TRUE
      ext.coeff.fit <- fitted$ext.coeff
    }
  }

  if (is.null(dim(coord))){
    if (length(coord) != ncol(data))
      stop("'data' and 'coord' don't match")
  }

  else if (nrow(coord) != ncol(data))
    stop("'data' and 'coord' don't match")

  if (!is.null(angles) & !missing(n.bins))
    stop("'It is not possible to pass 'n.bins' and 'angles' at the same time")

  if (any(!(which %in% c("mado", "ext"))))
    stop("'which' must be either 'mado', 'ext' or both")

  if (!(marge %in% c("mle", "emp", "model")))
    stop("'marge' must be either 'mle', 'emp' or 'model'")

  n.site <- ncol(data)
  n.obs <- nrow(data)
  n.pairs <- (n.site - 1) * n.site / 2
  dist <- distance(coord)

  if (!is.null(angles)){
    distVec <- distance(coord, vec = TRUE)
    n.angles <- length(angles)
    angles.coord <- atan2(distVec[,2], distVec[,1])

    col.angles <- rep(NA, n.site * (n.site - 1) / 2)
    idx.angles <- list()
    for (i in 2:n.angles){
      idx <- which((angles.coord < angles[i]) & (angles.coord >= angles[i-1]))
      idx.angles <- c(idx.angles, list(idx))
      col.angles[idx] <- i-1
    }
  }

  if (marge == "emp")
    data <- t(t(apply(data, 2, rank, na.last = "keep")) / (colSums(is.finite(data)) + 1))

  else if (marge == "mle"){
    for (i in 1:n.site){
      param <- gevmle(data[,i])
      data[,i] <- pgev(data[,i], param[1], param[2], param[3])
    }
  }

  else{
    param <- predict(fitted, std.err = FALSE)
    for (i in 1:n.site)
      data[,i] <- pgev(data[,i], param[i,"loc"], param[i,"scale"], param[i,"shape"])
  }

  data <- qgev(data, gev.param[1], gev.param[2], gev.param[3])

  mado <- .C("madogram", as.double(data), as.integer(n.obs),
             as.integer(n.site), mado = double(n.pairs),
             PACKAGE = "SpatialExtremes", NAOK = TRUE)$mado

  if (!missing(n.bins)){
    bins <- c(0, quantile(dist, 1:n.bins/(n.bins + 1)), max(dist))
    madoBinned <- rep(NA, length = n.bins + 1)

    for (k in 1:(n.bins + 1)){
      idx <- which((dist <= bins[k+1]) & (dist > bins[k]))

      if (length(idx)>0)
        madoBinned[k] <- mean(mado[idx])
    }

    mado <- madoBinned
    dist <- (bins[-1] + bins[-(n.bins+2)])/2
  }

  if (gev.param[3] == 0)
    ext.coeff <- exp(mado/gev.param[2])

  else
    ext.coeff <- gev2frech(gev.param[1] + mado / gamma(1 - gev.param[3]),
                           gev.param[1], gev.param[2], gev.param[3])

  if (length(which) == 2){
    op <- par(mfrow=c(1,2))
    on.exit(par(op))
  }

  if (missing(xlab))
    xlab <- "h"

  if (missing(ylab))
    ylab <- c(expression(nu(h)), expression(theta(h)))

  if (add){
    if (any(which == "mado")){

        if (!is.null(angles))
            points(dist, mado, col = col.angles, ...)

        else
            points(dist, mado, col = col[1], ...)

       if (fit.curves){
         if (gev.param[3] == 0)
           mado.fit <- function(h)
             gev.param[2] * log(ext.coeff.fit(h))
         else
           mado.fit <- function(h)
             gev.param[2] * gamma(1 - gev.param[3]) *
               (ext.coeff.fit(h)^gev.param[3] - 1) / gev.param[3]

         curve(mado.fit, from = xlim[1], to = xlim[2], add = TRUE, col = col[2], ...)
       }
    }

    if (any(which == "ext")){

        if (!is.null(angles))
            points(dist, ext.coeff, col = col.angles, ...)

        else
            points(dist, ext.coeff, col = col[1], ...)

      if (fit.curves)
        curve(ext.coeff.fit, from = xlim[1], to = xlim[2], add = TRUE, col = col[2], ...)
    }
  }

  else{
    if (any(which == "mado")){

        if (!is.null(angles))
            plot(dist, mado, xlab = xlab, ylab = ylab[1], col = col.angles, xlim = xlim, ...)

        else
            plot(dist, mado, xlab = xlab, ylab = ylab[1], col = col[1], xlim = xlim, ...)

      if (!missing(fitted)){
        if (gev.param[3] == 0)
          mado.fit <- function(h)
            gev.param[2] * log(ext.coeff.fit(h))

        else
          mado.fit <- function(h)
            gev.param[2] * gamma(1 - gev.param[3]) *
              (ext.coeff.fit(h)^gev.param[3] - 1) / gev.param[3]

        curve(mado.fit, from = xlim[1], to = xlim[2], add = TRUE, col = col[2], ...)
      }
    }

    if (any(which == "ext")){

        if (!is.null(angles))
            plot(dist, ext.coeff, xlab = xlab, ylab = ylab[2], col = col.angles,
                 xlim = xlim, ylim = c(1, max(2, ext.coeff, na.rm = TRUE)), ...)

        else
            plot(dist, ext.coeff, xlab = xlab, ylab = ylab[2], col = col[1],
                 xlim = xlim, ylim = c(1, max(2, ext.coeff, na.rm = TRUE)), ...)

      if (fit.curves)
        curve(ext.coeff.fit, from = xlim[1], to = xlim[2], add = TRUE, col = col[2], ...)
    }
  }

  invisible(cbind(dist = dist, madogram = mado, ext.coeff = ext.coeff))
}

fmadogram <- function(data, coord, fitted, n.bins, which = c("mado", "ext"),
                      xlab, ylab, col = c(1, 2), angles = NULL, marge = "emp",
                      add = FALSE, xlim = c(0, max(dist)), ...){

  if (missing(fitted) && missing(data) && missing(coord))
    stop("You must either specify a fitted model OR 'data' and 'coord'")

  fit.curves <- FALSE

  if (missing(fitted) & (marge == "model"))
    stop("'marge' can be set to 'model' only if you supplied 'fitted'")

  if (!missing(fitted)){
    data <- fitted$data
    coord <- fitted$coord

    if ((fitted$model == "Smith") && !fitted$iso)
      warning("The madogram is valid only for isotropic models. The fitted curves won't be plotted.")

    else{
      fit.curves <- TRUE
      ext.coeff.fit <- fitted$ext.coeff
    }
  }

  if (is.null(dim(coord))){
    if (length(coord) != ncol(data))
      stop("'data' and 'coord' don't match")
  }

  else if (nrow(coord) != ncol(data))
    stop("'data' and 'coord' don't match")

  if (!is.null(angles) & !missing(n.bins))
    stop("'It is not possible to pass 'n.bins' and 'angles' at the same time")

  if (any(!(which %in% c("mado", "ext"))))
    stop("'which' must be either 'mado', 'ext' or both")

  if (!(marge %in% c("mle", "emp", "model")))
    stop("'marge' must be either 'mle', 'emp' or 'model'")

  n.site <- ncol(data)
  n.obs <- nrow(data)
  n.pairs <- (n.site - 1) * n.site / 2
  dist <- distance(coord)

  if (!is.null(angles)){
    distVec <- distance(coord, vec = TRUE)
    n.angles <- length(angles)
    angles.coord <- atan2(distVec[,2], distVec[,1])

    col.angles <- rep(NA, n.site * (n.site - 1) / 2)
    idx.angles <- list()
    for (i in 2:n.angles){
      idx <- which((angles.coord < angles[i]) & (angles.coord >= angles[i-1]))
      idx.angles <- c(idx.angles, list(idx))
      col.angles[idx] <- i-1
    }
  }

  if (marge == "emp")
      data <- t(t(apply(data, 2, rank, na.last = "keep")) / (colSums(is.finite(data)) + 1))

  else if (marge == "mle"){
    for (i in 1:n.site){
      param <- gevmle(data[,i])
      data[,i] <- pgev(data[,i], param["loc"], param["scale"], param["shape"])
    }
  }

  else{
    param <- predict(fitted, std.err = FALSE)
    for (i in 1:n.site)
      data[,i] <- pgev(data[,i], param[i,"loc"], param[i,"scale"], param[i,"shape"])
  }

  fmado <- .C("madogram", as.double(data), as.integer(n.obs),
              as.integer(n.site), mado = double(n.pairs),
              PACKAGE = "SpatialExtremes", NAOK = TRUE)$mado

  if (!missing(n.bins)){
    bins <- c(0, quantile(dist, 1:n.bins/(n.bins + 1)), max(dist))
    fmadoBinned <- rep(NA, length = n.bins + 1)

    for (k in 1:(n.bins+1)){
      idx <- which((dist <= bins[k+1]) & (dist > bins[k]))

      if (length(idx)>0)
        fmadoBinned[k] <- mean(fmado[idx])
    }

    fmado <- fmadoBinned
    dist <- (bins[-1] + bins[-(n.bins+2)])/2
  }

  ext.coeff <- (1 + 2 * fmado) / (1 - 2 * fmado)

  if (length(which) == 2){
    op <- par(mfrow=c(1,2))
    on.exit(par(op))
  }

  if (missing(xlab))
    xlab <- "h"

  if (missing(ylab))
      ylab <- c(expression(nu[F](h)), expression(theta(h)))

  if (add){
    if (any(which == "mado")){
        if (!is.null(angles))
            points(dist, fmado, col = col.angles, ...)

        else
            points(dist, fmado, col = col[1], ...)

      if (fit.curves){
        fmado.fit <- function(h)
          (ext.coeff.fit(h) - 1) / (ext.coeff.fit(h) + 1) / 2

        curve(fmado.fit, from = xlim[1], to = xlim[2], add = TRUE, col = col[2], ...)
      }
    }

    if (any(which == "ext")){

        if (!is.null(angles))
            points(dist, ext.coeff, col = col.angles, ...)

        else
            points(dist, ext.coeff, col = col[1], ...)

      if(fit.curves)
        curve(ext.coeff.fit, from = xlim[1], to = xlim[2], add = TRUE, col = col[2], ...)
    }
  }

  else{
    if (any(which == "mado")){
        if (!is.null(angles))
            plot(dist, fmado, xlab = xlab, ylab = ylab[1], col = col.angles, xlim = xlim, ...)

        else
            plot(dist, fmado, xlab = xlab, ylab = ylab[1], col = col[1], xlim = xlim, ...)

      if (fit.curves){
        fmado.fit <- function(h)
          (ext.coeff.fit(h) - 1) / (ext.coeff.fit(h) + 1) / 2

        curve(fmado.fit, from = xlim[1], to = xlim[2], add = TRUE, col = col[2], ...)
      }
    }

    if (any(which == "ext")){

        if (!is.null(angles))
            plot(dist, ext.coeff, xlab = xlab, ylab = ylab[2], col = col.angles,
                 xlim = xlim, ylim = c(1, max(2, ext.coeff, na.rm = TRUE)), ...)

        else
            plot(dist, ext.coeff, xlab = xlab, ylab = ylab[2], col = col[1],
                 xlim = xlim, ylim = c(1, max(2, ext.coeff, na.rm = TRUE)), ...)

      if (fit.curves)
        curve(ext.coeff.fit, from = xlim[1], to = xlim[2], add = TRUE, col = col[2], ...)
    }
  }

  invisible(cbind(dist = dist, madogram = fmado, ext.coeff = ext.coeff))
}

lmadogram <- function(data, coord, n.bins, xlab, ylab, zlab, n.lambda = 11,
                      marge = "emp", col = terrain.colors(50, alpha = 0.5),
                      theta = 90, phi = 20, border = NA, ...){

  if (is.null(dim(coord))){
    if (length(coord) != ncol(data))
      stop("'data' and 'coord' don't match")
  }

  else if (nrow(coord) != ncol(data))
    stop("'data' and 'coord' don't match")

  if (!(marge %in% c("mle", "emp")))
    stop("'marge' must be either 'mle' or 'emp'")

  n.obs <- nrow(data)
  n.site <- ncol(data)
  n.pairs <- (n.site - 1) * n.site / 2
  lambda <- seq(0, 1, length = n.lambda)
  dist <- distance(coord)

  if (marge == "emp")
    data <- apply(data, 2, rank) / (n.obs + 1)

  else{
    for (i in 1:n.site){
      param <- gevmle(data[,i])
      data[,i] <- pgev(data[,i], param["loc"], param["scale"], param["shape"])
    }
  }

  lmado <- .C("lmadogram", as.double(data), as.integer(n.obs),
              as.integer(n.site), as.double(lambda), as.integer(n.lambda),
              lmado = double(n.pairs * n.lambda),
              PACKAGE = "SpatialExtremes")$lmado
  lmado <- matrix(lmado, nrow = n.lambda, ncol = n.pairs)

  if (anyDuplicated(dist)){
    ##If we have identical pairwise distances, we average over them
    unique.dist <- unique(dist)
    n.pairs <- length(unique.dist)
    lmadoBinned <- matrix(NA, ncol = n.pairs, nrow = n.lambda)

    for (i in 1:n.pairs){
      unique.idx <- which(dist == unique.dist[i])
      lmadoBinned[,i] <- rowMeans(lmado[,unique.idx, drop = FALSE])
    }

    dist <- unique.dist
    lmado <- lmadoBinned
  }

  dist <- sort(dist, index.return = TRUE)
  idx <- dist$ix
  dist <- dist$x
  lmado <- lmado[,idx]

  if (!missing(n.bins)){
    bins <- unique(unname(c(0, quantile(dist, 1:n.bins/(n.bins + 1)), max(dist))))
    n.bins <- length(bins) - 2
    lmadoBinned <- matrix(0, ncol = n.bins + 1, nrow = n.lambda)

    for (k in 1:(n.bins+1)){
      idx <- which((dist <= bins[k+1]) & (dist > bins[k]))

      if (length(idx) > 0)
        lmadoBinned[,k] <- rowMeans(lmado[,idx, drop = FALSE])
    }

    lmado <- lmadoBinned
    dist <- (bins[-1] + bins[-(n.bins+2)])/2
    n.pairs <-  n.bins+1
  }

  if (missing(xlab))
    xlab <- "lambda"

  if (missing(ylab))
    ylab <- "h"

  if (missing(zlab))
    zlab <- "lambda-madogram"

  zfacet <- lmado[-1,-1] + lmado[-1, -n.pairs] + lmado[-n.lambda,-1] +
    lmado[-n.lambda, -n.pairs]

  facetcol <- cut(zfacet, length(col))

  persp(lambda, dist, lmado, xlab = xlab, ylab = ylab, zlab = zlab,
        col = col[facetcol], theta = theta, phi = phi, border = border,
        ...)

  invisible(list(lambda = lambda, dist = dist, madogram = lmado))
}

variogram <- function(data, coord, n.bins, xlab, ylab, angles = NULL,
                      add = FALSE, xlim = c(0, max(dist)), ...){

  if (is.null(dim(coord))){
    if (length(coord) != ncol(data))
      stop("'data' and 'coord' don't match")
  }

  else if (nrow(coord) != ncol(data))
    stop("'data' and 'coord' don't match")

  if (!is.null(angles) & !missing(n.bins))
    stop("'It is not possible to pass 'n.bins' and 'angles' at the same time")

  n.site <- ncol(data)
  n.obs <- nrow(data)
  n.pairs <- (n.site - 1) * n.site / 2
  dist <- distance(coord)

  if (!is.null(angles)){
    distVec <- distance(coord, vec = TRUE)
    n.angles <- length(angles)
    angles.coord <- atan2(distVec[,2], distVec[,1])

    col <- rep(NA, n.site * (n.site - 1) / 2)
    idx.angles <- list()
    for (i in 2:n.angles){
      idx <- which((angles.coord < angles[i]) & (angles.coord >= angles[i-1]))
      idx.angles <- c(idx.angles, list(idx))
      col[idx] <- i-1
    }
  }

  vario <- .C("variogram", as.double(data), as.integer(n.obs),
             as.integer(n.site), vario = double(n.pairs),
             PACKAGE = "SpatialExtremes")$vario

  if (!missing(n.bins)){
    bins <- c(0, quantile(dist, 1:n.bins/(n.bins + 1)), max(dist))
    varioBinned <- rep(NA, length = n.bins + 1)

    for (k in 1:(n.bins + 1)){
      idx <- which((dist <= bins[k+1]) & (dist > bins[k]))

      if (length(idx)>0)
        varioBinned[k] <- mean(vario[idx])
    }

    vario <- varioBinned
    dist <- (bins[-1] + bins[-(n.bins+2)])/2
  }

  if (missing(xlab))
    xlab <- "h"

  if (missing(ylab))
    ylab <- expression(gamma(h))

  if (add)
    points(dist, vario, ...)

  else
    plot(dist, vario, xlab = xlab, ylab = ylab, xlim = xlim, ...)

  invisible(cbind(dist = dist, variogram = vario))
}
