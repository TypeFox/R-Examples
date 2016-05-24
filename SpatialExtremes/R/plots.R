extcoeff <- function(fitted, cov.mod, param, n = 200, xlab, ylab, ...){

  if (missing(fitted) && (missing(cov.mod) || missing(param)))
    stop("You must specify either 'fitted' either 'cov.mod' AND 'param'")

  if (!missing(fitted)){
    if (all(class(fitted) != "maxstab"))
      stop("The 'extcoeff' function is only available for object of class 'maxstab'")
  
    if (ncol(fitted$coord) > 2)
      stop("It's not possible to use this function when the coordinate space has a dimension > 2")
  
    model <- fitted$model
    extCoeff <- fitted$ext.coeff
    param <- fitted$param

    if (fitted$cov.mod != "caugen")
      param <- c(param, smooth2 = 1)
  }

  else{
    if (length(param) != 3)
      stop("You must specify all parameters in the Smith's or Schalther's model")
  
    if (cov.mod == "gauss"){
      model <- "Smith"
      names(param) <- c("cov11", "cov12", "cov22")
      Sigma <- matrix(c(param[1:2], param[2:3]), 2, 2)
      iSigma <- try(solve(Sigma), silent = TRUE)
      
      if (!is.matrix(iSigma) || (det(Sigma) <= 0) ||
          (param[1] <= 0))
        stop("You defined a non semi-definite positive matrix")
      
      extCoeff <- function(h)
        2 * pnorm(0.5 * sqrt(h %*% iSigma %*% h))
    }
    
    else{
      model <- "Schlather"
      names(param) <- c("nugget", "range", "smooth")
      extCoeff <- function(h)
        1 + sqrt(0.5 - 0.5 * (covariance(nugget = param[1], sill = 1 - param[1], range = param[2],
                                         smooth = param[3], smooth2 = param[3],
                                         cov.mod = cov.mod, plot = FALSE, dist = h)))
    }
  }
  
  ##Define an appropriate range for the x-y axis
  if (model == "Smith"){
    A <- matrix(c(param["cov11"], param["cov12"], param["cov12"],
                  param["cov22"]), 2, 2)
    
    eigen.values <- eigen(solve(A))
    eigen.vectors <- eigen.values$vectors
    eigen.values <- eigen.values$values

    r <- 25 / sqrt(2)
    axis1 <- eigen.vectors %*% c(r / sqrt(eigen.values[1]), 0)
    axis2 <- eigen.vectors %*% c(0, r / sqrt(eigen.values[2]))

    x.max <- max(abs(axis1[1]), abs(axis2[1]))
    y.max <- max(abs(axis1[2]), abs(axis2[2]))
    
    x.range <- 1.3 * c(-x.max, x.max)
    y.range <- 1.3 * c(-y.max, y.max)
    
  }

  if (model == "Schlather"){
    fun <- function(h) {
      theta <- extCoeff(h)

      if (theta < 1 + sqrt(0.5))
        abs(1.7 - extCoeff(h))

      else
        h^2
    }
    
    opt1 <- optimize(fun, c(0, 100 * param[2]))$minimum
    y.range <- x.range <- c(-abs(opt1), abs(opt1))
  }

  extcoeff.hat <- matrix(NA, nrow = n, ncol = n)

  xs <- seq(x.range[1], x.range[2], length = n)
  ys <- seq(y.range[1], y.range[2], length = n)

  if (model == "Smith"){
    if (!missing(fitted) && fitted$iso){
      for (i in 1:n)
        for (j in 1:n){
          h <- sqrt(xs[i]^2 + ys[j]^2)
          extcoeff.hat[i,j] <- extCoeff(h)
        }
    }

    else
      for (i in 1:n)
        for (j in 1:n)
          extcoeff.hat[i,j] <- extCoeff(c(xs[i], ys[j]))
  }

  if ((model == "Schlather")){
    for (i in 1:n)
      for (j in 1:n){
        h <- sqrt(xs[i]^2 + ys[j]^2)
        extcoeff.hat[i,j] <- extCoeff(h)
      }
  }

  if (missing(fitted))
    coord.names <- c("x", "y")

  else
    coord.names <- colnames(fitted$coord)

  if (missing(xlab))
    xlab <- coord.names[1]

  if (missing(ylab))
    ylab <- coord.names[2]
  
  contour(xs, ys, extcoeff.hat, xlab = xlab, ylab = ylab, ...)
}
    
map <- function(fitted, x, y, covariates = NULL, param = "quant",
                ret.per = 100, col = terrain.colors(64),
                plot.contour = TRUE, ...){

  if (!(param %in% c("loc", "scale", "shape", "quant")))
    stop("'param' should be one of 'loc', 'scale', 'shape' or 'quant'")

  if (ncol(fitted$coord) > 2)
    stop("It's not possible to use this function when the coordinate space has a dimension > 2")

  if (is.null(covariates) && !is.null(fitted$marg.cov))
    stop("Your model seems to make use of covariate but you supplied none")
  
  if (!is.null(covariates)){
    if (missing(x) || missing(y))
      stop("if 'covariates' is supplied 'x' and 'y' must be too")

    if (!is.array(covariates))
      stop("'covariates' must be an array - see the example")

    covariates.names <- dimnames(covariates)[[3]]
    model.names <- colnames(fitted$marg.cov)

    if (is.null(model.names))
      stop("Your fitted model doesn't seem to make use of covariates")

    if (!all(model.names %in% covariates.names))
      stop("Some required covariates are missing. Please check")

    dim.cov <- dim(covariates)
    
    if ((dim.cov[1] != length(x)) || (dim.cov[2] != length(y)))
      stop("'covariates' doesn't match with 'x' and 'y'")
  }

  else
    covariates.names <- NULL

  if (missing(x)){
    x.range <- range(fitted$coord[,1])
    x <- seq(x.range[1], x.range[2], length = 100)
  }

  if (missing(y)){
    y.range <- range(fitted$coord[,2])
    y <- seq(y.range[1], y.range[2], length = 100)
  }

  n.x <- length(x)
  n.y <- length(y)

  ans <- matrix(NA, nrow = n.x, ncol = n.y)

  if (param == "quant")
    fun <- function(x)
      .qgev(1 - 1/ret.per, x[,"loc"], x[,"scale"], x[,"shape"])

  else
    fun <- function(x)
      x[,param]

  new.data <- matrix(NA, nrow = n.x * n.y, ncol = 2 + length(covariates.names))

  for (i in 1:n.x)
    new.data[(n.y * (i - 1) + 1):(n.y * i),] <- cbind(x[i], y, covariates[i,,])

  colnames(new.data) <- c(colnames(fitted$coord), covariates.names)

  param.hat <- predict(fitted, new.data, std.err = FALSE)

  ans <- matrix(fun(param.hat), n.x, n.y, byrow = TRUE)

  image(x, y, ans, ..., col = col)

  if (plot.contour)
    contour(x, y, ans, add = TRUE)

  invisible(list(x = x, y = y, z = ans))
}

condmap <- function(fitted, fix.coord, x, y, covariates = NULL,
                    ret.per1 = 100, ret.per2 = ret.per1,
                    col = terrain.colors(64), plot.contour = TRUE, ...){

  if (ncol(fitted$coord) > 2)
    stop("It's not possible to use this function when the coordinate space has a dimension >= 2")

  if (is.null(covariates) && !is.null(fitted$marg.cov))
    stop("Your model seems to make use of covariate but you supplied none")
  
  if (!is.null(covariates)){
    if (missing(x) || missing(y))
      stop("if 'covariates' is supplied 'x' and 'y' must be too")

    if (!is.array(covariates))
      stop("'covariates' must be an array - see the example")

    covariates.names <- dimnames(covariates)[[3]]
    model.names <- colnames(fitted$marg.cov)

    if (is.null(model.names))
      stop("Your fitted model doesn't seem to make use of covariates")

    if (!all(model.names %in% covariates.names))
      stop("Some required covariates are missing. Please check")

    dim.cov <- dim(covariates)
    
    if ((dim.cov[1] != length(x)) || (dim.cov[2] != length(y)))
      stop("'covariates' doesn't match with 'x' and 'y'")
  }

  else
    covariates.names <- NULL

  if (missing(x)){
    x.range <- range(fitted$coord[,1])
    x <- seq(x.range[1], x.range[2], length = 100)
  }

  if (missing(y)){
    y.range <- range(fitted$coord[,2])
    y <- seq(y.range[1], y.range[2], length = 100)
  }

  if (any(x == fix.coord[1]) & any(y == fix.coord[2]))
    stop("'fix.coord' belongs to your grid. This is not possible!")
  
  n.x <- length(x)
  n.y <- length(y)

  ans <- double(n.x * n.y)
  new.data <- matrix(NA, nrow = n.x * n.y, ncol = 2 + length(covariates.names))

  for (i in 1:n.x)
    new.data[(n.y * (i - 1) + 1):(n.y * i),] <- cbind(x[i], y, covariates[i,,])

  colnames(new.data) <- c(colnames(fitted$coord), covariates.names)
  param <- predict(fitted, new.data, std.err = FALSE)

  ##z1: quantile related to ret.per1
  z1 <- - 1 / log(1 - 1/ret.per1)

  if (fitted$model == "Smith"){
    Sigma <- matrix(c(fitted$param[1:2], fitted$param[2:3]), 2)
    iSigma <- solve(Sigma)

    cond.prob <- function(z2){
      c1 <- 0.5 * a + log(z2/z1) / a
      c2 <- a - c1
      - 1/ret.per2 + ret.per1 * (1 / ret.per1 - exp(-1/z2) +
                                   exp(-pnorm(c1) / z1 - pnorm(c2) / z2))
    }

    for (i in 1:(n.x * n.y)){
      delta <- fix.coord - new.data[i,1:2]
      a <- sqrt(delta %*% iSigma %*% delta)
      ans[i] <- uniroot(cond.prob, c(1e-4, 1e5))$root
    }
  }    

  if (fitted$model == "Schlather"){
    cond.prob <- function(z2)
      ret.per1 * (exp(-0.5 * (1/z1 + 1/z2) *
                      (1 + sqrt(1 - 2 * (fitted$cov.fun(h) + 1) *
                                z1 * z2 / (z1 + z2)^2))) +
                  1 / ret.per1 - exp(-1/z2)) - 1 / ret.per2

    for (i in 1:(n.x * n.y)){
      h <- sqrt(sum((fix.coord - new.data[i,1:2])^2))
      ans[i] <- uniroot(cond.prob, c(1e-4, 1e5))$root
    }
  }    

  if (fitted$model == "Geometric"){
    cond.prob <- function(z2){
      c1 <-  0.5 * a + log(z2/z1)  /a
      c2 <- a - c1
      
      - 1/ret.per2 + ret.per1 * (1 / ret.per1 - exp(-1/z2) +
                                 exp(-pnorm(c1) / z1 - pnorm(c2) / z2))
    }

    for (i in 1:(n.x * n.y)){
      h <- sqrt(sum((fix.coord - new.data[i,1:2])^2))
      a <- sqrt(2 * fitted$par["sigma2"] * (1 - fitted$cov.fun(h)))
      ans[i] <- uniroot(cond.prob, c(1e-4, 1e5))$root
    }
  }

  if (fitted$model == "Brown-Resnick"){
    cond.prob <- function(z2){
      c1 <-  0.5 * a + log(z2/z1)  /a
      c2 <- a - c1
      
      - 1/ret.per2 + ret.per1 * (1 / ret.per1 - exp(-1/z2) +
                                 exp(-pnorm(c1) / z1 - pnorm(c2) / z2))
    }

    for (i in 1:(n.x * n.y)){
      h <- sqrt(sum((fix.coord - new.data[i,1:2])^2))
      a <- (h / fitted$param["range"])^(0.5 * fitted$param["smooth"])
      ans[i] <- uniroot(cond.prob, c(1e-4, 1e5))$root
    }
  }
  
  ans <- matrix(.frech2gev(ans, param[,"loc"], param[,"scale"], param[,"shape"]),
                n.x, n.y, byrow = TRUE)
  
  image(x, y, ans, ..., col = col)

  if (plot.contour)
    contour(x, y, ans, add = TRUE)

  invisible(list(x = x, y = y, z = ans))
}

map.latent <- function(fitted, x, y, covariates = NULL, param = "quant",
                       ret.per = 100, col = terrain.colors(64),
                       plot.contour = TRUE, fun = mean, level = 0.95,
                       show.data = TRUE, control = list(nlines = 500), ...){

  if (!(param %in% c("loc", "scale", "shape", "quant")))
    stop("'param' should be one of 'loc', 'scale', 'shape' or 'quant'")

  if (ncol(fitted$coord) > 2)
    stop("It's not possible to use this function when the coordinate space has a dimension > 2")

  if (is.null(covariates) && !is.null(fitted$marg.cov))
    stop("Your model seems to make use of covariate but you supplied none")
  
  if (!is.null(covariates)){
    if (missing(x) || missing(y))
      stop("if 'covariates' is supplied 'x' and 'y' must be too")

    if (!is.array(covariates))
      stop("'covariates' must be an array - see the example")

    covariates.names <- dimnames(covariates)[[3]]
    model.names <- colnames(fitted$marg.cov)

    if (is.null(model.names))
      stop("Your fitted model doesn't seem to make use of covariates")

    if (!all(model.names %in% covariates.names))
      stop("Some required covariates are missing. Please check")

    dim.cov <- dim(covariates)
    
    if ((dim.cov[1] != length(x)) || (dim.cov[2] != length(y)))
      stop("'covariates' doesn't match with 'x' and 'y'")
  }

  else
    covariates.names <- NULL

  if (missing(x)){
    x.range <- range(fitted$coord[,1])
    x <- seq(x.range[1], x.range[2], length = 100)
  }

  if (missing(y)){
    y.range <- range(fitted$coord[,2])
    y <- seq(y.range[1], y.range[2], length = 100)
  }

  n.x <- length(x)
  n.y <- length(y)

  n.chain <- nrow(fitted$chain.loc)
  n.loccoeff <- ncol(fitted$loc.dsgn.mat)
  n.scalecoeff <- ncol(fitted$scale.dsgn.mat)
  n.shapecoeff <- ncol(fitted$shape.dsgn.mat)

  ans <- array(NA, c(n.x, n.y, n.chain))

  new.data <- matrix(NA, nrow = n.x * n.y, ncol = 2 + length(covariates.names))

  for (i in 1:n.x)
    new.data[(n.y * (i - 1) + 1):(n.y * i),] <- cbind(x[i], y, covariates[i,,])

  colnames(new.data) <- c(colnames(fitted$coord), covariates.names)

  loc.dsgn.mat <- modeldef(new.data, fitted$loc.form)$dsgn.mat
  scale.dsgn.mat <- modeldef(new.data, fitted$scale.form)$dsgn.mat
  shape.dsgn.mat <- modeldef(new.data, fitted$shape.form)$dsgn.mat

  if (param == "loc"){
    for (i in 1:n.chain){
      loc <- matrix(loc.dsgn.mat %*% fitted$chain.loc[i, 1:n.loccoeff], n.x, n.y,
                    byrow = TRUE)
      res <- as.numeric(fitted$chain.loc[i,-(1:(n.loccoeff+3))] -
                        fitted$loc.dsgn.mat %*% fitted$chain.loc[i,1:n.loccoeff])
      ans[,,i] <- loc + condrgp(1, cbind(x, y), fitted$coord, res,
                                cov.mod = fitted$cov.mod[1],
                                sill = fitted$chain.loc[i,"sill"], range = fitted$chain.loc[i,"range"],
                                smooth = fitted$chain.loc[i,"smooth"], grid = TRUE,
                                control = control)$cond.sim
    }
  }

  else if (param == "scale"){
    for (i in 1:n.chain){
      scale <- matrix(scale.dsgn.mat %*% fitted$chain.scale[i, 1:n.scalecoeff], n.x, n.y,
                        byrow = TRUE)
      res <- as.numeric(fitted$chain.scale[i,-(1:(n.scalecoeff+3))] -
                        fitted$scale.dsgn.mat %*% fitted$chain.scale[i,1:n.scalecoeff])

      n.trials <- 1
      flag <- FALSE
      while (!flag){
        scale <- scale + condrgp(1, cbind(x, y), fitted$coord, res,
                                 cov.mod = fitted$cov.mod[2],
                                 sill = fitted$chain.scale[i,"sill"],
                                 range = fitted$chain.scale[i,"range"],
                                 smooth = fitted$chain.shape[i,"smooth"], grid = TRUE,
                                 control = control)$cond.sim
        
        flag <- all(scale > 0)
        n.trials <- n.trials + 1

        if (n.trials >= 100){
          ans[,,i] <- NA
          break
        }
      }

      ans[,,i] <- scale

    }
  }

  else if (param == "shape"){
    for (i in 1:n.chain){
      shape <- matrix(shape.dsgn.mat %*% fitted$chain.shape[i, 1:n.shapecoeff], n.x, n.y,
                      byrow = TRUE)
      res <- as.numeric(fitted$chain.shape[i,-(1:(n.shapecoeff+3))] -
                        fitted$shape.dsgn.mat %*% fitted$chain.shape[i,1:n.shapecoeff])
    
      ans[,,i] <- shape + condrgp(1, cbind(x, y), fitted$coord, res,
                                  cov.mod = fitted$cov.mod[1],
                                  sill = fitted$chain.shape[i,"sill"], range = fitted$chain.shape[i,"range"],
                                  smooth = fitted$chain.shape[i,"smooth"], grid = TRUE,
                                  control = control)$cond.sim
    }
  }

  else {
    for (i in 1:n.chain){
      loc <- matrix(loc.dsgn.mat %*% fitted$chain.loc[i, 1:n.loccoeff], n.x, n.y,
                    byrow = TRUE)
      res <- as.numeric(fitted$chain.loc[i,-(1:(n.loccoeff+3))] -
                        fitted$loc.dsgn.mat %*% fitted$chain.loc[i,1:n.loccoeff])
      
      loc <- loc + condrgp(1, cbind(x, y), fitted$coord, res,
                           cov.mod = fitted$cov.mod[1],
                           sill = fitted$chain.loc[i,"sill"], range = fitted$chain.loc[i,"range"],
                           smooth = fitted$chain.loc[i,"smooth"], grid = TRUE,
                           control = control)$cond.sim
      
      flag <- FALSE
      scale <- matrix(scale.dsgn.mat %*% fitted$chain.scale[i, 1:n.scalecoeff], n.x, n.y,
                        byrow = TRUE)
      res <- as.numeric(fitted$chain.scale[i,-(1:(n.scalecoeff+3))] -
                        fitted$scale.dsgn.mat %*% fitted$chain.scale[i,1:n.scalecoeff])

      n.trials <- 1
      while (!flag){        
        scale <- scale + condrgp(1, cbind(x, y), fitted$coord, res,
                                 cov.mod = fitted$cov.mod[2],
                                 sill = fitted$chain.scale[i,"sill"],
                                 range = fitted$chain.scale[i,"range"],
                                 smooth = fitted$chain.shape[i,"smooth"], grid = TRUE,
                                 control = control)$cond.sim
        
        flag <- all(scale > 0)
        n.trials <- n.trials + 1

        if (n.trials >= 100)
          break
      }

      if (n.trials < 100){
        shape <- matrix(shape.dsgn.mat %*% fitted$chain.shape[i, 1:n.shapecoeff], n.x, n.y,
                        byrow = TRUE)
        res <- as.numeric(fitted$chain.shape[i,-(1:(n.shapecoeff+3))] -
                          fitted$shape.dsgn.mat %*% fitted$chain.shape[i,1:n.shapecoeff])
        
        shape <- shape + condrgp(1, cbind(x, y), fitted$coord, res,
                                 cov.mod = fitted$cov.mod[1],
                                 sill = fitted$chain.shape[i,"sill"], range = fitted$chain.shape[i,"range"],
                                 smooth = fitted$chain.shape[i,"smooth"], grid = TRUE,
                                 control = control)$cond.sim
        
        ans[,,i] <- .qgev(1 - 1 / ret.per, loc, scale, shape)
      }

      else
        ans[,,i] <- NA
    }
  }
  
  prob.low <- (1 - level) /2
  prob.up <- 1 - prob.low
  post.sum <- ci.low <- ci.up <- matrix(NA, n.x, n.y)

  for (i in 1:n.x){
    post.sum[i,] <- apply(ans[i,,], 1, fun, na.rm = TRUE)
    ci.low[i,] <- apply(ans[i,,], 1, quantile, prob.low, na.rm = TRUE)
    ci.up[i,] <- apply(ans[i,,], 1, quantile, prob.up, na.rm = TRUE)
  }

  op <- par(no.readonly = TRUE)
  layout(matrix(1:6,1), widths = rep(c(1, .4), 3))
  
  for (i in 1:3){
    par(mar = rep(0.5, 4))

    dummy <- switch(as.character(i), "1" = ci.low, "2" = post.sum, "3" = ci.up)
    breaks <- seq(min(dummy), max(dummy), length = length(col) + 1)
    image(x, y, dummy, ..., col = col, breaks = breaks, xaxt = "n",
          yaxt = "n")
  
    if (plot.contour)
      contour(x, y, dummy, add = TRUE)

    if (show.data)
      points(fitted$coord)

    mar <- c(0.5, 1, 0.5, 3)
    par(las = 1, pty = "m", mar = mar)

    plot.new()
    plot.window(xlim = c(0,1), ylim = range(breaks), xaxs = "i", yaxs = "i")
    rect(0, breaks[-length(breaks)], 1, breaks[-1], col = col, border = NA)
    axis(4, at = pretty(dummy))
    box()
  }

  on.exit(par(op))
  invisible(list(x = x, y = y, post.sum = post.sum, ci.low = ci.low, ci.up = ci.up))
}
