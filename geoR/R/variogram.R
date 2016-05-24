"variog" <-
  function(geodata, coords = geodata$coords, data = geodata$data, 
       uvec = "default", breaks = "default",
       trend = "cte", lambda = 1,
       option = c("bin", "cloud", "smooth"),
       estimator.type = c("classical", "modulus"), 
       nugget.tolerance, max.dist, pairs.min = 2,
       bin.cloud = FALSE, direction = "omnidirectional", tolerance = pi/8,
       unit.angle = c("radians","degrees"), angles = FALSE, 
       messages, ...) 
{
  if(missing(geodata)) geodata <- list(coords=coords, data = data)
  call.fc <- match.call()
  if(missing(messages))
    messages.screen <- as.logical(ifelse(is.null(getOption("geoR.messages")), TRUE, getOption("geoR.messages")))
  else messages.screen <- messages
  keep <- list(...)
  if(is.null(keep$keep.NA)) keep.NA <- FALSE
  else keep.NA <- keep$keep.NA
  ##
  ## Directional variogram - setting angles
  ##
  unit.angle <- match.arg(unit.angle)
  if(mode(direction) == "numeric"){
    if(length(direction) > 1) stop("only one direction is allowed")
    if(length(tolerance) > 1) stop("only one tolerance value is allowed")
    if(unit.angle == "degrees"){
      ang.deg <- direction
      ang.rad <- (ang.deg * pi)/180
      tol.deg <- tolerance
      tol.rad <- (tol.deg * pi)/180
    }
    else{
      ang.rad <- direction
      ang.deg <- (ang.rad * 180)/pi
      tol.rad <- tolerance
      tol.deg <- (tol.rad * 180)/pi
    }
    if(ang.rad > pi | ang.rad < 0)
      stop("direction must be an angle in the interval [0,pi[ radians")
    if(tol.rad > pi/2 | tol.rad < 0)
      stop("tolerance must be an angle in the interval [0,pi/2] radians")
    if(tol.deg >= 90){
      direction <- "omnidirectional"
      cat("variog: computing omnidirectional variogram\n")
    }
    else{
      if(messages.screen){
        cat(paste("variog: computing variogram for direction = ", round(ang.deg, digits=3), " degrees (", round(ang.rad, digits=3), " radians)\n", sep=""))
        cat(paste("        tolerance angle = ", round(tol.deg, digits=3), " degrees (", round(tol.rad, digits=3), " radians)\n", sep=""))
      }
    }
  }
  else if(messages.screen) cat("variog: computing omnidirectional variogram\n")
  ##   
  ##
  coords <- as.matrix(coords)
  data <- as.matrix(data)
  if(nrow(coords) != nrow(data)) stop("coords and data have incompatible dimensions") 
  data.var <- apply(data, 2, var)
  n.data <- nrow(coords)
  n.datasets <- ncol(data)
  data <- drop(data)
  ##
  ## variogram estimator
  ##
  option <- match.arg(option)
  estimator.type <- match.arg(estimator.type)
  ##
  ## transformation
  ##
  if (abs(lambda - 1) > 0.0001) {
    if (abs(lambda) < 0.0001) data <- log(data)
    else data <- ((data^lambda) - 1)/lambda
  }
  ##
  ## trend removal
  ##
  xmat <- unclass(trend.spatial(trend = trend, geodata = geodata))
  if (nrow(xmat) != n.data) 
    stop("coords and trend have incompatible sizes")
  if (trend != "cte") {
    if (is.vector(data)) {
      temp.fit <- lm(data ~ xmat + 0)
      beta.ols <- temp.fit$coeff
      data <- temp.fit$residuals
      temp.fit <- NULL
      names(data) <- NULL
    }
    else {
      only.res <- function(y, x)
        lm(y ~ xmat + 0)$residuals
      data <- apply(data, 2, only.res, x = xmat)
      only.beta <- function(y, x)
        lm(y ~ xmat + 0)$coef
      beta.ols <- apply(data, 2, only.beta, x = xmat)
    }
  }
  else beta.ols <- colMeans(as.matrix(data))
  ##
  ## Defining bins
  ##
  u <- as.vector(dist(as.matrix(coords)))
  if(missing(nugget.tolerance) || nugget.tolerance < 1e-11){
    nugget.tolerance <- 1e-12
    nt.ind <- FALSE
  }
  else{
    if(mode(nugget.tolerance) != "numeric") stop("nugget.tolerance must be numeric")
    nt.ind <- TRUE
  }
  min.dist <- min(u)
  if(min.dist < nugget.tolerance) nt.ind <- TRUE
  ##
  ## directional
  ##
  if(direction != "omnidirectional" | angles){
    u.ang <- .C("tgangle",
                as.double(as.vector(coords[,1])),
                as.double(as.vector(coords[,2])),
                as.integer(dim(coords)[1]),
                res = as.double(rep(0, length(u))),
                PACKAGE = "geoR")$res
    if(any(is.na(u.ang)))
      stop("NA returned in angle calculations maybe due to co-located data")
    u.ang <- atan(u.ang)
    u.ang[u.ang < 0] <- u.ang[u.ang < 0] + pi
  }
  ##
  if (option == "bin" && bin.cloud == FALSE && direction == "omnidirectional") {
    if (missing(max.dist)) umax <- max(u)  
    else umax <- max(u[u < max.dist])
    dbins <- .define.bins(max.dist = umax, uvec = uvec, breaks = breaks, nugget.tolerance = nugget.tolerance)
    uvec <- dbins$uvec ; bins.lim <- dbins$bins.lim
    nbins <- length(bins.lim) - 1
    if (missing(max.dist)) max.dist <- max(bins.lim)
    if(bins.lim[1] < 1e-16) bins.lim[1] <- -1
    bin.f <- function(data) {
      cbin <- vbin <- sdbin <- rep(0, nbins)
      .C("binit", as.integer(n.data),
         as.double(as.vector(coords[, 1])),
         as.double(as.vector(coords[, 2])), as.double(as.vector(data)), 
         as.integer(nbins), as.double(as.vector(bins.lim)), 
         as.integer(estimator.type == "modulus"), as.double(max.dist), 
         cbin = as.integer(cbin), vbin = as.double(vbin), 
         as.integer(TRUE), sdbin = as.double(sdbin),
         PACKAGE = "geoR")[c("vbin", "cbin", "sdbin")]
    }
    result <- array(unlist(lapply(as.data.frame(data), bin.f)), 
                    dim = c(nbins, 3, n.datasets))
    indp <- (result[, 2, 1] >= pairs.min)
    result[!indp,1,] <- NA
    if(bins.lim[1] < 0) bins.lim[1] <- 0
    if(!nt.ind){
      uvec <- uvec[-1]
      indp <- indp[-1]
      bins.lim <- bins.lim[-1]
      result <- result[-1,,,drop=FALSE]
    }
    if(keep.NA)
      result <- list(u = uvec, v = result[, 1, ], 
                     n = result[, 2, 1], sd = result[, 3, ], 
                     bins.lim = bins.lim, 
                     ind.bin = indp)
    else
      result <- list(u = uvec[indp], v = result[indp, 1, ], 
                     n = result[indp, 2, 1], sd = result[indp, 3, ],
                     bins.lim = bins.lim, ind.bin = indp)
  }
  else {
    data <- as.matrix(data)
    v <- matrix(0, nrow = length(u), ncol = n.datasets)
    for (i in 1:n.datasets) {
      v[, i] <- as.vector(dist(data[, i]))
      if (estimator.type == "modulus") 
        v[, i] <- v[, i,drop=FALSE]^(0.5)
      else v[, i] <- (v[, i,drop=FALSE]^2)/2
    }
    if (!missing(max.dist)) {
      v <- v[u <= max.dist,,drop=FALSE]
      if(direction != "omnidirectional")
        u.ang <- u.ang[u <= max.dist]
      u <- u[u <= max.dist]
    }
    if(direction != "omnidirectional"){
      ang.lower <- ang.rad - tol.rad
      ang.upper <- ang.rad + tol.rad
     if(ang.lower >= 0 & ang.upper < pi)
        ang.ind <- (!is.na(u.ang) & ((u.ang >= ang.lower) & (u.ang <= ang.upper)))
      if(ang.lower < 0)
        ang.ind <- (!is.na(u.ang) & ((u.ang < ang.upper) | (u.ang > (pi + ang.lower))))
      if(ang.upper >= pi)
        ang.ind <- (!is.na(u.ang) & ((u.ang > ang.lower) | (u.ang < (ang.upper - pi))))
      ##ang.ind <- ((u.ang >= ang.rad - tol.rad)&(u.ang <= ang.rad + tol.rad))
      v <- v[ang.ind,,drop=FALSE]
      u <- u[ang.ind]
    }
    data <- drop(data)
    v <- drop(v)
    if (option == "cloud"){
      result <- list(u = u, v = v)
      if(angles) result$angles <- u.ang
    }
    if (option == "bin") {
      if (missing(max.dist)) umax <- max(u)
      else umax <- max(u[u < max.dist])
      if(bin.cloud == 'diff') dd <- diffpairs(coords,data)$diff
      else dd <- 0
      result <- .rfm.bin(cloud = list(u = u, v = v, d = dd),
                        estimator.type = estimator.type, 
                        uvec = uvec, breaks = breaks, nugget.tolerance = nugget.tolerance, 
                        bin.cloud = bin.cloud, max.dist = umax, keep.NA = keep.NA)
      if(keep.NA){
        if (pairs.min > 0) {
          indp <- (result$n < pairs.min)
          if(!nt.ind){
            for(i in 1:5) result[[i]] <- result[[i]][-1]
            indp <- indp[-1]
          }
          if (is.matrix(result$v)) {
            result$v[indp, ] <- result$sd[indp, ] <- NA
          }
          else {
            result$v[indp] <- result$sd[indp] <- NA
          }
        }
        result$ind.bin <- indp
      }
      else{
        if (pairs.min > 0) {
#          indp <- (result$n >= pairs.min)
          if(!nt.ind){
            for(i in 1:5) result[[i]] <- result[[i]][-1]
#            indp <- indp[-1]
          }
          indp <- (result$n >= pairs.min)
          if (is.matrix(result$v)) {
            result$v <- result$v[indp, ]
            result$sd <- result$sd[indp, ]
          }
          else {
            result$v <- result$v[indp]
            result$sd <- result$sd[indp]
          }
          result$u <- result$u[indp]
          result$n <- result$n[indp]
        }
        result$ind.bin <- indp
      }
    }
    if (option == "smooth") {
      if (is.matrix(v)) stop("smooth not yet available for more than one data-set")
      temp <- ksmooth(u, v, ...)
      result <- list(u = temp[[1]], v = temp[[2]])
    }
    if(missing(max.dist)) max.dist <- max(u)
  }
  if(nt.ind){
    if(!exists(".variog4.nomessage",where=1)) cat("variog: co-locatted data found, adding one bin at the origin\n")
    if(all(result$u[1:2] < 1e-11)) result$u[2] <- sum(result$bins.lim[2:3])/2
  }
  result <- c(result, list(var.mark = data.var, beta.ols = beta.ols,
                           output.type = option, max.dist = max.dist, 
                           estimator.type = estimator.type, n.data = n.data,
                           lambda = lambda, trend = trend, pairs.min = pairs.min))
  result$nugget.tolerance <- nugget.tolerance
  if(direction != "omnidirectional") result$direction <- ang.rad
  else result$direction <- "omnidirectional"
  if(direction != "omnidirectional") result$tolerance <- tol.rad
  else result$tolerance <- "none" 
  result$uvec <- uvec
  result$call <- call.fc
  oldClass(result) <- "variogram"
  return(result)
}

"variog4" <-
  function (geodata, coords = geodata$coords, data = geodata$data, 
            uvec = "default", breaks = "default", trend = "cte",
            lambda = 1, option = c("bin", "cloud", "smooth"),
            estimator.type = c("classical", "modulus"), 
            nugget.tolerance, max.dist, pairs.min = 2,
            bin.cloud = FALSE, direction = c(0, pi/4, pi/2, 3*pi/4),
            tolerance = pi/8, unit.angle = c("radians", "degrees"),
            messages, ...) 
{
  if(missing(geodata)) geodata <- list(coords = coords, data = data)
  if(missing(messages))
    messages.screen <- as.logical(ifelse(is.null(getOption("geoR.messages")), TRUE, getOption("geoR.messages")))
  else messages.screen <- messages
  if(missing(nugget.tolerance)) nugget.tolerance <- 1e-12
  u <- as.vector(dist(as.matrix(coords)))
  if(length(direction) != 4)
    stop("argument direction must be a vector with 4 values. For different specification use the function variog()")
  if(length(tolerance) != 1)
    stop("only one value can be provided to the argument tolerance. For different specification use the function variog()")
  res <- list()
  unit.angle <- match.arg(unit.angle)
  if(unit.angle == "radians") dg <- direction * 180/pi
  else dg <- direction
  if (missing(max.dist)) umax <- max(u)  
  else umax <- max(u[u < max.dist])
  assign(".variog4.nomessage", TRUE, pos=1)
  for(angle in direction){
    res[[as.character(round(dg[which(direction == angle)], digits=1))]] <-
      variog(geodata=geodata,
             uvec=uvec, breaks = breaks, trend = trend,
             lambda = lambda, option = option,
             estimator.type = estimator.type,
             nugget.tolerance = nugget.tolerance,
             max.dist = max.dist,
             pairs.min = pairs.min,
             bin.cloud = bin.cloud,
             direction = angle,
             tolerance = tolerance,
             unit.angle = unit.angle,
             messages = messages.screen, keep.NA = TRUE)
    NULL
  }
  if (exists(".variog4.nomessage", where=1)) remove(".variog4.nomessage", pos=1, inherits = TRUE)
  res$omnidirectional <- variog(geodata=geodata,
                                uvec=uvec, breaks = breaks,
                                trend = trend,
                                lambda = lambda, option = option,
                                estimator.type = estimator.type,
                                nugget.tolerance = nugget.tolerance,
                                max.dist = max.dist,
                                pairs.min = pairs.min,
                                bin.cloud = bin.cloud,
                                direction = "omnidirectional",
                                tolerance = tolerance,
                                unit.angle = unit.angle,
                                messages = messages.screen,
                                keep.NA = TRUE 
                                )
  oldClass(res) <- "variog4"
  return(res) 
}




"plot.variog4" <-
  function (x, omnidirectional = FALSE, same.plot = TRUE,
            legend = TRUE,...)
{
  fc.call <- match.call()
  ymax <- max(c(x[[1]]$v, x[[2]]$v, x[[3]]$v, x[[4]]$v), na.rm=TRUE)
  n.o <- names(x)[1:4]
  GP <- list(...)
  argsnames <- c(names(formals(plot.default)), names(par()))
  names(GP) <- argsnames[charmatch(names(GP), argsnames)]
  if(is.null(GP$xlab)) GP$xlab <- fc.call$xlab <- "distance"
  if(is.null(GP$ylab)) GP$ylab <- fc.call$ylab <- "semivariance"
  if (same.plot) {
    xx <- x[[5]]$u
    yy <- cbind(x[[1]]$v, x[[2]]$v, x[[3]]$v, x[[4]]$v)
    if (omnidirectional) yy <- cbind(x[[5]]$v, yy)
    ## OFC, changed : if (is.null(GP$lty)) GP$lty <- 1:5
    if (is.null(GP$lty)){
      if(omnidirectional) GP$lty <- fc.call$lty <- 1:5
      else GP$lty <- fc.call$lty <- 1:4
    }
    if (is.null(GP$lwd)) GP$lwd <- 1
    ## OFC, changed : if (is.null(GP$col)) GP$col <- 1:5
    if (is.null(GP$col)){
      if(omnidirectional) GP$col <- 1:5
      else GP$col <- 1:4
    }
    if (is.null(GP$pch)) GP$pch <- NULL
    if (is.null(GP$type)) GP$type <- "l"
    matplot(x = xx, y = yy, type = GP$type, lty=GP$lty, lwd=GP$lwd,
            col=GP$col, pch=GP$pch, xlab=GP$xlab, ylab=GP$ylab,
            xlim = c(0, max(xx)), ylim=c(0,max(yy, na.rm=TRUE)))
    #if (is.null(GP$xlim)) fc.call$xlim <- c(0, max(xx))
    #if (is.null(GP$xlim)) fc.call$ylim <- c(0,max(yy, na.rm=TRUE)))
    #fc.call[[1]] <- matplot
    #eval(fc.call, sys.frame(parent.frame()))
    if (legend) {
      if (omnidirectional) {
        legend(0, ymax,
               legend = c("omnid.",
                 substitute(a * degree, list(a = n.o[1])),
                 substitute(a * degree, list(a = n.o[2])),
                 substitute(a * degree, list(a = n.o[3])),
                 substitute(a * degree, list(a = n.o[4])),
                 expression()),
               lty = GP$lty, lwd = GP$lwd, col = GP$col)
      }
      else {
        legend(0, ymax,
               legend = c(substitute(a * degree,
                 list(a = n.o[1])), substitute(a * degree, list(a = n.o[2])),
                 substitute(a * degree, list(a = n.o[3])), substitute(a *
                                               degree, list(a = n.o[4])), expression()),
               lty = GP$lty, lwd = GP$lwd, col = GP$col)
      }
    }
  }
  else {
    temp.mf <- par()$mfrow
    par(mfrow = c(2, 2))
    if (is.null(GP$lty)) {
      GP$lty <- rep(1, 4)
      if (omnidirectional)
        GP$lty <- c(GP$lty, 2)
    }
    else {
      if (length(GP$lty) == 1)
        if (omnidirectional)
          GP$lty <- rep(GP$lty, 5)
        else GP$lty <- rep(GP$lty, 4)
      if (length(GP$lty) == 2)
        if (omnidirectional)
          GP$lty <- c(rep(GP$lty[1], 4), GP$lty[2])
        else GP$lty <- c(rep(GP$lty, 4))
      if (length(GP$lty) == 4 & omnidirectional)
        GP$lty <- c(rep(GP$lty, 2))
    }
    if (is.null(GP$lwd)) {
      GP$lwd <- rep(1, 4)
      if (omnidirectional)
        GP$lwd <- c(GP$lwd, 1)
    }
    else {
      if (length(GP$lwd) == 1)
        if (omnidirectional)
          GP$lwd <- rep(GP$lwd, 5)
        else GP$lwd <- rep(GP$lwd, 4)
      if (length(GP$lwd) == 2)
        if (omnidirectional)
          GP$lwd <- c(rep(GP$lwd[1], 4), GP$lwd[2])
        else GP$lwd <- rep(GP$lwd, 4)
      if (length(GP$lwd) == 4 & omnidirectional)
        GP$lwd <- c(rep(GP$lwd, 1))
    }
    if (is.null(GP$col)) {
      GP$col <- rep(1, 4)
      if (omnidirectional) GP$col <- c(GP$col, 1)
    }
    else {
      if (length(GP$col) == 1)
        if (omnidirectional) GP$col <- rep(GP$col, 5)
        else GP$col <- rep(GP$col, 4)
      if (length(GP$col) == 2)
        if (omnidirectional)
          GP$col <- c(rep(GP$col[1], 4), GP$col[2])
        else GP$col <- rep(GP$col, 2)
      if (length(GP$col) == 4 & omnidirectional)
        GP$col <- c(rep(GP$col, 1))
    }
    if (is.null(GP$pch)) {
      GP$pch <- rep(1, 4)
      if (omnidirectional)
        GP$pch <- c(GP$pch, 1)
    }
    else {
      if (length(GP$pch) == 1)
        if (omnidirectional)
          GP$pch <- rep(GP$pch, 5)
        else GP$pch <- rep(GP$pch, 4)
      if (length(GP$pch) == 2)
        if (omnidirectional)
          GP$pch <- c(rep(GP$pch[1], 4), GP$pch[2])
        else GP$pch <- rep(GP$pch, 2)
      if (length(GP$pch) == 4 & omnidirectional)
        GP$pch <- c(rep(GP$pch, 2))
    }
    if (is.null(GP$type)) {
      GP$type <- rep("l", 4)
      if (omnidirectional)
        GP$type <- c(GP$type, "l")
    }
    else {
      if (length(GP$type) == 1)
        if (omnidirectional)
          GP$type <- rep(GP$type, 5)
        else GP$type <- rep(GP$type, 4)
      if (length(GP$type) == 2 & omnidirectional)
        GP$type <- c(rep(GP$type[1], 4), GP$type[2])
      if (length(GP$type) == 4 & omnidirectional)
        GP$type <- c(rep(GP$type, 2))
    }
    for (i in 1:4) {
      plot.default(x[[i]]$u, x[[i]]$v,
                   xlim = c(0, max(x[[i]]$u)), ylim = c(0, ymax),
                   type = GP$type[i],
           col = GP$col[i], lwd = GP$lwd[i], lty = GP$lty[i],
           pch = GP$pch[i], xlab=GP$xlab, ylab=GP$ylab)
      if (omnidirectional) {
        lines(x$omnidirectional, type = GP$type[5],
              col = GP$col[5], lwd = GP$lwd[5], lty = GP$lty[5])
        legend(0, ymax, legend = c(substitute(a * degree,
                          list(a = n.o[i])), "omn.", expression()),
               lty = c(GP$lty[i], GP$lty[5]),
               col=c(GP$col[i],  GP$col[5]),
               lwd=c(GP$lwd[i],  GP$lwd[5]))
      }
      else title(main = substitute(a * degree, list(a = n.o[i])),
                 cex = 1.3)
    }
    par(mfrow = temp.mf)
  }
  return(invisible())
}

".rfm.bin" <-
  function (cloud, l = 13, uvec = "default", breaks = "default", nugget.tolerance, 
            estimator.type = c("classical", "modulus"), bin.cloud = FALSE,
            max.dist, keep.NA = FALSE)
{
  estimator.type <- match.arg(estimator.type)
  dbins <- .define.bins(max.dist = max(cloud$u), uvec = uvec, breaks = breaks, nugget.tolerance = nugget.tolerance)
  uvec <- dbins$uvec ; bins.lim <- dbins$bins.lim
  nbins <- nc <- length(bins.lim) - 1
  if(is.null(max.dist)) max.dist <- max(bins.lim)
  min.dist <- min(cloud$u)
  if (!is.matrix(cloud$v)) {
    vbin <- rep(0, nc)
    nbin <- rep(0, nc)
    sdbin <- rep(0, nc)
    if (bin.cloud == TRUE | bin.cloud == 'diff') bins.clouds <- list()
    for (i in 1:nc) {
      ind <- (cloud$u > bins.lim[i]) & (cloud$u <= bins.lim[i+1])
      vbin[i] <- mean(cloud$v[ind])
      nbin[i] <- sum(ind)
      if (estimator.type == "modulus") 
        vbin[i] <- ((vbin[i])^4)/(0.914 + (0.988/nbin[i]))
      if (nbin[i] > 0) sdbin[i] <- sqrt(var(cloud$v[ind]))
      else sdbin[i] <- NA
      if (bin.cloud == TRUE) bins.clouds[[i]] <- cloud$v[ind]
      if (bin.cloud == 'diff') bins.clouds[[i]] <- cloud$d[ind]
      NULL
    }
    if (uvec[1] == 0) uvec[1] <- (bins.lim[1] + bins.lim[2])/2
    if (min.dist == 0) {
      ind <- (cloud$u == 0)
      n.zero <- sum(ind)
      v.zero <- mean(cloud$v[ind])
      if (bin.cloud == TRUE | bin.cloud == 'diff') {
        bins.clouds[2:(length(bins.clouds) + 1)] <- bins.clouds[1:nc]
        if(bin.cloud == 'diff') bins.clouds[[1]] <- cloud$d[ind]
        else bins.clouds[[1]] <- cloud$v[ind]
      }
      if (estimator.type == "modulus") 
        v.zero <- ((v.zero)^4)/(0.914 + (0.988/n.zero))
      if (n.zero > 0) 
        sd.zero <- sqrt(var(cloud$v[ind]))
      else sd.zero <- NA
      uvec <- c(0, uvec)
      vbin <- c(v.zero, vbin)
      nbin <- c(n.zero, nbin)
      sdbin <- c(sd.zero, sdbin)
    }
#    if(keep.NA == FALSE){
#      u <- uvec[!is.na(vbin)]
#      v <- vbin[!is.na(vbin)]
#      n <- nbin[!is.na(vbin)]
#      sd <- sdbin[!is.na(vbin)]
#    }
#    else{
      u <- uvec
      v <- vbin
      n <- nbin
      sd <- sdbin
#    }
    if (bin.cloud == TRUE | bin.cloud == 'diff') 
      bins.clouds <- bins.clouds[!is.na(vbin)]
  }
  else {
    if (bin.cloud == TRUE | bin.cloud == 'diff') 
      stop("option bins.cloud=TRUE allowed only for 1 variable")
    nvcols <- ncol(cloud$v)
    vbin <- matrix(0, nrow = nc, ncol = nvcols)
    nbin <- rep(0, nc)
    sdbin <- matrix(0, nrow = nc, ncol = nvcols)
    for (i in 1:nc) {
      ind <- (cloud$u >= bins.lim[i]) & (cloud$u < bins.lim[i+1])
      nbin[i] <- sum(ind)
      for (j in 1:nvcols) {
        vbin[i, j] <- mean(cloud$v[ind, j])
        if (estimator.type == "modulus") 
          vbin[i, j] <- ((vbin[i, j])^4)/(0.914 + (0.988/nbin[i]))
        if (nbin[i] > 0) 
          sdbin[i, j] <- sqrt(var(cloud$v[ind, j]))
        else sdbin[i, j] <- NA
      }
      NULL
    }
    if (uvec[1] == 0) 
      uvec[1] <- (bins.lim[1] + bins.lim[2])/2
    if (min.dist == 0) {
      v.zero <- rep(0, nvcols)
      n.zero <- rep(0, nvcols)
      sd.zero <- rep(0, nvcols)
      for (j in 1:nvcols) {
        ind <- (cloud$u == 0)
        n.zero[j] <- sum(ind)
        v.zero[j] <- mean(cloud$v[ind, j])
        if (estimator.type == "modulus") 
          v.zero[j] <- ((v.zero[j])^4)/(0.914 + (0.988/n.zero[j]))
        if (n.zero[j] > 0) 
          sd.zero[j] <- sqrt(var(cloud$v[ind, j]))
        else sd.zero[j] <- NA
        uvec <- c(0, uvec)
        vbin <- rbind(v.zero, vbin)
        nbin <- c(n.zero, nbin)
        sdbin <- rbind(sd.zero, sdbin)
      }
    }
    if(keep.NA == FALSE){
      u <- uvec[!is.na(vbin[, 1])]
      n <- nbin[!is.na(vbin[, 1])]
      v <- matrix(0, nrow = length(u), ncol = nvcols)
      sd <- matrix(0, nrow = length(u), ncol = nvcols)
    }
    else{
      u <- uvec
      n <- nbin
      v <- matrix(0, nrow = length(u), ncol = nvcols)
      sd <- matrix(0, nrow = length(u), ncol = nvcols)
    }
    for (j in 1:nvcols) {
      if(keep.NA == FALSE){
        v[, j] <- vbin[!is.na(vbin[, j]), j]
        sd[, j] <- sdbin[!is.na(vbin[, j]), j]
      }
      else{
        v[, j] <- vbin[, j]
        sd[, j] <- sdbin[, j]
      }
    }
  }
  if (nugget.tolerance > 1e-12) {
    u[1] <- 0
  }
  result <- list(u = u, v = v, n = n, sd = sd, bins.lim = bins.lim, output.type = "bin")
  if (!is.matrix(cloud$v) && (bin.cloud == TRUE | bin.cloud == 'diff'))
    result$bin.cloud <- bins.clouds
  ##  if (!is.null(class(cloud))) 
  if (length(class(cloud)) > 0) 
  oldClass(result) <- class(cloud)
  return(result)
}

"plot.variogram" <-
  function (x, max.dist, vario.col = "all", scaled = FALSE,  
            var.lines = FALSE,  envelope.obj = NULL,
            pts.range.cex, bin.cloud = FALSE,  ...) 
{
  fc.call <- match.call()
  fc.call$max.dist <- fc.call$vario.col <- fc.call$scaled <- 
    fc.call$pts.range.cex <-  fc.call$bin.cloud <-
      fc.call$envelope.obj <- NULL
  Ldots <- list(...)
  argsnames <- c(names(formals(plot.default)), names(par()))
  names(Ldots) <- argsnames[charmatch(names(Ldots), argsnames)]
  if(missing(max.dist)) max.dist <- max(x$u)
  if(is.null(Ldots$xlim)) fc.call$xlim <- c(0, max.dist) 
  if (bin.cloud == TRUE && all(is.na(x$bin.cloud))) 
    stop("plot.variogram: object must be a binned variogram with option bin.cloud=TRUE")
  if (bin.cloud == TRUE && any(!is.na(x$bin.cloud))) 
    boxplot(x$bin.cloud, varwidth = TRUE, 
            xlab = "distance", ylab = paste(x$estimator.type, "variogram"))
  else {
    if(!missing(pts.range.cex)){
      cex.min <- min(pts.range.cex)
      cex.max <- max(pts.range.cex)
      if(cex.min != cex.max){
        pts.prop <- TRUE
        sqn <- sqrt(x$n[x$u <= max.dist])
        pts.cex <- cex.min + ((sqn - min(sqn)) * (cex.max - cex.min) / (max(sqn) - min(sqn)))
      }
      else pts.prop <- FALSE
    }
    else pts.prop <- FALSE 
    u <- x$u[x$u <= max.dist]
    v <- x$v
    if(is.vector(v) | length(v) == length(x$u))
      v <- matrix(v, ncol=1)
    v <- v[x$u <= max.dist,, drop=FALSE]
    if(vario.col == "all")
      vario.col <- 1:dim(v)[2]
    else
      if(mode(vario.col) != "numeric" | any(vario.col > ncol(v)))
        stop("argument vario.col must be equals to \"all\" or a vector indicating the column numbers to be plotted")
    v <- v[, vario.col, drop=FALSE]
    if (scaled) v <- t(t(v)/x$var.mark[vario.col])
    if(is.null(list(...)$ylim)){ 
      ymax <- max(v)
      if (!is.null(envelope.obj)) 
        ymax <- max(c(envelope.obj$v.upper, ymax))
      fc.call$ylim <- c(0, ymax)
    }
    if(ncol(v) == 1){
      fc.call[[1]] <- as.name("plot")
      fc.call$x <- data.frame(distance=u, semivariance = as.vector(v))
      if(pts.prop) fc.call$cex <- pts.cex
      eval(fc.call, sys.frame(sys.parent()))
    }
    else{
      fc.call[[1]] <- as.name("matplot")
      fc.call$x <- u
      fc.call$y <- v
      if(is.null(Ldots$xlab)) fc.call$xlab <- "distance"
      if(is.null(Ldots$ylab)) fc.call$ylab <- "semivariance"
      if(is.null(Ldots$ty)){
        if (x$output.type == "bin") fc.call$type <- "p"
        if (x$output.type == "smooth") fc.call$type <- "l"
        if (x$output.type == "cloud") fc.call$type <- "p"
      }
#      if(is.null(Ldots$col)) fc.call$col <- 1:6
#      if(is.null(Ldots$lty)) fc.call$lty <- 1:5
#      if(is.null(Ldots$lwd)) $lwd <- 1
#      if(is.null(Ldots$pch)) Ldots$pch <- NULL
#      if(is.null(Ldots$cex)) Ldots$cex <- NULL
#      if(is.null(Ldots$add)) Ldots$add <- FALSE
#      
#      matplot(x=u, y= v, xlim = c(0, max.dist), ylim = Ldots$ylim, 
#              xlab = Ldots$xlab, ylab = Ldots$ylab, type = Ldots$type,
#              add = Ldots$add, pch = Ldots$pch,
#              lty = Ldots$lty, lwd = Ldots$lwd, col = Ldots$col)
      eval(fc.call, sys.frame(sys.parent()))
    }
    if (var.lines) {
      if (scaled) abline(h = 1, lty = 3)
      else abline(h = x$var.mark, lty = 3)
    }
    if (!is.null(envelope.obj)) {
      lines(u, envelope.obj$v.lower, lty = 4)
      lines(u, envelope.obj$v.upper, lty = 4)
    }
  }
  return(invisible())
}

"lines.variomodel" <-
  function (x, ...)
{
  UseMethod("lines.variomodel")
}

#"lines.variomodel.list" <-
#  function (x, ...)
#{
#do.call
#}

"lines.variomodel.default" <-
  function (x, cov.model, cov.pars, nugget, kappa,
            max.dist, scaled = FALSE, ...)
{
  if(missing(x)) x <- seq(0, max.dist, l=100)
  if(is.list(x))
    for(i in names(x)) assign(i, x[[i]])
  if(mode(x) == "numeric" && missing(max.dist)) max.dist <- max(x, na.rm=TRUE)
  ## reading and checking model/other components
  fn.env <- sys.frame(sys.nframe())
  .check.cov.model(cov.model=cov.model, cov.pars=cov.pars, kappa=kappa,
                   env=fn.env, output=FALSE)
  if(missing(nugget)){
    if(missing(x) || is.null(x$nugget)) 
      stop("argument nugget must be provided")
    else nugget <- x$nugget
  }
  if(missing(max.dist)){
    if(missing(x) || is.null(x$max.dist)) 
      stop("argument max.dist must be provided")
    else max.dist <- x$max.dist
  }
  ## computing the total sill for single or nested model 
  if (is.vector(cov.pars)) 
    sill.total <- nugget + cov.pars[1]
  else sill.total <- nugget + sum(cov.pars[, 1])
  ## checking whether to scale the variogram 
  if (scaled){
    if(is.vector(cov.model))
      cov.pars[1] <-  cov.pars[1]/sill.total
    else cov.pars[,1] <-  cov.pars[,1]/sill.total
    sill.total <- 1
  }
  ## defining a function to plot the variogram curve
  gamma.f <- function(x){
    if(any(cov.model == c("linear", "power"))){
      if(is.vector(cov.pars)){
        if(cov.model == "linear") cov.pars[2] <- 1
        return(nugget + cov.pars[1] * (x^cov.pars[2]))
      }
      else{
        if(length(cov.model) == 1) cov.model <- rep(cov.model, nrow(cov.pars))
        if(length(cov.model) != nrow(cov.pars)) stop("cov.model and cov.pars have incompatible dimentions")
        vals <- rep(nugget, length(x))
        for (i in 1:nrow(cov.pars)){
          if(any(cov.model[i] == c("linear", "power"))){
            if(cov.model[i] == "linear") cov.pars[i,2] <- 1
            vals <- vals + (cov.pars[i,1] * (x^cov.pars[i,2]))
          }
          else vals <- vals + cov.pars[i,1] - cov.spatial(x, cov.model = cov.model[i],
                                                          kappa = kappa[i], cov.pars = cov.pars[i,])
        }
        return(vals)
      }
    }
    else
      return(sill.total -
             cov.spatial(x, cov.model = cov.model,
                         kappa = kappa, cov.pars = cov.pars))
    
  }
  ## ploting the curve
  curve(gamma.f(x), from=0, to=max.dist, add=TRUE, ...)
  return(invisible())
}

"lines.variomodel.variofit" <-
  "lines.variomodel.likGRF" <-
  function (x, max.dist, scaled = FALSE, ...)
{
  my.l <- list()
  if(missing(max.dist)){
    my.l$max.dist <- x$max.dist
    if (is.null(my.l$max.dist)) 
      stop("argument max.dist needed for this object")
  }
  else
    my.l$max.dist <- max.dist
  if (any(x$cov.model == c("matern","powered.exponential",
            "cauchy", "gencauchy", "gneiting.matern"))) 
    my.l$kappa <- x$kappa
  else kappa <- NULL
  if (is.vector(x$cov.pars)) 
    my.l$sill.total <- x$nugget + x$cov.pars[1]
  else my.l$sill.total <- x$nugget + sum(x$cov.pars[, 1])
  my.l$nugget <- x$nugget
  my.l$cov.pars <- x$cov.pars
  my.l$cov.model <- x$cov.model
  if (scaled){
    if(is.vector(x$cov.model))
      my.l$cov.pars[1] <-  my.l$cov.pars[1]/my.l$sill.total
    else my.l$cov.pars[,1] <-  my.l$cov.cov.pars[,1]/my.l$sill.total
    my.l$sill.total <- 1
  }
  gamma.f <- function(x, my.l){
    if(any(my.l$cov.model == c("linear", "power")))
      return(my.l$nugget + my.l$cov.pars[1] * (x^my.l$cov.pars[2]))
    else
      return(my.l$sill.total -
        cov.spatial(x, cov.model = my.l$cov.model,
                    kappa = my.l$kappa,
                    cov.pars = my.l$cov.pars))
  }
  curve(gamma.f(x,my.l=my.l), from=0, to=my.l$max.dist, add=TRUE, ...)
  return(invisible())
}

"lines.variogram" <-
  function (x, max.dist, type = "o", scaled = FALSE, pts.range.cex, ...) 
{
  if(missing(max.dist)) max.dist <- max(x$u)
  if(!missing(pts.range.cex)){
    cex.min <- min(pts.range.cex)
    cex.max <- max(pts.range.cex)
    if(cex.min != cex.max){
      pts.prop <- TRUE
      sqn <- sqrt(x$n[x$u <= max.dist])
      pts.cex <- cex.min + ((sqn - min(sqn)) * (cex.max - cex.min) / (max(sqn) - min(sqn)))
    }
    else pts.prop <- FALSE
  }
  else pts.prop <- FALSE 
  if (scaled) 
    x$v <- x$v/x$var.mark
  if (!is.matrix(x$v)){
    if(pts.prop)
      lines(x$u[x$u <= max.dist], x$v[x$u <= max.dist], 
            type = type, cex = pts.cex, ...)
    else
      lines(x$u[x$u <= max.dist], x$v[x$u <= max.dist], 
            type = type, ...)
  }
  else {
    for (j in 1:ncol(x$v)){
      if(pts.prop)
        lines(x$u[x$u <= max.dist], 
              x$v[x$u <= max.dist, j], type = type, cex = pts.cex, ...)
      else
        lines(x$u[x$u <= max.dist], 
              x$v[x$u <= max.dist, j], type = type, ...)
    }
  }
  return(invisible())
}

"variog.mc.env" <-
  function (geodata, coords = geodata$coords, data = geodata$data,
            obj.variog, nsim = 99, save.sim = FALSE, messages) 
{
  call.fc <- match.call()
  if(missing(geodata)) geodata <- list(coords=coords, data=data)
  if(missing(messages))
    messages.screen <- as.logical(ifelse(is.null(getOption("geoR.messages")), TRUE, getOption("geoR.messages")))
  else messages.screen <- messages
  ##
  ## Checking input
  ##
  obj.variog$v <- NULL
  if((is.matrix(data) | is.data.frame(data)))
    if(ncol(data) > 1)
    	stop("envelops can be computed for only one data set at once")
  if(!is.null(obj.variog$estimator.type))
    estimator.type <- obj.variog$estimator.type
  else estimator.type <- "classical"
  ##
  ## transformation
  ##
  if (abs(obj.variog$lambda - 1) > 0.0001) {
    if (abs(obj.variog$lambda) < 0.0001) 
      data <- log(data)
    else data <- ((data^obj.variog$lambda) - 1)/obj.variog$lambda
  }
  ##
  ## trend removal
  ##
  xmat <- unclass(trend.spatial(trend = obj.variog$trend,
                                geodata = geodata))
  if (obj.variog$trend != "cte") {
    if (is.vector(data)) {
      data <- lm(data ~ xmat + 0)$residuals
      names(data) <- NULL
    }
    else {
      only.res <- function(y, x) {
        lm(y ~ xmat + 0)$residuals
      }
      data <- apply(data, 2, only.res, x = xmat)
    }
  }
  ##
  ## generating several "data-sets" by permutating data values
  ##
  if (messages.screen) 
    cat(paste("variog.env: generating", nsim,
              "simulations by permutating data values\n"))
  simula <- list(coords = coords)
  n.data <- length(data)
  perm.f <- function(i, data, n.data){return(data[sample(1:n.data)])}
  simula$data <- apply(as.matrix(1:nsim),1, perm.f, data=data, n.data=n.data)
  ##
  ## computing empirical variograms for the simulations
  ##
  if (messages.screen) 
    cat(paste("variog.env: computing the empirical variogram for the", 
              nsim, "simulations\n"))
  nbins <- length(obj.variog$bins.lim) - 1
  ##
  ##
  if(obj.variog$direction == "omnidirectional"){
    bin.f <- function(sim){
      cbin <- vbin <- sdbin <- rep(0, nbins)  
      temp <- .C("binit",
                 as.integer(obj.variog$n.data),
                 as.double(as.vector(coords[,1])),
                 as.double(as.vector(coords[,2])),
                 as.double(as.vector(sim)),
                 as.integer(nbins),
                 as.double(as.vector(obj.variog$bins.lim)),
                 as.integer(estimator.type == "modulus"),
                 as.double(max(obj.variog$u)),
                 as.double(cbin),
                 vbin = as.double(vbin),
                 as.integer(FALSE),
                 as.double(sdbin),
                 PACKAGE = "geoR")$vbin
      return(temp)
    }
    simula.bins <- apply(simula$data, 2, bin.f)
  }
  else{
    variog.vbin <- function(x, ...)
      variog(geodata = geodata, data = x, 
             uvec = obj.variog$uvec,
             estimator.type = obj.variog$estimator.type, 
             nugget.tolerance = obj.variog$nugget.tolerance, max.dist = obj.variog$max.dist,
             pairs.min = obj.variog$pairs.min,
             direction = obj.variog$direction, tolerance=obj.variog$tolerance,
             messages.screen = FALSE, ...)$v 
    simula.bins <- apply(simula$data, 2, variog.vbin)
  }
  simula.bins <- simula.bins[obj.variog$ind.bin,]
  if(save.sim == FALSE) simula$data <- NULL
  ##
  ## computing envelops
  ##
  if (messages.screen) 
    cat("variog.env: computing the envelops\n")
  limits <- apply(simula.bins, 1, range)
  res.env <- list(u = obj.variog$u, v.lower = limits[1, ],
                  v.upper = limits[2,])
  if(save.sim) res.env$simulations <- simula$data
  res.env$call <- call.fc
  oldClass(res.env) <- "variogram.envelope"
  return(res.env)
}

"lines.variogram.envelope" <-
  function(x, lty=3, ...)
{
  lines(x$u, x$v.lower, ...)
  lines(x$u, x$v.upper, ...)
  return(invisible())
}

".define.bins" <-
  function(max.dist, uvec = "default", breaks = "default", nugget.tolerance)
{
  if(all(breaks ==  "default")){
    if (all(uvec == "default")) uvec <- 13
    if (mode(uvec) == "numeric"){
      if(length(uvec) == 1){
        bins.lim <- seq(0, max.dist, l = uvec+1)
        bins.lim <- c(0, nugget.tolerance, bins.lim[bins.lim >  nugget.tolerance])
        uvec <- 0.5 * (bins.lim[-1] + bins.lim[-length(bins.lim)])
      }
      else{
        uvec <- c(0, uvec)
        nvec <- length(uvec)
        d <- 0.5 * diff(uvec[2:nvec])
        bins.lim <- c(0, (uvec[2:(nvec - 1)] + d), (d[nvec - 2] + uvec[nvec]))
        bins.lim <- c(0, nugget.tolerance, bins.lim[bins.lim >  nugget.tolerance])
      }
    }
    else stop("argument uvec can only take a numeric vector")
  }
  else{
    if(mode(breaks) != "numeric") stop("argument breaks can only take a numeric vector")
    else bins.lim <- breaks
    bins.lim <- c(0, nugget.tolerance, bins.lim[bins.lim >  nugget.tolerance])
    uvec <- 0.5 * (bins.lim[-1] + bins.lim[-length(bins.lim)])
  }
  return(list(uvec = uvec, bins.lim = bins.lim))
}
