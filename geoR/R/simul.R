## comments are Previous version before changes in RandomFields 2.*
"geoR2RF" <-
  function (cov.model, cov.pars, nugget = 0, kappa, aniso.pars)
{
  cov.model <- match.arg(cov.model, choices =  geoRCovModels)
  if(missing(aniso.pars)) aniso.pars <- NULL
  if(missing(kappa)) kappa <- NULL
  if (length(cov.pars) != 2)
    stop("cov.pars must be an vector of size 2 with values for the parameters sigmasq and phi")
  RFmodel <- switch(cov.model, matern = "whittle", exponential = "exponential",
                    gaussian = "gauss", spherical = "spherical", circular = "circular",
                    cubic = "cubic", wave = "wave", power = "not compatible",
                    powered.exponential = "stable", cauchy = "cauchy",
                    gencauchy = "gencauchy", gneiting = "gneiting",
                    gneiting.matern = "not compatible", pure.nugget = "nugget")
  if (RFmodel == "not compatible") {
    warning("geoR cov.model not compatible with RandomFields model")
    return(RFmodel)
  }
  if (any(RFmodel == "gencauchy")) kappa <- rev(kappa)
  if (!any(RFmodel == c("gencauchy", "whittle", "stable"))) kappa <- NULL
#  if(is.null(aniso.pars))
#    return(list(list(model=RFmodel, var=cov.pars[1], kappa=kappa, scale=cov.pars[2]),
#                "+",
#                list(model="nugget", var=nugget)))
#    else{
#    mat <- solve(matrix(c(cos(aniso.pars[1]), sin(aniso.pars[1]),
#                          -sin(aniso.pars[1]), cos(aniso.pars[1])), ncol=2)) %*%
#                            diag(c(aniso.pars[2], 1)/cov.pars[2])
#    return(list(list(model=RFmodel, var=cov.pars[1], kappa=kappa, aniso=mat), 
#                "+",
#                list(model="nugget", var=nugget, aniso=diag(1,2))))
  if(is.null(aniso.pars))
    model <- list("+",
                list("$",  var=cov.pars[1], scale=cov.pars[2],
                if (is.null(kappa)) list(RFmodel) 
 		else list(RFmodel, k=kappa)), 
                list("$", var=nugget, list("nugget")))
  else{
    mat <- solve(matrix(c(cos(aniso.pars[1]), sin(aniso.pars[1]),
                          -sin(aniso.pars[1]), cos(aniso.pars[1])), ncol=2)) %*%
                            diag(c(aniso.pars[2], 1)/cov.pars[2])
    model <- list("+",
                list("$", var=cov.pars[1], aniso=mat,
                if (is.null(kappa)) list(RFmodel) 
                else
                list(RFmodel, k=kappa)),
                list("$", var=nugget, list("nugget")))
  }
#  Print(model)
  return(model)
}

"grf" <-
  function (n, grid = "irreg", nx, ny, xlims = c(0, 1), ylims = c(0, 1),
            borders, nsim = 1, cov.model = "matern",
            cov.pars = stop("missing covariance parameters sigmasq and phi"), 
            kappa = 0.5, nugget = 0, lambda = 1, aniso.pars = NULL, mean = 0, 
            method, RF = TRUE, messages) 
{
  ## checking input options
  call.fc <- match.call()
  if (missing(messages)) 
    messages.screen <- as.logical(ifelse(is.null(getOption("geoR.messages")), 
                                         TRUE, getOption("geoR.messages")))
  else messages.screen <- messages
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    warning(".Random.seed not initialised. Creating it with by calling runif(1)")
    runif(1)
  }
  rseed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  ## model setup
  cov.model <- match.arg(cov.model, choices =  geoRCovModels)
  if (cov.model == "stable") cov.model <- "powered.exponential"
  if (cov.model == "matern" && kappa == 0.5) cov.model <- "exponential"
  tausq <- nugget
  if (is.vector(cov.pars)) {
    sigmasq <- cov.pars[1]
    phi <- cov.pars[2]
    nst <- 1
  }
  else {
    sigmasq <- cov.pars[, 1]
    phi <- cov.pars[, 2]
    nst <- nrow(cov.pars)
  }
  sill.total <- tausq + sum(sigmasq)
  messa <- .grf.aux1(nst, nugget, sigmasq, phi, kappa, cov.model)
  ## 1 or 2-D simulation
  results <- list(coords=NULL, data=NULL)
  ##
  if ((!missing(nx) && nx == 1) | (!missing(ny) && ny == 1) | 
      diff(xlims) == 0 | diff(ylims) == 0) {
    sim1d <- TRUE
    if (messages.screen) cat("simulations in 1D\n")
  }
  else sim1d <- FALSE
  ## setting locations for simulation 
  if(mode(grid) == "character"){
    grid <- match.arg(grid, choices = c("irreg", "reg"))
    attr(results, "grid") <- grid    
    if(!missing(borders) & !sim1d){
     #if(!require(splancs, quietly=TRUE))
     #   stop("package splancs is required to simulate within borders provided by the user")
      results$borders <- borders
      if(grid == "irreg") results$coords  <- splancs::csr(poly=borders, npoints=n)
      else{
        if(!missing(nx) && !missing(ny)){
          bb <- bbox(borders)
          results$coords  <- splancs::gridpts(poly=borders, xs = diff(bb[1,])/(nx-1),
                                              ys=diff(bb[1,])/(nx-1))
        }
        else results$coords  <- splancs::gridpts(poly=borders, npts=n)
        xgrid <- round(sort(unique(results$coords[,1])), digits=12)
        ygrid <- round(sort(unique(results$coords[,2])), digits=12)
        attr(results, "xgrid") <- c(range(xgrid),unique(diff(xgrid)))
        attr(results, "ygrid") <- c(range(ygrid),unique(diff(ygrid)))
        names(attr(results, "xgrid")) <- c("xmin", "xmax", "xstep")
        names(attr(results, "ygrid")) <- c("ymin", "ymax", "ystep")
      }
    }
  }
  else{
    results$coords <- as.matrix(grid)
  	x1vals <- sort(unique(round(results$coords[,1], digits=12)))
  	x2vals <- sort(unique(round(results$coords[,2], digits=12)))
  	if(length(unique(diff(x1vals)) == 1) && length(unique(diff(x2vals)) == 1)) 
	    attr(results, "grid") <- "reg"
    if (messages.screen) cat("grf: simulation on a set of locations provided by the user\n")
  }
  if (!is.matrix(results$coords) & !is.data.frame(results$coords)) {
    ## isto aqui nunca eh chamado e pode ser eliminado???
    if (missing(nx)) {
      if (sim1d) nx <- ifelse(diff(xlims) == 0, 1, n)
      else nx <- ifelse((mode(grid)=="character" && grid=="reg"), round(sqrt(n)), n)
    }
    if (missing(ny)) {
      if (sim1d) ny <- ifelse(diff(ylims) == 0, 1, n)
      else ny <- ifelse((mode(grid)=="character" && grid=="reg"), round(sqrt(n)), n)
    }
    if (mode(grid) == "character" && grid == "irreg") {
        results$coords <- cbind(x = runif(nx, xlims[1], xlims[2]), 
                                y = runif(ny, ylims[1], ylims[2]))
      if (messages.screen) 
        cat(paste("grf: simulation(s) on randomly chosen locations with ", 
                  n, " points\n"))
    }
    else {
      xpts <- seq(xlims[1], xlims[2], length = nx)
      ypts <- seq(ylims[1], ylims[2], length = ny)
      xspacing <- ifelse(length(xpts) == 1, 0, diff(xpts[1:2])) 
      yspacing <- ifelse(length(ypts) == 1, 0, diff(ypts[1:2]))
      results$coords <- as.matrix(expand.grid(x = xpts, y = ypts))
      equal.spacing <- ifelse(abs(xspacing - yspacing) < 1e-12, TRUE, FALSE)
      if (messages.screen) 
        cat(paste("grf: generating grid ", nx, " * ", 
                  ny, " with ", (nx * ny), " points\n"))
      attr(results, "xgrid") <- c(xmin = xlims[1], xmax = xlims[2], xstep = xspacing)
      attr(results, "ygrid") <- c(ymin = ylims[1], ymax = ylims[2], ystep = yspacing)
    }
    if(!sim1d & missing(borders)){
      lbor <- as.matrix(expand.grid(xlims, ylims))
      results$borders <- lbor[chull(lbor),]
    }
  }
  n <- nrow(results$coords)
  if (length(unique(round(results$coords[, 1], digits = 12))) == 1 |
      length(unique(round(results$coords[, 2], digits = 12))) == 1) 
    sim1d <- TRUE
  else sim1d <- FALSE
  if (!RF && !is.null(aniso.pars)) {
    if (length(aniso.pars) != 2 | mode(aniso.pars) != "numeric") 
      stop("anisotropy parameters must be provided as a numeric vector with two elements: the rotation angle (in radians) and the anisotropy ratio (a number greater than 1)")
    if (messages.screen) 
      cat("grf: transforming to the isotropic space \n")
    results$coords <- coords.aniso(coords = results$coords, 
                                   aniso.pars = aniso.pars)
  }
  if (missing(method)) {
    method <- "cholesky"
  #  if (n > 500 && RF && require(RandomFields, quietly=TRUE)) 
    if (n > 500 && RF) 
      method <- "RF"
  }
  method <- match.arg(method, choices = c("cholesky", "svd", 
                                "eigen", "RF", "circular.embedding"))
  #if (method == "circular.embedding") {
  #  if (require(RandomFields, quietly=TRUE)) {
  #    method <- "RF"
  #    if (messages.screen) 
  #      warning("method \"circular.embedding\" now uses algorithm from the package RandomFields")
  #  }
  #  else stop("Option for method \"circular.embedding\" requires the instalation of the package RandomFields")
  #}
  if (messages.screen) {
    cat(messa$nst)
    cat(messa$nugget)
    cat(messa$cov.structures)
    if (method == "RF") 
      cat("grf: simulation using the function GaussRF from package RandomFields \n")
    else cat(paste("grf: decomposition algorithm used is: ", 
                   method, "\n"))
  }
  if (all(phi == 0)) 
    results$data <- matrix(rnorm((n * nsim), mean = 0, sd = sqrt(sill.total)), 
                           nrow = n, ncol = nsim)
  else {
    if (method == "RF") {
      RandomFields::RFoldstyle(old=TRUE)
      #require(RandomFields, quietly=TRUE)
      assign("setRF", geoR2RF(cov.model = cov.model, cov.pars = cov.pars, 
                       nugget = nugget, kappa = kappa, aniso.pars=aniso.pars), pos=1)
      if (!exists("xpts") || is.null(xpts)){
        results$data <- RandomFields::GaussRF(x = results$coords[, 1],y = results$coords[, 2],
                                model = get("setRF", pos=1), grid = FALSE, n = nsim)
      }
      else{
        results$data <- drop(matrix(RandomFields::GaussRF(x = xpts, y = ypts, model = get("setRF", pos=1),
                                            grid = TRUE, n = nsim), ncol = nsim))
      }
    }
    else
      results$data <- drop(crossprod(varcov.spatial(coords = results$coords, 
                                                    cov.model = cov.model, kappa = kappa,
                                                    nugget = nugget, 
                                                    cov.pars = cov.pars,
                                                    only.decomposition = TRUE,
                                                    func.inv = method)$sqrt.varcov, 
                                     matrix(rnorm((n * nsim)), nrow = n, ncol = nsim)))
  }
  if (length(mean) != 1 & length(mean) != dim(as.matrix(results$data))[1] & length(mean) != length(results$data)) 
    stop("the mean must be a scalar or a vector of the same size as the data")
  results$data <- results$data + mean
  if (lambda != 1) {
    if (lambda != 0) 
      results$data <- (results$data * lambda + 1)^(1/lambda)
    else results$data <- exp(results$data)
    messa$transformation <- paste("grf: Data transformed (Box-Cox), for lambda =", 
                                  lambda)
    if (messages.screen) 
      cat(messa$transformation)
    cat("\n")
  }
  if (!RF && !is.null(aniso.pars)) {
    if (messages.screen) 
      cat("grf: back-transforming to the anisotropic space \n")
    results$coords <- coords.aniso(coords = results$coords, 
                                   aniso.pars = aniso.pars, reverse = TRUE)
  }
  else {
    aniso.pars <- "no anisotropy parameters provided/used"
  }
  if (messages.screen) 
    cat(paste("grf: End of simulation procedure. Number of realizations:", 
              nsim, "\n"))
  results[c("cov.model", "nugget", "cov.pars", "kappa", "lambda", "aniso.pars", "method",
            ".Random.seed", "messages", "call")] <- 
              list(cov.model, nugget, cov.pars, kappa = kappa, lambda = lambda, aniso.pars, 
                   method, rseed,  messa, call.fc)
  attr(results, "borders") <- call.fc$borders
  attr(results, "sp.dim") <- ifelse(sim1d, "1d", "2d")
  oldClass(results) <- c("grf", "geodata", "variomodel")
  return(results)
}

              
".grf.aux1" <-
  function (nst, nugget, sigmasq, phi, kappa, cov.model) 
{
  cov.nst <- paste("grf: process with ", nst, " covariance structure(s)\n")
  cov.nugget <- paste("grf: nugget effect is: tausq=", nugget,"\n")
  cov.message <- NULL
  for (i in 1:nst) {
    if (phi[i] == 0) 
      cov.message[i] <- paste("grf: covariance model", i, "is a pure nugget effect\n")
    else {
      if(any(cov.model == c("matern","powered.exponential", 
          "cauchy", "gencauchy", "gneiting.matern"))) 
        cov.message[i] <- paste("grf: covariance model ", 
                                i, " is: ", cov.model, "(sigmasq=", sigmasq[i], 
                                ", phi=", phi[i], ", kappa = ", kappa, ")\n", sep = "")
      else cov.message[i] <- paste("grf: covariance model ", 
                                   i, " is: ", cov.model, "(sigmasq=", sigmasq[i], 
                                   ", phi=", phi[i], ")\n", sep = "")
    }
  }
  return(list(nst = cov.nst, nugget = cov.nugget, cov.structures = cov.message))
}

"lines.variomodel.grf" <-
  function (x, max.dist, n = 100, ...) 
{
  if(missing(max.dist)) max.dist <- max(dist(x$coords))
  distance <- seq(0, max.dist, length = n)
  if (is.vector(x$cov.pars)) 
    sill.total <- x$nugget + x$cov.pars[1]
  else sill.total <- x$nugget + sum(x$cov.pars[, 1])
  gamma <- sill.total - cov.spatial(distance, cov.model = x$cov.model, 
                                  kappa = x$kappa, cov.pars = x$cov.pars)
  lines(distance, gamma, ...)
  return(invisible())
}

".geoR_fullGrid" <-
  function(x, borders)
{
  xgrid <- seq(attr(x, "xgrid")["xmin"], attr(x, "xgrid")["xmax"],
               by=attr(x, "xgrid")["xstep"])
  ygrid <- seq(attr(x, "ygrid")["ymin"], attr(x, "ygrid")["ymax"],
               by=attr(x, "ygrid")["ystep"])
  data2plot <- rep(NA, length(xgrid)*length(ygrid))
  data2plot[.geoR_inout(expand.grid(xgrid, ygrid), borders)] <- x$data
  return(data2plot)
}

"image.grf" <-
  function (x, sim.number = 1, borders, x.leg, y.leg, ...) 
{
  ##
  ## this seems to cause problems overlapping maps
  ##op <- par(no.readonly=TRUE)
  ##on.exit(par(op))
  ##
  x1vals <- sort(unique(round(x$coords[,1], digits=12)))
  x2vals <- sort(unique(round(x$coords[,2], digits=12)))
  #if(length(unique(x1vals) == 1) && length(unique(x2vals) == 1)) reggrid <- TRUE
  nx <- length(x1vals)
  ny <- length(x2vals)
  ldots <- match.call(expand.dots = FALSE)$...
  if(is.vector(x$data)){
    if(sim.number != 1) stop("there is just one simulation in the object")
  }
  else
    x$data <- x$data[,sim.number]
  n <- length(x$data)
  ##
  ## Plotting simulations in 1-D
  ##
  if(attr(x, 'sp.dim') == "1d" | nx == 1 | ny == 1)
    do.call("plot.1d", c(list(x = x, x1vals = x1vals),
                         .ldots.set(ldots, type="plot.1d", data="simulation")))
  else{
    ##
    ## Plotting simulations in 2-D
    ##
    ## Checking for rectangular grid
    ##
    if(is.null(attr(x, "grid")) || attr(x, "grid") != "reg")
      stop("cannot produce image plot, probably you've got an irregular grid of locations")
    ##
    ## Preparing image plot elements
    ##
    if(missing(borders))
      borders <-  eval(attr(x, "borders"), envir= attr(x, "parent.env"))
    if (!is.null(borders)){
      if(nx * ny == n)
        x$data[!.geoR_inout(x$coords, borders)] <- NA
      else
        x$data <- .geoR_fullGrid(x = x, borders = borders)
    }
    do.call("image", c(list(x=x1vals, y=x2vals, z=matrix(x$data, ncol=ny)),
                       .ldots.set(ldots, type="image", data="simulation")))
    ##
    ## Adding the legend (if the case)
    ##
    if(!missing(x.leg) && !missing(y.leg)){
      if(is.null(ldots$col)) ldots$col <- heat.colors(12)
      do.call("legend.krige", c(list(x.leg=x.leg, y.leg=y.leg,
                                     values = x$data), ldots))
    }
    if(!is.null(borders)) polygon(borders)
  }
  return(invisible())
}

"persp.grf" <- 
  function(x, sim.number = 1, borders, ...)
{
  x1vals <- unique(round(x$coords[,1], digits=12))
  x2vals <- unique(round(x$coords[,2], digits=12))
  nx <- length(x1vals)
  ny <- length(x2vals)
  ldots <- match.call(expand.dots = FALSE)$...
  if(is.vector(x$data)){
    if(sim.number != 1) stop("there is just one simulation in the object")
  }
  else
    x$data <- x$data[,sim.number]
  n <- length(x$data)
  if(attr(x, 'sp.dim') == "1d" | nx == 1 | ny == 1)
    do.call("plot.1d", c(list(x = x, x1vals = x1vals),
                       .ldots.set(ldots, type="plot.1d", data="simulation")))
  else{
    if(is.null(attr(x, "grid")) || attr(x, "grid") != "reg")
      stop("cannot produce persp plot, probably you've got an irregular grid of locations")
    if(missing(borders))
      borders <-  eval(attr(x, "borders"), envir= attr(x, "parent.env"))
    if (!is.null(borders)){
      if(nx * ny == n)
        x$data[!.geoR_inout(x$coords, borders)] <- NA
      else
        x$data <- .geoR_fullGrid(x = x, borders = borders)
    }
    do.call("persp", c(list(x=x1vals, y=x2vals, z=matrix(x$data, ncol = ny)),
                       .ldots.set(ldots, type="persp", data="simulation")))
  }
  return(invisible())
}

"contour.grf" <- 
  function(x, sim.number = 1, borders, filled = FALSE, ...)
{
  x1vals <- unique(round(x$coords[,1], digits=12))
  x2vals <- unique(round(x$coords[,2], digits=12))
  nx <- length(x1vals)
  ny <- length(x2vals)
  ldots <- match.call(expand.dots = FALSE)$...
  if(is.vector(x$data)){
    if(sim.number != 1) stop("there is just one simulation in the object")
  }
  else
    x$data <- x$data[,sim.number]
  n <- length(x$data)
  if(attr(x, 'sp.dim') == "1d" | nx == 1 | ny == 1)
    do.call("plot.1d", c(list(x = x, x1vals = x1vals),
                       .ldots.set(ldots, type="plot.1d", data="simulation")))
  else{
    if(is.null(attr(x, "grid")) || attr(x, "grid") != "reg")
      stop("cannot produce contour plot, probably you've got an irregular grid of locations")
    if(missing(borders))
      borders <-  eval(attr(x, "borders"), envir= attr(x, "parent.env"))
    if (!is.null(borders)){
      if(nx * ny == n)
        x$data[!.geoR_inout(x$coords, borders)] <- NA
      else
        x$data <- .geoR_fullGrid(x = x, borders = borders)
    }
    if(filled)
      ldots.contour <- .ldots.set(ldots, type="filled.contour",
                                 data="prediction")
    else
      ldots.contour <- .ldots.set(ldots, type="contour",
                                 data="prediction")
    if(filled){
      if(is.null(ldots.contour$plot.axes)){
        ldots.contour$plot.axes <- quote({
          axis(1)
          axis(2)
          if(!is.null(coords.data)) points(coords.data, pch=20)
          if(!is.null(borders)) polygon(borders, lwd=2)
        })
      }
      do.call("filled.contour", c(list(x=x1vals, y=x2vals, z=matrix(x$data, ncol = ny)),
                                  ldots.contour))
    }
    else{
      do.call("contour", c(list(x=x1vals, y=x2vals, z=matrix(x$data, ncol = ny)),
                           ldots.contour))
    }
    if(!is.null(borders)) polygon(borders)
  }
  return(invisible())
}

"plot.grf" <-
  function (x, model.line = TRUE, plot.locations = FALSE, ...) 
{
  nsim <- ncol(x$data)
  if (plot.locations){
    points.geodata(x, pt.divide="equal", xlab = "Coord X", ylab = "Coord Y")
    if(is.null(list(...)$ask)){
      ask.now <- par()$ask
      par(ask = TRUE)
      on.exit(par(ask=ask.now)) 
    }
  }
  if (is.vector(x$cov.pars)) 
    sill.total <- x$nugget + x$cov.pars[1]
  else sill.total <- x$nugget + sum(x$cov.pars[, 1])
  if (x$lambda != 1){
    if (x$lambda == 0) data <- log(x$data)
    else data <- ((x$data^x$lambda)-1)/x$lambda
  }
  else
    data <- x$data          
  sim.bin <- variog(x, data=data)
  ldots <- list(...)
  if(is.null(ldots$ylim))
    plot(sim.bin, ylim=c(0, max(c(sill.total, sim.bin$v))),...)
  else
    plot(sim.bin, ...)    
  if (model.line){
    model <- list(nugget = x$nugget, cov.pars = x$cov.pars, 
                  kappa = x$kappa, max.dist = max(sim.bin$u),
                  cov.model = x$cov.model)
    if(is.null(ldots$lty))
      lines.variomodel(model, lty=2)
    else
      lines.variomodel(model, ...)
  }
  return(invisible())
}

"print.grf" <-
  function(x, ...)
{
  print.default(x, ...)
}

"lines.grf" <- function(x, ...){
  if(attr(x, "sp.dim") != "1d")
    stop("can only be used for simulations in  1-D")
  if(is.matrix(x$data))
    matplot(x$coords[,1], x$data, add=TRUE, ...)
  else
    lines(x$coords[,1], x$data, ...)
  return(invisible())
}
