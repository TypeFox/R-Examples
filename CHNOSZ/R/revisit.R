# CHNOSZ/revisit.R
# 20090415 functions related to diversity calculations
# 20100929 merged draw.diversity and revisit

optimal.index <- function(z, objective) {
  # for a vector, returns the index of the optimum value
  # for a matrix, returns the x, y coordinates of the optimum
  # find the minimum or maximum? look at attribute of objfun
  objfun <- get.objfun(objective)
  optimum <- attributes(objfun)$optimum
#  # do we care about the sign of the index?
#  if(tolower(target) %in% c("sd", "sd.log", "cv", "cv.log", "rmsd", "cvrmsd")) 
#    doabs <- TRUE else doabs <- FALSE
#  if(doabs) z <- abs(z)
  # the value of the optimum
  if(optimum=="minimal") optval <- z[which.min(z)]
  else optval <- z[which.max(z)]
  # the index of the optimum
  # (or indices if multiple instances of the optimum)
  ret.val <- which(z==optval, arr.ind=TRUE)
  return(ret.val)
}

extremes <- function(z, objective) {
  # are we interested in a maximum or minimum?
  objfun <- get.objfun(objective)
  optimum <- attributes(objfun)$optimum
#  # do we care about the sign of the index?
#  if(tolower(target) %in% c("sd", "sd.log", "cv", "cv.log", "rmsd", "cvrmsd")) 
#    doabs <- TRUE else doabs <- FALSE
#  if(doabs) z <- abs(z)
  # takes a matrix, returns the y as f(x) and x as f(y)
  # trajectories of the optimum
  y <- x <- numeric()
  xres <- ncol(z)
  yres <- nrow(z)
  if(optimum=="minimal") {
    for(i in 1:xres) y <- c(y, which.min(z[i,]))
    for(i in 1:yres) x <- c(x, which.min(z[,i]))
  } else {
    for(i in 1:xres) y <- c(y, which.max(z[i,]))
    for(i in 1:yres) x <- c(x, which.max(z[,i]))
  }
  # stop if we missed some
  if(length(x)!=xres) stop("optima not found for all y")
  if(length(y)!=yres) stop("optima not found for all x")
  return(list(x=x, y=y))
}

revisit <- function(eout, objective = "CV", loga2 = NULL, loga0 = NULL, ispecies = NULL,
  col = par("fg"), yline = 2, ylim = NULL, cex = par("cex"),
  lwd = par("lwd"), mar = NULL, side = 1:4, xlim = NULL, labcex = 0.6,
  pch = 1, main = NULL, plot.it = NULL, add = FALSE, plot.optval = TRUE,
  style.2D = "contour") {
  # given a list of logarithms of activities of species
  # (as vectors or matrices or higher dimensional arrays) 
  # calculate a diversity index or thermodynamic property 
  # with the same dimensions   20090316 jmd v1
  # eout can be the output from equilibrate (enables plotting)
  # or simply a list of logarithms of activity 
  # (each list entry must have the same dimensions)

  # if the entries have the same lengths they are assumed to be logarithms of activity
  ud <- unique(sapply(eout, length))
  if(length(ud)==1) {
    # eout is list of logarithms of activity
    if(missing(plot.it)) plot.it <- FALSE
    if(plot.it) stop("can't make a plot if 'eout' is not the output from equilibrate()")
    loga1 <- eout
    eout.is.eout <- FALSE
  } else {
    # test if eout is the output from equilibrate()
    if(!"loga.equil" %in% names(eout))
      stop(paste("the list provided in 'eout' is not the output from equilibrate()"))
    if(missing(plot.it)) plot.it <- TRUE
    loga1 <- eout$loga.equil
    eout.is.eout <- TRUE
  }

  # take a subset (or all) of the species
  if(is.null(ispecies)) ispecies <- 1:length(loga1)
  loga1 <- loga1[ispecies]
  # the dimensions 
  dim1 <- dim(as.array(loga1[[1]]))
  # the number of dimensions
  nd <- ifelse(identical(dim1, 1L), 0, length(dim1))
  msgout(paste("revisit: calculating", objective, "in", nd, "dimensions\n"))

  # get the objective function
  objfun <- get.objfun(objective)
  # the arguments to the function (a1, [a2]) or (loga1, [loga2], [Astar])
  objargs <- names(formals(objfun))

  # these objectives only depend on the activities (a1/loga1):
  # shannon, SD, CV, QQR
  if(any(grepl("a1", objargs))) {
    # vectorize the list entries: a1/loga1
    loga1 <- sapply(loga1, c)
    # for 0-D case we need to keep loga1 as a 1-row matrix (sapply simplifies to vector)
    if(nd==0) loga1 <- t(loga1)
    # convert infinite values to NA
    loga1[is.infinite(loga1)] <- NA
    # if we have loga0, calculate the base-2 log ratio (loga1/loga0)
    base <- 10
    if(!is.null(loga0)) {
      loga1 <- t(t(loga1) - loga0) * log2(10)
      base <- 2
    }
    # remove logarithms if necessary
    if(any(grepl("loga1", objargs))) a1 <- loga1
    else a1 <- base^loga1
  }

  # these objectives also depend on reference activities (a2/loga2):
  # RMSD, CVRMSD, spearman, pearson, DGtr
  if(any(grepl("a2", objargs))) {
    # check that all needed arguments are present
    if(is.null(loga2)) stop(paste("loga2 must be supplied for", objective))
    # if loga2 is a single value, expand it into the dimensions of loga1
    if(length(loga2)==1) loga2 <- rep(loga2, length.out=ncol(loga1))
    # check that loga2 has the same length as loga1
    if(!identical(ncol(loga1), length(loga2))) stop(paste("loga2 has different length (", 
      length(loga2), ") than list in eout (", ncol(loga1), ")", sep=""))
    # remove logarithms if necessary
    if(any(grepl("loga2", objargs))) a2 <- loga2
    else a2 <- base^loga2
  }


  # construct array of values: Astar (for DGtr)
  if(any(grepl("Astar", objargs))) {
    Astar <- eout$Astar[ispecies]
    # one row for each condition
    Astar <- sapply(Astar, as.vector)
    # for 0-D case we want a 1-row matrix (sapply simplifies to vector)
    if(nd==0) Astar <- t(Astar)
  }

  # calculation of the objective function
  # the symbol "H" is reminiscent of the first implemented target, shannon entropy
  if(length(objargs) == 1) H <- objfun(a1)
  else if(length(objargs) == 2) H <- objfun(a1, a2)
  else if(length(objargs) == 3) H <- objfun(a1, a2, Astar)

  # replace dims
  dim(H) <- dim1

  ## now on to plotting + assembling return values
  # for zero or more than two dimensions we'll just return the values
  # of the objective function and the index of the optimal value
  iopt <- optimal.index(H, objective)
  ret.val <- list(H=H, iopt=iopt)
  # get information about the x-axis
  if(eout.is.eout & nd > 0) {
    xname <- eout$vars[1]
    # the x-values
    xs <- eout$vals[[1]]
    xrange <- range(xs)
  } else xs <- NA

  # make plots and assemble return values
  if(nd==0) {
    # a 0-D plot
    if(plot.it) {
      if(objective=="qqr") {
        # make a q-q plot for qqr
        qqnorm(loga1, col=col, pch=pch, main=NA)
        qqline(loga1)
      } else if(any(grepl("a2", objargs))) {
        # plot the points for a referenced objective
        ylab <- "loga1"
        xlab <- "loga2"
        if(is.null(xlim)) xlim <- extendrange(loga2)
        if(is.null(ylim)) ylim <- extendrange(loga1)
        plot(loga2, loga1, xlab=xlab, ylab=ylab, pch=pch, col=col, xlim=xlim, ylim=ylim)
        # add a 1:1 line
        lines(range(loga2), range(loga2), col="grey")
        # add a lowess line
        if(!is.null(lwd)) {
          ls <- loess.smooth(loga2, loga1, family="gaussian")
          lines(ls$x, ls$y, col="red", lwd=lwd)
        }
      } else plot.it <- FALSE
      # add a title
      if(missing(main)) main <- paste(objective, "=", round(H,3)) 
      if(plot.it) title(main=main)
    }
  } else if(nd==1) {
    # locate the optimal value
    ixopt <- c(iopt)
    xopt <- xs[ixopt]
    optimum <- H[ixopt]
    ret.val <- list(H=H, ixopt=ixopt, xopt=xopt, optimum=optimum)
    # a 1-D plot
    if(plot.it) {
      if(is.null(ylim)) ylim <- extendrange(H, f=0.075)
      if(is.null(xlim)) xlim <- xrange
      # format the objective name if it's DGtr
      if(objective=="DGtr") ylab <- expr.property("DGtr/2.303RT")
      else ylab <- objective
      if(!add) thermo.plot.new(xlim=xlim, ylim=ylim, xlab=axis.label(xname),
        ylab=ylab, yline=yline, cex=cex, lwd=lwd, mar=mar, side=side)
      # plot the values
      lines(xs, c(H), col=col)
      # indicate the optimal value
      if(plot.optval) abline(v=xopt, lty=2)
    } 
  } else if(nd==2) {
    # a 2-D plot
    # information about the y-axis
    if(eout.is.eout) {
      yname <- eout$vars[2]
      ys <- eout$vals[[2]]
      yrange <- range(ys)
    } else ys <- NA
    # locate the optimal value
    ixopt <- iopt[, 1]
    iyopt <- iopt[, 2]
    optimum <- H[ixopt, iyopt]
    ret.val <- list(H=H, ixopt=ixopt, iyopt=iyopt, xopt=xs[ixopt], yopt=ys[iyopt], optimum=optimum) 
    if(plot.it) {
      # start the plot
      if(is.null(xlim)) xlim <- xrange
      if(is.null(ylim)) ylim <- yrange
      if(!add) thermo.plot.new(xlim=xlim, ylim=ylim, xlab=axis.label(xname), 
        ylab=axis.label(yname), yline=yline, side=side, cex=cex, mar=mar)
      if(style.2D=="contour") contour(xs, ys, H, add=TRUE, labcex=labcex)
      else if(style.2D=="image") {
        # get colors based on number of species
        nspecies <- length(loga1)
        if(missing(col)) col <- heat.colors(nspecies)
        image(xs, ys, H, add=TRUE, col=col)
      } else stop(paste("2D plot style", style.2D, "not one of 'contour' or 'image'"))
      if(plot.optval) {
        # plot the location(s) of the extremum
        points(xs[ixopt], ys[iyopt], pch=8, cex=2)
        # show trajectories of the extrema
        iexts <- extremes(H, objective)
        # take out large jumps
        yext <- ys[iexts$y]
        yext.1 <- c(yext[2:length(yext)], yext[length(yext)])
        yext.2 <- c(yext[1], yext[1:length(yext)-1])
        yext[abs(yext.1-yext)/abs(diff(range(ys))) > 0.1] <- NA
        yext[abs(yext.2-yext)/abs(diff(range(ys))) > 0.1] <- NA
        lines(xs, yext, lty=3, col="blue")
        xext <- xs[iexts$x]
        xext.1 <- c(xext[2:length(xext)], xext[length(xext)])
        xext.2 <- c(xext[1], xext[1:length(xext)-1])
        xext[abs(xext.1-xext)/abs(diff(range(xs))) > 0.1] <- NA
        xext[abs(xext.2-xext)/abs(diff(range(xs))) > 0.1] <- NA
        lines(xext, ys, lty=3, col="seagreen")
      }
    }
  }

  # return the results
  return(invisible(ret.val))
}


