plotGOF.LR <- function(
  x,
  beta.subset = "all",
  main = "Graphical Model Check",
  xlab,
  ylab,
  tlab = "item",
  xlim,
  ylim,
  type = "p",
  pos = 4,
  conf = NULL,
  ctrline = NULL,
  asp = 1,
  x_axis = TRUE,
  y_axis = TRUE,
  set_par = TRUE,
  reset_par = TRUE,
  ...
){
# graphical model check
# beta.subset...plot only a subset of beta-parameters; either "all" or an index vector
# x...object of class LR (from LRtest)
# tlab ... labelling: "item" abbreviated beta parameter name, "number" number from beta par list,
#            "identify" interactive, "none"
# pos ... (where the textlabel appears)
# conf ... confidence ellipses: NULL or
#              list(gamma=0.95, col="red", ia=TRUE, lty="dashed", which=all items in beta.subset)
# ctrline ... control lines (confidence bands): NULL or list(gamma=0.95,lty="solid", col="blue")
# ...     additional graphic parameters

  # save current options() and par() values and restore them on exit
  old_options <- options(locatorBell = FALSE)
  if(set_par){ old_par <- par(mar=c(4,4,3,0)+.5, no.readonly = TRUE) }
  on.exit({
    if(set_par && reset_par){ par(old_par) }
    options(old_options)
  })

  if(length(x$likgroup) > 2L) warning("Only the parameters for the first two subgroups are plotted!")

  if(missing(xlab)) xlab <- paste0("Beta for Group: ", x$spl.gr[1L])
  if(missing(ylab)) ylab <- paste0("Beta for Group: ", x$spl.gr[2L])

  nparg1 <- length(x$betalist[[1L]])
  nparg2 <- length(x$betalist[[2L]])
  if(nparg1 != nparg2) stop("Unequal number of parameters in the subgroups! Plot cannot be produced, choose another split in LRtest!")



  beta1 <- -x$betalist[[1L]] # -1 to obtain difficulty parameters
  beta2 <- -x$betalist[[2L]]

  if(is.character(beta.subset)) {
    if(beta.subset == "all"){
      beta.subset <- seq_along(beta1)
      #textlab <- names(beta1)
      switch(EXPR = tlab,
        item     = textlab <- gsub("^beta\\ I", "", names(beta1)), # remove "beta I" from names
        number   = textlab <- seq_along(beta1),
        identify = labs    <- gsub("^beta\\ I", "", names(beta1))
      )
    } else {
      textlab <- beta.subset
    }
  } else {
    switch(EXPR = tlab,
      item     = textlab <- gsub("^beta\\ I", "", names(beta1)[beta.subset]), # remove "beta I" from names
      number   = textlab <- beta.subset,
      identify = labs    <- gsub("^beta\\ I", "", names(beta1)[beta.subset])
    )
  }

# se's needed for ellipses and control lines
  if(is.list(conf) || is.list(ctrline)){
    if(any(is.na(unlist(x$selist)))) {
      warning('Confidence ellipses or control lines cannot be plotted.\n  LR object without standard errors.\n  Use option "se = TRUE" in LRtest()')
      conf <- ctrline <- NULL
    } else {
      s1 <- x$selist[[1L]]
      s2 <- x$selist[[2L]]
      v1 <- s1^2
      v2 <- s2^2
      suspicious.se <- any(cbind(s1, s2)[beta.subset] > 10)
      if(suspicious.se) warning("Suspicious size of standard error(s).\n  Check model specification, split criterion, data.")
    }

    #if(any(abs(cbind(beta1,beta2)[beta.subset])>8)){
    #   warning("Suspicious size of parameter estimate(s).\n  Check model specification, split criterion, data.")
    #   if(is.null(conf)) conf$ia <- FALSE
  }

###   confidence ellipses   ####################################################
###   COMPUTATIONS   ###########################################################
  if(is.list(conf)){
    # simple ellipse to replace the function ellipse() from the car package
    simple_ellipse <- function(center, a, b, n = 200L, border_col){
      angle_t <- seq(0, 2*pi, length.out = n)[-1L]
      polygon(center[1L] + a * cos(angle_t), center[2L] + b * sin(angle_t), lwd = 1.0, border = border_col)
    }

    # select items for which ellipses are drawn  ## rh 2011-05-31
    if(is.null(conf$which)) conf$which <- beta.subset#seq_along(beta.subset)
    if(!all(conf$which %in% beta.subset)) stop("Incorrect item number(s) for which ellipses are to be drawn")
    if(is.null(conf$col)){
      conf$c <- rep("red", length.out = length(beta1))
    } else if(!is.null(conf$which)){
      # conf$c <- rep(NA,length.out=length(beta.subset))
      conf$c <- rep(NA, length.out = length(conf$which))
      if(length(conf$c)!=length(conf$which)){
        stop('"which" and "col" must have the same length in specification of "conf"')
      } else {
        conf$c[conf$which] <- conf$col
      }
    }
    conf$col <- conf$c

    if(is.null(conf$gamma)) conf$gamma <- 0.95
    if(is.null(conf$lty)) conf$lty <- "dotted"
    if(is.null(conf$ia)) conf$ia <- FALSE

    z <- qnorm((1.0-conf$gamma)/2.0, lower.tail = FALSE)

    ci1u <- beta1 + z*s1
    ci1l <- beta1 - z*s1
    ci2u <- beta2 + z*s2
    ci2l <- beta2 - z*s2
  }
################################################################################

###   95% control lines (Wright)   #############################################
###   COMPUTATIONS   ###########################################################
  if(is.list(ctrline)){
    if(is.null(ctrline$gamma)) ctrline$gamma <- 0.95
    if(is.null(ctrline$col))   ctrline$col <- "blue"
    if(is.null(ctrline$lty))   ctrline$lty <- "solid"

    z <- qnorm((1.0 - ctrline$gamma)/2.0, lower.tail = FALSE)

    d      <- (beta1 + beta2)/2
    se.d   <- sqrt(v1 + v2)
    d      <- sort(d)
    se.d   <- se.d[order(d)]
    upperx <- d - z*se.d/2
    uppery <- d + z*se.d/2
  }
################################################################################

  if(!exists("ci1l", inherits = FALSE)) ci1l <- NA
  if(!exists("ci1u", inherits = FALSE)) ci1u <- NA
  if(!exists("ci2l", inherits = FALSE)) ci2l <- NA
  if(!exists("ci2u", inherits = FALSE)) ci2u <- NA

  if(!exists("upperx", inherits = FALSE)) upperx <- NA
  if(!exists("uppery", inherits = FALSE)) uppery <- NA

  if(missing(xlim)){
    xlim <- range(beta1[beta.subset], ci1l, ci1u, upperx, uppery, na.rm = TRUE)
  }
  if(missing(ylim)){
    ylim <- range(beta2[beta.subset], ci2l, ci2u, upperx, uppery, na.rm = TRUE)
  }

  plot.new()
  plot.window(xlim = xlim, ylim = ylim, asp = asp)
  title(main = main, xlab = xlab, ylab = ylab)
  if(x_axis) axis(1)
  if(y_axis) axis(2)

  abline(0, 1)

# confidence ellipses - if not interactive
  if(is.list(conf) && !conf$ia){
    # non-interactive: plot of all ellipses at once
    x <- beta1
    y <- beta2
    for(i in beta.subset){
      if(i %in% conf$which){
        segments( x0 = c(x[i], ci1l[i]), y0 = c(ci2l[i], y[i]),
                  x1 = c(x[i], ci1u[i]), y1 = c(ci2u[i], y[i]),
                  col = conf$col[i], lty = conf$lty )
        simple_ellipse( center = c(x[i],y[i]),
                        a = abs(diff(c(ci1u[i],ci1l[i])))/2,
                        b = abs(diff(c(ci2u[i],ci2l[i])))/2,
                        n = 200L, border_col = conf$col[i] )
      }
    }
  }

# 95% control lines (Wright) - plotting
  if(is.list(ctrline)){
    lines(upperx, uppery, col = ctrline$col, lty = ctrline$lty)
    lines(uppery, upperx, col = ctrline$col, lty = ctrline$lty)
  }

  if(exists("textlab", inherits = FALSE)){
    text(beta1[beta.subset], beta2[beta.subset], labels = textlab, pos = pos, ...)
  }

  points(x = beta1[beta.subset], y = beta2[beta.subset], type = type, ...)

  if(exists("labs", inherits = FALSE)){
    xycoords <- cbind(beta1[beta.subset], beta2[beta.subset])
    nothing  <- identify(xycoords, labels = labs, atpen = TRUE, offset = 1)
  }

  box()



# interactive confidence ellipses
  if(is.list(conf) && conf$ia){
    identifyEll <- function(x, y, ci1u, ci1l, ci2u,ci2l, v1, v2, conf, n=length(x), ...){
    ## source: example from help("identify")
    ## a function to use identify to select points, and overplot the
    ## points with a confidence ellipse as they are selected
      xy <- xy.coords(x, y)
      x <- xy$x
      y <- xy$y
      sel <- rep(FALSE, length(x))
      res <- integer(0)
      while(sum(sel) < n){
        ans <- identify(x[!sel], y[!sel], n=1, plot=FALSE, ...)
        if(!length(ans)) break
        ans <- which(!sel)[ans]
        i <- ans
        segments( x0 = c(x[i], ci1l[i]), y0 = c(ci2l[i], y[i]),
                  x1 = c(x[i], ci1u[i]), y1 = c(ci2u[i], y[i]),
                  col = conf$col[i], lty = conf$lty )
        simple_ellipse( center = c(x[i],y[i]),
                        a = abs(diff(c(ci1u[i],ci1l[i])))/2,
                        b = abs(diff(c(ci2u[i],ci2l[i])))/2,
                        n = 200L, border_col = conf$col[i] )
        sel[ans] <- TRUE
        res <- c(res, ans)
      }
      #res
    }
    identifyEll(beta1[beta.subset],beta2[beta.subset],
                        ci1u[beta.subset], ci1l[beta.subset], ci2u[beta.subset], ci2l[beta.subset],
                        v1[beta.subset], v2[beta.subset], conf)
  }

}
