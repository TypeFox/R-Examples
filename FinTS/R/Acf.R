Acf <- function(x, lag.max = NULL,
         type = c("correlation", "covariance", "partial"),
         plot = TRUE, na.action = na.fail, demean = TRUE, ...){
  dots <- list(...)
#  
  ACF <- acf(x, lag.max = lag.max, type=type, 
         plot = FALSE, na.action = na.action, demean = demean)
  class(ACF) <- c("Acf", class(ACF))
#
  dots$x <- ACF
  if(!('main' %in% names(dots)))
    dots$main <- deparse(substitute(x))
#  
  do.call('plot', dots) 
#
  invisible(ACF)
}

plot.Acf <- function(x, ci = 0.95, type = "h", xlab = "Lag", ylab = NULL,
          ylim = NULL, main = NULL,
          ci.col = "blue", ci.type = c("white", "ma"),
          max.mfrow = 6, ask = Npgs > 1 && dev.interactive(),
          mar = if(nser > 2) c(3,2,2,0.8) else par("mar"),
          oma = if(nser > 2) c(1,1.2,1,1) else par("oma"),
          mgp = if(nser > 2) c(1.5,0.6,0) else par("mgp"),
          xpd = par("xpd"), cex.main = if(nser > 2) 1 else par("cex.main"),
          verbose = getOption("verbose"), acfLag0=FALSE, 
          ...){
    ci.type <- match.arg(ci.type)
    if((nser <- ncol(x$lag)) < 1)stop("x$lag must have at least 1 column")
    if (is.null(ylab))
        ylab <- switch(x$type,
                       correlation = "ACF",
                       covariance = "ACF (cov)",
                       partial = "Partial ACF")
    if (is.null(snames <- x$snames))
        snames <- paste("Series ", if (nser == 1) x$series else 1:nser)

    with.ci <- ci > 0 && x$type != "covariance"
    with.ci.ma <- with.ci && ci.type == "ma" && x$type == "correlation"
    if(with.ci.ma && x$lag[1,1,1] != 0) {
        warning("can use ci.type=\"ma\" only if first lag is 0")
        with.ci.ma <- FALSE
    }
    clim0 <- if (with.ci) qnorm((1 + ci)/2)/sqrt(x$n.used) else c(0, 0)

    Npgs <- 1 ## we will do [ Npgs x Npgs ] pages !
    nr <- nser
    if(nser > 1) { ## at most m x m (m := max.mfrow)  panels per page
        sn.abbr <- if(nser > 2) abbreviate(snames) else snames

        if(nser > max.mfrow) {
            ##  We need more than one page: The plots are laid out
            ##  such that we can manually paste the paper pages and get a
            ##  nice square layout with diagonal !
            ## NB: The same applies to pairs() where we'd want several pages
            Npgs <- ceiling(nser / max.mfrow)
            nr <- ceiling(nser / Npgs)  # <= max.mfrow
        }
        opar <- par(mfrow = rep(nr, 2), mar = mar, oma = oma, mgp = mgp,
                    ask = ask, xpd = xpd, cex.main = cex.main)
        on.exit(par(opar))
        if(verbose) { # FIXME: message() can be suppressed but not str()
            message("par(*) : ", appendLF=FALSE, domain = NA)
            str(par("mfrow","cex", "cex.main","cex.axis","cex.lab","cex.sub"))
        }
    }
#   Drop acfLag0?
    nlags <- dim(x$lag)[1]
    istrt <- 1
    if((x$type=="correlation") && !acfLag0) istrt <- 2
    indx <- (istrt:nlags)
#       
    if (is.null(ylim)) {
        ## Calculate a common scale
      ylim <- range(x$acf[indx, 1:nser, 1:nser], na.rm = TRUE)
      if (with.ci) ylim <- range(c(-clim0, clim0, ylim))
      if (with.ci.ma) {
        for (i in 1:nser) {
          clim <- clim0 * sqrt(cumsum(c(1, 2*x$acf[-1, i, i]^2)))
          ylim <- range(c(-clim, clim, ylim))
        }
      }
    }

    for (I in 1:Npgs) for (J in 1:Npgs) {
        ## Page [ I , J ] : Now do   nr x nr  'panels' on this page
      iind <- (I-1)*nr + 1:nr
      jind <- (J-1)*nr + 1:nr
      if(verbose)
        message("Page [",I,",",J,"]: i =",
                paste(iind,collapse=","),"; j =",
                paste(jind,collapse=","), domain = NA)
      for (i in iind) for (j in jind)
        if(max(i,j) > nser) {
          frame(); box(col = "light gray")
          ## the above is EXTREMELY UGLY; should have a version
          ## of frame() that really does advance a frame !!
        }
        else {
          clim <- if (with.ci.ma && i == j)
            clim0 * sqrt(cumsum(c(1, 2*x$acf[-1, i, j]^2))) else clim0
#           plot
          plot(x$lag[indx, i, j], x$acf[indx, i, j], type = type,
               xlab = xlab, ylab = if(j==1) ylab else "", ylim = ylim, ...)
          abline(h = 0)
          if (with.ci && ci.type == "white")
            abline(h = c(clim, -clim), col = ci.col, lty = 2)
          else if (with.ci.ma && i == j) {
            lines(x$lag[indx, i, j], clim, col = ci.col, lty = 2)
            lines(x$lag[indx, i, j], -clim, col = ci.col, lty = 2)
          }
          title(if (!is.null(main)) main else
                if (i == j) snames[i]
                else paste(sn.abbr[i], "&", sn.abbr[j]),
                line = if(nser > 2) 1 else 2)
        }
        if(Npgs > 1) {                  # label the page
          mtext(paste("[",I,",",J,"]"), side=1, line = -0.2, adj=1,
                col = "dark gray", cex = 1, outer = TRUE)
        }
    }
}
