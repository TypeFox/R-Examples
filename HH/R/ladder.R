"ladder" <-
function(formula.in, data=sys.parent(),
         main.in="Ladders of Powers",
         panel.in=panel.cartesian,
         xlab=deparse(formula.in[[3]]),
         ylab=deparse(formula.in[[2]]),
         scales=list(alternating=if.R(s=TRUE, r=FALSE),
           labels=FALSE, ticks=FALSE, cex=.6),
         par.strip.text=list(cex=.6),
         cex=.5, pch=16, between=list(x=.3, y=.3),
         dsx=xlab,
         dsy=ylab,
         ladder.function=ladder.f,
         strip.number=if.R(r=2, s=1),
         strip.names,
         strip.style=1,
         strip,
         oma=c(0,0,0,0),  ## S-Plus
         axis3.line=.61,
         layout=c(length(tmp$x.power), length(tmp$y.power)),
         axis.key.padding = 10, ## R right axis
         key.axis.padding = 10, ## R top axis
         useOuter=TRUE, ## R useOuterStrips(combineLimits(result))
         ...) {
  ## plot all possible powers of y variable against x variable
  if ((length(formula.in[[2]]) != 1) || (length(formula.in[[3]]) != 1))
    stop("The left-hand side and right-hand side of the formula must each have exactly one variable.")
  tmp <- ladder3(data[xlab], data[ylab], dsx, dsy, ladder.function)
  if.R(r=if (useOuter) strip.number <- 2,
       s={})
  formula.out <- switch(as.character(strip.number),
                        "1"=formula(Y ~ X | group),
                        "0"=formula(Y ~ X | x * y),
                        "2"=formula(Y ~ X | x * y),
                        stop("strip.number must be 0, 1, or 2."))
  if (missing(strip.names)) strip.names <-
    switch(as.character(strip.number),
           "1"=c(FALSE,FALSE),
           "0"=c(TRUE,TRUE),
           "2"=c(TRUE,TRUE),
           stop("strip.number must be 0, 1, or 2."))
  if (missing(strip)) {
    strip.function <- strip.ladder
    if.R(s={
      strip.function$strip.names <- strip.names
      strip.function$style <- strip.style
    },r={
      sf.args <- formals(strip.function)
      sf.args$strip.names <- strip.names
      sf.args$style <- strip.style
      formals(strip.function) <- sf.args
    })} else strip.function <- strip
  if (strip.number==0) strip.function <- FALSE
  if.R(r={}, s=old.par <- par(oma=oma))
  if.R(r=scales$relation <- "free",
       s={})
  xyplot.args <- list(formula.out, data=tmp$data,
                      panel=panel.in,
                      xlab=xlab, ylab=ylab,
                      x.label=tmp$x, y.label=tmp$y,
                      strip=strip.function,
                      scales=scales,
                      par.strip.text=par.strip.text,
                      between=between,
                      cex=cex, pch=pch,
                      layout=layout,
                      main=main.in,
                      axis3.line=axis3.line,
                      ...)
  xyplot.R.args <- list(par.settings =
                        list(layout.widths=
                             list(axis.key.padding=axis.key.padding),
                             layout.heights=
                             list(key.axis.padding=key.axis.padding)))
  if.R(r=xyplot.args <- c(xyplot.args, xyplot.R.args),
       s={})
  result <- do.call("xyplot", xyplot.args)
  if (useOuter)
    if.R(r=result <- useOuterStrips(combineLimits(result)),
         s={})
  result
}

"ladder3" <-
function(x, y,
         dsx=deparse(substitute(x)),
         dsy=deparse(substitute(y)),
         ladder.function=ladder.f) {
  ## construct two ladder.f data.frames and construct all possible pairs
  if (length(x) != length(y)) stop("x and y must have the same length.")
  lfx <- ladder.function(x)
  n <- nrow(lfx)
  ncx <- ncol(lfx)
  xxx <- rep(as.list(lfx), rep(ncx,ncx))
  names(xxx) <- rep(names(lfx), rep(ncx,ncx))
  lfy <- ladder.function(y)
  yyy <- rep(as.list(lfy), ncx)
  names(yyy) <- rep(names(lfy), ncx)
  gx <- paste(dsx, "^", names(xxx), sep="")
  gx <- ordered(gx, levels=unique(gx))
  gy <- paste(dsy, "^", names(yyy), sep="")
  gy <- ordered(gy, levels=unique(gy))
  ggg <- paste(dsy, "^", names(yyy),
               " ~ ",
               dsx, "^", names(xxx),
               sep="")

  ggo <- ggg
  dim(ggo) <- c(length(names(lfy)), length(names(lfx)))
  ggo <- as.vector(t(ggo))
  ggg <- ordered(ggg, levels=ggo)

  result <- data.frame(unlist(xxx), unlist(yyy),
                       rep(gx, rep(n, ncx*ncx)),
                       rep(gy, rep(n, ncx*ncx)),
                       rep(ggg, rep(n, ncx*ncx)))
  names(result) <- c("X", "Y", "x", "y", "group")

  list(data=result, x.power=names(lfx), y.power=names(lfy))
}

ladder.f <- function(x, name.prefix="") {
  ## construct a data.frame, one column per power
  if (any(x <= 0)) warning('Non-positive values in argument to ladder.f.  Consider using a start value as "ladder.f(x+(min(x)+.5))".')
  result <- data.frame(-1/x, -1/sqrt(x), log(x), sqrt(x), x, x^2)
  names(result) <- paste(name.prefix, c(-1, -0.5, 0, 0.5, 1, 2), sep="")
  result
}

ladder.fstar <- function(x, name.prefix="") {
  ## construct a data.frame, one column per power.  Use the scaled Box--Cox formulas
  if (any(x <= 0)) warning('Non-positive values in argument to ladder.fstar.  Consider using a start value as "ladder.f(x+(min(x)+.5))".')
  result <- data.frame((1/x - 1)/(-1), (1/sqrt(x)-1)/(-.5),
                        log(x), (sqrt(x)-1)/.5, x-1, (x^2 - 1)/2)
  names(result) <- paste(name.prefix, c(-1, -0.5, 0, 0.5, 1, 2), sep="")
  result
}

"strip.ladder" <-
function(which.given,
         which.panel,
         var.name,
         factor.levels,
         shingle.intervals,
         par.strip.text=trellis.par.get("add.text"),
         strip.names=c(TRUE,TRUE),
         style=1,
         ...) {
  strip.default(which.given=which.given,
                which.panel=which.panel,
                var.name=var.name,
                factor.levels=factor.levels,
                shingle.intervals=shingle.intervals,
                par.strip.text=par.strip.text,
                strip.names=strip.names,
                style=style,
                ...)
}
