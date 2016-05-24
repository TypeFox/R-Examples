#' @title Generic methods for objects with class \code{'wrsp'}.
#' 
#' @description
#' An object with class 'wrsp' is a list containing the
#' response information, and the
#' mechanical, hydraulic, and material properties used to
#' generate the response for a sealed well.
#' 
#' @details
#' The response information is a
#' matrix with frequency, complex response
#' [\eqn{\omega}, \eqn{Z_\alpha (\omega)}]
#' where the units of \eqn{\omega} will be as they were input.
#' The amplitude of \eqn{Z}
#' is in meters per strain, 
#' and the phase is in radians.
#' 
#' @name wrsp-methods
#' @aliases wrsp
#' @rdname wrsp-methods
#' @docType methods
#'
#' @seealso 
#' \code{\link{well_response}}
#' 
#' \code{\link{kitagawa-package}}
#' 
#' @author A. J. Barbour <andy.barbour@@gmail.com>
#' 
#' @param x 'wrsp' object
#' @param object 'wrsp' object
#' @param n numeric; the number of \code{head} and \code{tail} to print
#' @param series character; the series to plot (amplitude or phase)
#' @param pch point character, as in \code{\link{par}}
#' @param xlims limits for x-axis (applies to both amp and phs frames)
#' @param ylims optional list of limits for y-axis (i.e., \code{list(amp=c(..),phs=c(...))})
#' @param logamp logical; should the amplitude be in log10 units
#' @param ... optional arguments
#' 
#' @examples
#' W <- well_response(1:10, T.=1, S.=1, Vw.=1, Rs.=1, Ku.=1, B.=1)
#' str(W)
#' print(W)
#' print(summary(W))
#' #
#' # Plot the response
#' plot(rnorm(10), xlim=c(-1,11), ylim=c(-2,2))
#' lines(W)
#' lines(W, "phs", col="red")
#' points(W)
#' points(W, "phs")
#' #
#' Wdf <- as.data.frame(W)
#' plot(Mod(wellresp) ~ omega, Wdf) # amplitude
#' plot(Arg(wellresp) ~ omega, Wdf) # phase
#' #
#' # or use the builtin method plot.wrsp
#' plot(W)
#' # change limits:
#' plot(W, xlims=c(-1,1), ylims=list(amp=c(5,8), phs=185*c(-1,1)))
NULL

#' @rdname wrsp-methods
#' @aliases as.data.frame.wrsp
#' @method as.data.frame wrsp
#' @S3method as.data.frame wrsp
as.data.frame.wrsp <- function(x, ...){
  WR <- x[["Response"]]
  df <- as.data.frame(WR)
  #names(df) <- "n.wrsp"
  return(df)
}
#' @rdname wrsp-methods
#' @aliases data.frame.wrsp
#' @method data.frame wrsp
#' @S3method data.frame wrsp
data.frame.wrsp <- as.data.frame.wrsp

#' @rdname wrsp-methods
#' @aliases print.wrsp
#' @method print wrsp
#' @S3method print wrsp
print.wrsp <- function(x, n=3, ...){
  stopifnot(is.wrsp(x))
  message("Sealed well response:")
  WR <- x[["Response"]]
  WR <- as.data.frame(WR)
  print(head(WR, n))
  message("\t...")
  print(tail(WR, n))
}

#' @rdname wrsp-methods
#' @aliases summary.wrsp
#' @method summary wrsp
#' @S3method summary wrsp
summary.wrsp <- function(object, ...){
  stopifnot(is.wrsp(object))
  WR <- object[["Response"]]
  toret <- summary.default(WR)
  class(toret) <- "summary.wrsp"
  return(toret)
}

#' @rdname wrsp-methods
#' @aliases print.summary.wrsp
#' @method print summary.wrsp
#' @S3method print summary.wrsp
print.summary.wrsp <- function(x, ...){
  message("Sealed well response summary:")
  print.summaryDefault(x)
}

#' @rdname wrsp-methods
#' @aliases lines.wrsp
#' @method lines wrsp
#' @S3method lines wrsp
lines.wrsp <- function(x, series=c("amp","phs"), ...){
  stopifnot(is.wrsp(x))
  series <- match.arg(series)
  WR <- x[["Response"]]
  x <- WR[,1]
  y <- WR[,2]
  y <- switch(series, amp=Mod(y), phs=Arg(y))
  graphics::lines(x, y, ...)
}

#' @rdname wrsp-methods
#' @aliases points.wrsp
#' @method points wrsp
#' @S3method points wrsp
points.wrsp <- function(x, series=c("amp","phs"), pch="+", ...){
  stopifnot(is.wrsp(x))
  series <- match.arg(series)
  WR <- x[["Response"]]
  x <- WR[,1]
  y <- WR[,2]
  y <- switch(series, amp=Mod(y), phs=Arg(y))
  graphics::points(x, y, pch=pch, ...)
}

#' @rdname wrsp-methods
#' @aliases plot.wrsp
#' @method plot wrsp
#' @S3method plot wrsp
plot.wrsp <- function(x, 
                      xlims=c(-3,1), 
                      ylims=list(amp=NULL, phs=185*c(-1,1)), logamp=TRUE, ...){
  #
  stopifnot(is.wrsp(x))
  WR <- x[["Response"]]
  au <- x[["Response.units"]]
  fu <- x[["Omega"]][["Units"]]
  mdl <- x[["Model"]][["Model"]]
  # enforce freq units in Hz
  fc <- switch(fu, rad_per_sec=2*pi, Hz=1)
  stopifnot(!is.null(fc))
  frq <- log10(WR[,1] / fc)
  amp <- Mod(WR[,2])
  if (logamp){
    amp <- log10(amp)
    au <- paste("log10",au)
  }
  phs <- Arg(WR[,2])*180/pi
  ##
  origpar <- par(no.readonly = TRUE)
  par(mar=c(2,4,1,1), 
      oma=c(2,0.1,1,0.1), 
      tcl=-0.3,
      mgp=c(2.5, 0.5, 0), las=1)
  layout(matrix(c(1,2), ncol=1), heights=c(0.5,0.5))
  # amplitude
  alims <- ylims[["amp"]]
  if (is.null(alims)) alims <- range(pretty(amp))
  plot(0,0,col=NA,
       ylim=alims,
       yaxs="i", ylab=sprintf("[%s]",au),
       xlim=xlims,
       xaxs="i", xaxt="n", xlab=""
  )
  lines(frq, amp, type="l", lwd=1.5, ...)
  log10_ticks()
  mtext(sprintf("Sealed well-response (%s)",mdl), font=2, line=1.0, cex=1.0)
  mtext("(a) Amplitude", adj=0.015, font=4, line=0.1, cex=0.8)
  # phase shift
  plot(0,0,col=NA,
       ylim=ylims[["phs"]],
       yaxs="i", yaxt="n", ylab="[degrees]",
       xlim=xlims,
       xaxs="i", xaxt="n", xlab=""
  )
  abline(h=c(-1,-0.5,0,0.5,1)*180, col="grey80", lty=2)
  lines(frq, phs, type="l", lwd=1.5, ...)
  lbls <- ats <- seq(-180,180,by=30)
  lbls[seq_along(lbls)%%2==0] <- ""
  axis(2, at=ats, labels=lbls)
  log10_ticks()
  mtext("(b) Phase", adj=0.015, font=4, line=0.1, cex=0.8)
  mtext("Frequency [Hz]", side=1, line=2)
  on.exit(par(origpar))
  return(invisible(NULL))
}

#' @details \code{\link{kitplot}} was previously a standalone function, but
#' is now simply a reference to \code{\link{plot.wrsp}}.
#' @rdname wrsp-methods
#' @export
#' @family PlotUtilities
kitplot <- function(x, ...) UseMethod("kitplot")
#' @rdname wrsp-methods
#' @aliases kitplot.wrsp
#' @method kitplot wrsp
#' @S3method kitplot wrsp
kitplot.wrsp <- plot.wrsp

#' @title Generic methods for objects with class \code{'owrsp'}.
#' 
#' @description
#' An object with class 'owrsp' is a list containing the
#' response information, and the
#' mechanical, hydraulic, and material properties used to
#' generate the response for an open well.
#' 
#' @details
#' The response information is a
#' matrix with frequency, complex response
#' [\eqn{\omega}, \eqn{Z_\alpha (\omega)}]
#' where the units of \eqn{\omega} will be as they were input.
#' The amplitude of \eqn{Z}
#' is in meters per strain, 
#' and the phase is in radians.
#' 
#' @name owrsp-methods
#' @aliases owrsp
#' @rdname owrsp-methods
#' @docType methods
#'
#' @seealso 
#' \code{\link{open_well_response}}
#' 
#' \code{\link{kitagawa-package}}
#' 
#' @author A. J. Barbour <andy.barbour@@gmail.com>
#' 
#' @param x 'owrsp' object
#' @param object 'owrsp' object
#' @param n numeric; the number of \code{head} and \code{tail} to print
#' @param series character; the series to plot (amplitude or phase)
#' @param pch point character, as in \code{\link{par}}
#' @param xlims limits for x-axis (applies to both amp and phs frames)
#' @param ylims optional list of limits for y-axis (i.e., \code{list(amp=c(..),phs=c(...))})
#' @param logamp logical; should the amplitude be in log10 units
#' @param ... optional arguments
#' 
#' @examples
#' S. <- 1e-5  	# Storativity [nondimensional]
#' T. <- 1e-4		# Transmissivity [m**2 / s]
#' frq <- 1/(1:200)
#' # Defaults to the Rojstaczer formulation
#' W <- open_well_response(frq, T. = T., S. = S., Rs. = 0.12, freq.units="Hz")
#' # (warning message about missing 'z')
#' W <- open_well_response(frq, T. = T., S. = S., Rs. = 0.12, freq.units="Hz", z=1)
#' str(W)
#' print(W)
#' print(summary(W))
# #
# # Plot the response
#' plot(rnorm(10), xlim=c(-1,11), ylim=c(-2,2))
#' lines(W)
#' lines(W, "phs", col="red")
#' points(W)
#' points(W, "phs")
#' #
#' Wdf <- as.data.frame(W)
#' plot(Mod(wellresp) ~ omega, Wdf) # amplitude
#' plot(Arg(wellresp) ~ omega, Wdf) # phase
# #
# # or use the builtin method
#' plot(W)
#' # change limits:
#' plot(W, xlims=c(-4,0), ylims=list(amp=c(-7,-3), phs=185*c(-1,1)))
NULL

#' @rdname owrsp-methods
#' @aliases as.data.frame.owrsp
#' @method as.data.frame owrsp
#' @S3method as.data.frame owrsp
as.data.frame.owrsp <- function(x, ...){
  WR <- x[["Response"]]
  df <- as.data.frame(WR)
  #names(df) <- "n.owrsp"
  return(df)
}
#' @rdname owrsp-methods
#' @aliases data.frame.owrsp
#' @method data.frame owrsp
#' @S3method data.frame owrsp
data.frame.owrsp <- as.data.frame.owrsp

#' @rdname owrsp-methods
#' @aliases print.owrsp
#' @method print owrsp
#' @S3method print owrsp
print.owrsp <- function(x, n=3, ...){
  stopifnot(is.owrsp(x))
  message("Open well response:")
  WR <- x[["Response"]]
  WR <- as.data.frame(WR)
  print(head(WR, n))
  message("\t...")
  print(tail(WR, n))
}

#' @rdname owrsp-methods
#' @aliases summary.owrsp
#' @method summary owrsp
#' @S3method summary owrsp
summary.owrsp <- function(object, ...){
  stopifnot(is.owrsp(object))
  WR <- object[["Response"]]
  toret <- summary.default(WR)
  class(toret) <- "summary.owrsp"
  return(toret)
}

#' @rdname owrsp-methods
#' @aliases print.summary.owrsp
#' @method print summary.owrsp
#' @S3method print summary.owrsp
print.summary.owrsp <- function(x, ...){
  message("Open well response summary:")
  print.summaryDefault(x)
}

#' @rdname owrsp-methods
#' @aliases lines.owrsp
#' @method lines owrsp
#' @S3method lines owrsp
lines.owrsp <- function(x, series=c("amp","phs"), ...){
  stopifnot(is.owrsp(x))
  series <- match.arg(series)
  WR <- x[["Response"]]
  x <- WR[,1]
  y <- WR[,2]
  y <- switch(series, amp=Mod(y), phs=Arg(y))
  graphics::lines(x, y, ...)
}

#' @rdname owrsp-methods
#' @aliases points.owrsp
#' @method points owrsp
#' @S3method points owrsp
points.owrsp <- function(x, series=c("amp","phs"), pch="+", ...){
  stopifnot(is.owrsp(x))
  series <- match.arg(series)
  WR <- x[["Response"]]
  x <- WR[,1]
  y <- WR[,2]
  y <- switch(series, amp=Mod(y), phs=Arg(y))
  graphics::points(x, y, pch=pch, ...)
}

#' @rdname owrsp-methods
#' @aliases plot.owrsp
#' @method plot owrsp
#' @S3method plot owrsp
plot.owrsp <- function(x, 
                      xlims=c(-3,1), 
                      ylims=list(amp=NULL, phs=185*c(-1,1)), logamp=TRUE, ...){
  #
  stopifnot(is.owrsp(x))
  WR <- x[["Response"]]
  au <- x[["Response.units"]]
  fu <- x[["Omega"]][["Units"]]
  mdl <- x[["Model"]][["Model"]]
  # enforce freq units in Hz
  fc <- switch(fu, rad_per_sec=2*pi, Hz=1)
  stopifnot(!is.null(fc))
  frq <- log10(WR[,1] / fc)
  amp <- Mod(WR[,2])
  if (logamp){
    amp <- log10(amp)
    au <- paste("log10",au)
  }
  phs <- Arg(WR[,2])*180/pi
  ##
  origpar <- par(no.readonly = TRUE)
  par(mar=c(2,4,1,1), 
      oma=c(2,0.1,1,0.1), 
      tcl=-0.3,
      mgp=c(2.5, 0.5, 0), las=1)
  layout(matrix(c(1,2), ncol=1), heights=c(0.5,0.5))
  # amplitude
  alims <- ylims[["amp"]]
  if (is.null(alims)) alims <- range(pretty(amp))
  plot(0,0,col=NA,
       ylim=alims,
       yaxs="i", ylab=sprintf("[%s]",au), 
       xlim=xlims,
       xaxs="i", xaxt="n", xlab=""
  )
  lines(frq, amp, type="l", lwd=1.5, ...)
  log10_ticks()
  mtext(sprintf("Open well-response (%s)",mdl), font=2, line=1.0, cex=1.0)
  mtext("(a) Amplitude", adj=0.015, font=4, line=0.1, cex=0.8)
  # phase shift
  plot(0,0,col=NA,
       ylim=ylims[["phs"]],
       yaxs="i", yaxt="n", ylab="[degrees]",
       xlim=xlims,
       xaxs="i", xaxt="n", xlab=""
  )
  abline(h=c(-1,-0.5,0,0.5,1)*180, col="grey80", lty=2)
  lines(frq, phs, type="l", lwd=1.5, ...)
  lbls <- ats <- seq(-180,180,by=30)
  lbls[seq_along(lbls)%%2==0] <- ""
  axis(2, at=ats, labels=lbls)
  log10_ticks()
  mtext("(b) Phase", adj=0.015, font=4, line=0.1, cex=0.8)
  mtext("Frequency [Hz]", side=1, line=2)
  on.exit(par(origpar))
  return(invisible(NULL))
}