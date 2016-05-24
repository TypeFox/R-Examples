"plot.seas.norm" <-
function(x, start = 1, rep = 0, ylim, varwidth = FALSE, normwidth = FALSE,
         leg, add.alt = FALSE, main, ylab, ...) {
  orig <- as.character(substitute(x))[[1]]
  if (!inherits(x, "seas.norm"))
    stop(gettextf("%s is not a %s object", sQuote(orig), sQuote("seas.norm")))
  fun <- x$fun
  precip.norm <- x$precip.norm
  width <- x$width
  seas <- x$seas
  seas$bin <- rownames(seas)
  ann <- x$ann
  var <- x$var
  if (missing(main))
    main <- .seastitle(id=x$id, orig=orig, name=x$name, fun=fun,
                       range=x$year.range, ...)
  xlab <- .seasxlab(width, x$start.day)
  n.bins <- length(x$bins)
  num <- n.bins + rep
  sel <- 1:num
  if (rep != 0 || start != 1) {
    sel <- (start:(start-1+num) - 1) %% n.bins + 1
    seas <- seas[sel,]
  }
  if (x$a.cut) {
    active <- seas$active
    active[is.na(active)] <- 0
    active[active == 0] <- 0.2  # to show at least a sliver of data
    if (varwidth) {  # normalize width of bars
      maxf <- max(seas$active, na.rm=TRUE)
      if (normwidth) {
        if (is.logical(normwidth))
          active <- active / max(active, na.rm=TRUE)
        else {
          med.days <- x$bin.lengths
          active <- active * med.days / normwidth  # assumed numeric
          if (normwidth < maxf)
            warning(
              paste(
                gettextf("%s < maximum width value; not a good plot",
                         sQuote("normwidth")),
                gettextf("this should be no less than %.1f", round(maxf, 1)),
                sep="\n"))
        }
      }
    } else active <- 1
  } else {
    active <- 1
  }
  unit <- x$units
  if (!is.null(unit) && unit != "") {
    alt.unit <- ifelse(unit=="mm", "in", ifelse(unit=="in", "mm", ""))
    if (is.null(fun) || fun %in% c("mean", "sd", "median", "min", "max", "var")) {
      dly.unit <- gettextf("%s/day", unit)
      alt.unit <- gettextf("%s/day", alt.unit)
    } else alt.unit <- ylab <- ""
    if (!is.null(fun) && fun == "var") {
      squareSym <- iconv("\262", "latin1", "")
      dly.unit <- sprintf("(%s)%s", unit, squareSym)
      alt.unit <- sprintf("(%s)%s", alt.unit, squareSym)
      aly.ylab <- alt.unit
    }
  } else unit <- dly.unit <- NULL
  if (missing(ylab))
    ylab <- if (precip.norm) dly.unit
      else .seasylab(orig, x$long.name, dly.unit)
  if (missing(ylim)) {
    maxy <- if (precip.norm)
      max(rowSums(seas[,c("rain", "snow")]), na.rm=TRUE)
      else
        max(seas[,var], na.rm=TRUE)
    ylim <- c(0, maxy) * 1.04
  } else if (length(ylim) == 1) {
    ylim <- c(0, ylim)
  }
  mar <- par("mar")
  par(yaxs="i", xaxs="r")  # keep it this way
  if (add.alt) {
    mar[4] <- mar[2]
    bty <- "u"
  } else {
    bty <- "l"
  }
  par(mar=mar, bty=bty)
  frame()
  plot.window(c(0.5, num + 0.5), ylim=ylim)
  .seasmonthgrid(width, x$bin.lengths, start, rep, x$start.day)
   lx <- 1:num - active / 2
   rx <- 1:num + active / 2
  bot <- 1:num * 0
  border <- if (varwidth) FALSE else TRUE
  lwd <- par("lwd")
  if (precip.norm) {
    # snow boxes
    op <- getOption("seas.snow")
    rect(lx, bot, rx, seas$snow, border=border,
         col=op$col, density=op$density, angle=op$angle, lwd=op$lwd * lwd)
    # rain boxes
    op <- getOption("seas.rain")
    rect(lx, seas$snow, rx, seas$snow + seas$rain, border=border,
         col=op$col, density=op$density, angle=op$angle, lwd=op$lwd * lwd)
  } else {  # precipitation only
    op <- getOption("seas.precip")
    rect(lx, bot, rx, seas[,var], border=border,
         col=op$col, density=op$density, angle=op$angle, lwd=op$lwd * lwd)
  }
  if (missing(leg))
    leg <- ifelse(is.null(fun) || fun %in% c("mean", "median"), TRUE, FALSE)
  if (is.logical(leg) && leg == TRUE) {
    leg.x <- 0.5 + num * 0.02
    leg.y <- max(ylim) * 0.98
  } else if (substitute(leg) == "locator") {
    xy <- locator(1)
    leg.x <- xy$x
    leg.y <- xy$y
    leg <- TRUE
  } else if (all(is.finite(leg)) && length(leg) == 2) {
    leg.x <- leg[1]
    leg.y <- leg[2]
    leg <- TRUE
  } else {
    leg <- FALSE
  }
  if (leg) {
    annrate <- paste(unit, gettext("year"), sep="/")
    if (!precip.norm)
      leg.text <- paste(gettext("Total"), var, round(ann[1, var], 1), annrate)
    else
      leg.text <- c(gettextf("Rain %.1f", ann$rain),
                    gettextf("Snow %.1f", ann$snow),
                    sprintf("%s %.1f %s",
                            gettext("Total"), ann$precip, annrate))
    if (x$a.cut)  # attach number of active days before rest of legend
      leg.text <- c(gettextf("%s days with %s", round(ann$active, 1),
                             ifelse(precip.norm, gettext("precipitation"),
                                    var)), leg.text)
    text(leg.x, leg.y, paste(leg.text, collapse="\n"), adj=c(0, 1))
  }
  na <- seas$na
  na.h = max(ylim) / 150
  na[na < 0.05] <- NA  # don't plot red box if there is < 5% data missing
  rect(1:num - na/2, 0, 1:num + na / 2, na.h,
       col=getOption("seas.na")$col, border=FALSE)
  box()
  abline(h=0)
  axis(1, 1:num, seas$bin, lwd=lwd)
  axis(2, lwd=lwd)
  title(main=main, xlab=xlab, ylab=ylab)
  if (add.alt) {
    mm2in <- function(v)(v / 25.4)
    in2mm <- function(v)(v * 25.4)
    if (unit == "mm") {
      alt.ax <- pretty(mm2in(ylim))
      axis(side=4, at=in2mm(alt.ax), labels=alt.ax, srt=90, lwd=lwd)
      mtext(alt.unit, side=4, line=2.8)
    } else if (unit=="in") {
      alt.ax <- pretty(in2mm(ylim))
      axis(side=4, at=mm2in(alt.ax), labels=alt.ax, srt=90, lwd=lwd)
      mtext(alt.unit, side=4, line=2.8)
    }
  }
}
