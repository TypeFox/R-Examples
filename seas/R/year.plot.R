"year.plot" <-
function(x, start.day = 1, precip.only = FALSE, precip.ylim,
         temp.ylim, na.cut = 10, ...) {
  orig <- as.character(substitute(x))[[1]]
  if (precip.only)
    var.precip <- "precip"
  else
    var.precip <- c("precip", "rain", "snow")
  sc.precip <- seas.df.check(x, orig, var.precip)
  var.precip <- sc.precip$var
  sc.temp <- seas.df.check(x, orig, "t_mean")
  main <- .seastitle(id=sc.precip$id, name=sc.precip$name,
                     orig=orig, range=sc.precip$year.range, ...)
  x$ann <- mkann(x, start.day)
  years <- levels(x$ann)
  n.years <- length(years)
  precip <- array(dim=c(length(years), 3),
                  dimnames=list(years, c("rain", "snow", "precip")))
  if (!precip.only) {
    precip[,2] <- tapply(x$rain, x$ann, sum, na.rm=TRUE)
    precip[,1] <- tapply(x$snow, x$ann, sum, na.rm=TRUE)
  }
  precip[,3] <- tapply(x$precip, x$ann, sum, na.rm=TRUE)
  ndays <- year.length(2000, attr(x$ann, "calendar"))
  sum.is.num <- function(d) return(sum(is.finite(d)))
  precip.na <- tapply(x$precip, x$ann, sum.is.num)
  temp.na <- tapply(x$t_mean, x$ann, sum.is.num)
  precip.na[is.na(precip.na)] <- 0
  temp.na[is.na(temp.na)] <- 0
  precip.na <- ndays - precip.na
  temp.na <- ndays - temp.na
  precip.na[precip.na < na.cut] <- NA
  temp.na[temp.na < na.cut] <- NA
  precip.na <- sqrt(precip.na / ndays)
  temp.na <- sqrt(temp.na / ndays)
  if (missing(precip.ylim))
    precip.ylim <- c(0, max(precip, na.rm=TRUE))
  else if (length(precip.ylim) == 1)
    precip.ylim <- c(0, precip.ylim)
  if (missing(temp.ylim))
    temp.ylim <- range(x$t_mean, na.rm=TRUE)
  if (precip.only)
    pl <- t(precip[,3])
  else
    pl <- t(precip[,1:2])
  if (is.null(sc.precip$units))
    sc.precip$units <- "mm"
  if (is.null(sc.precip$long.name))
    sc.precip$long.name <- "precipitation"
  sc.precip$units <- gettextf("%s/year", sc.precip$units)
  sc.precip$ylab <- .seasylab(orig, sc.precip$long.name, sc.precip$units)
  ylab1 <- sc.precip$ylab
  if (precip.only) {
    op <- getOption("seas.precip")
    col <- op$col
    density <- op$density
    angle <- op$angle
  } else {
    op.rain <- getOption("seas.rain")
    op.snow <- getOption("seas.snow")
    col <- c(op.snow$col, op.rain$col)
    density <- c(op.snow$density, op.rain$density)
    angle <- c(op.snow$angel, op.rain$angle)
  }
  old.par <- par(no.readonly=TRUE)
  on.exit(par(old.par))
  #mar=c(2.1, 4.1, 4.1, 2.1)
  #mai=c(0.3, 0.8, 0.3, 0.1)
  par(mfrow=c(2, 1), bty="n")
  par(mar=c(1.1, 4.1, 4.6, 0.6))
  if (inherits(start.day, c("POSIXct", "Date")))
    start.day <- as.integer(format(start.day, "%j"))
  year.lab <- levels(x$ann)
  if (start.day != 1) {
    year.lab <- sub("_", "\n", year.lab)
    mgp <- par("mgp")
    mgp[2] <- mgp[2] + 0.5
    par(mgp=mgp)
  }
  barplot(pl, space=0, names.arg=year.lab, ylab=ylab1, ylim=precip.ylim,
          col=col, density=density, angle=angle)
  title(main=main)
  axis(1, 1:n.years -0.5, labels=FALSE)
  na.h <- diff(precip.ylim) / 80
  op.na <- getOption("seas.na")
  rect(1:n.years - 0.5 - precip.na / 2, 0, 1:n.years - 0.5 + precip.na / 2,
       na.h, col=op.na$col, border=FALSE)
  na.temp <- is.na(tapply(x$t_mean, x$ann, mean))
  for (i in names(na.temp[na.temp])) {
    x <- rbind(x, NA)
    x$ann[nrow(x)] <- i
    x$t_mean[nrow(x)] <- 0  # put 0 in NULL bins to plot boxplots properly
  }
  degSymb <- iconv("\260", "latin1", "")
  if (is.null(sc.temp$units))
    sc.temp$units <- sprintf("%sC", degSymb)
  if (is.null(sc.temp$long.name))
    sc.temp$long.name <- "temperature"
  ylab2 <- .seasylab(long.name=sc.temp$long.name, units=sc.temp$units)
  par(mar=c(1.6, 4.1, 2.1, 0.6))
  boxplot(t_mean ~ ann, x, varwidth=TRUE, ylim=temp.ylim, ylab=ylab2,
          xaxt="n", bty="n", col=getOption("seas.temp")$col[1])
  axis(1, 1:n.years, FALSE)
  axis(3, 1:n.years, FALSE)
  abline(h=0)
  na.h <- diff(temp.ylim) / 80
  rect(1:n.years - temp.na / 2, -na.h, 1:n.years + temp.na / 2, na.h,
       col=op.na$col, border=FALSE)
}
