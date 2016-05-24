calendarHeat <- function(values, ncolors=99,color="r2b",date.form = "%Y-%m-%d") {
  varname=names(values)
  dates=time(values)
  Yr=unique(lubridate::year(dates))
  if (length(Yr)>=6) {
    start=grep(as.character(max(Yr)-5),dates)[1]
    values<-values[start:length(values),]}
  else {values<-values}

  dates=time(values)
  if (class(dates) == "character" | class(dates) == "factor" ) {
    dates <- strptime(dates, date.form)
  }
  caldat <- data.frame(value = values, dates = dates)

  min.date <- as.Date(paste(format(min(dates), "%Y"),
                            "-1-1",sep = ""))
  max.date <- as.Date(paste(format(max(dates), "%Y"),
                            "-12-31", sep = ""))
  dates.f <- data.frame(date.seq = seq(min.date, max.date, by="days"))

  caldat <- data.frame(date.seq = seq(min.date, max.date, by="days"), value = NA)
  dates <- as.Date(dates)
  caldat$value[match(dates, caldat$date.seq)] <- values

  caldat$dotw <- as.numeric(format(caldat$date.seq, "%w"))
  caldat$woty <- as.numeric(format(caldat$date.seq, "%U")) + 1
  caldat$yr <- as.factor(format(caldat$date.seq, "%Y"))
  caldat$month <- as.numeric(format(caldat$date.seq, "%m"))
  yrs <- as.character(unique(caldat$yr))
  d.loc <- as.numeric()
  for (m in min(yrs):max(yrs)) {
    d.subset <- which(caldat$yr == m)
    sub.seq <- seq(1,length(d.subset))
    d.loc <- c(d.loc, sub.seq)
  }
  caldat <- cbind(caldat, seq=d.loc)

  #color styles
  r2b <- c("#0571B0", "#92C5DE", "#F7F7F7", "#F4A582", "#CA0020") #red to blue
  r2g <- c("#D61818", "#FFAE63", "#FFFFBD", "#B5E384")   #red to green
  w2b <- c("#045A8D", "#2B8CBE", "#74A9CF", "#BDC9E1", "#F1EEF6")   #white to blue

  assign("col.sty", get(color))
  calendar.pal <- colorRampPalette((col.sty), space = "Lab")
  def.theme <- lattice.getOption("default.theme")
  cal.theme <-
    function() {
      theme <-
        list(
          strip.background = list(col = "transparent"),
          strip.border = list(col = "transparent"),
          axis.line = list(col="transparent"),
          par.strip.text=list(cex=0.8))
    }
  lattice.options(default.theme = cal.theme)
  yrs <- (unique(caldat$yr))
  nyr <- length(yrs)
  print(cal.plot <- levelplot(value~woty*dotw | yr, data=caldat,
                              as.table=TRUE,
                              aspect=.12,
                              layout = c(1, nyr%%7),
                              between = list(x=0, y=c(1,1)),
                              strip=TRUE,
                              main = paste("Calendar Heat Map of ", varname, sep = ""),
                              scales = list(
                                x = list(
                                  at= c(seq(2.9, 52, by=4.42)),
                                  labels = month.abb,
                                  alternating = c(1, rep(0, (nyr-1))),
                                  tck=0,
                                  cex = 0.7),
                                y=list(
                                  at = c(0, 1, 2, 3, 4, 5, 6),
                                  labels = c("Sunday", "Monday", "Tuesday", "Wednesday", "Thursday","Friday", "Saturday"),
                                  alternating = 1,
                                  cex = 0.6,
                                  tck=0)),
                              xlim =c(0.4, 54.6),
                              ylim=c(6.6,-0.6),
                              cuts= ncolors - 1,
                              col.regions = (calendar.pal(ncolors)),
                              xlab="" ,
                              ylab="",
                              colorkey= list(col = calendar.pal(ncolors), width = 0.6, height = 0.5),
                              subscripts=TRUE
  ) )
  panel.locs <- trellis.currentLayout()
  for (row in 1:nrow(panel.locs)) {
    for (column in 1:ncol(panel.locs))  {
      if (panel.locs[row, column] > 0)
      {
        trellis.focus("panel", row = row, column = column,
                      highlight = FALSE)
        xyetc <- trellis.panelArgs()
        subs <- caldat[xyetc$subscripts,]
        dates.fsubs <- caldat[caldat$yr == unique(subs$yr),]
        y.start <- dates.fsubs$dotw[1]
        y.end   <- dates.fsubs$dotw[nrow(dates.fsubs)]
        dates.len <- nrow(dates.fsubs)
        adj.start <- dates.fsubs$woty[1]

        for (k in 0:6) {
          if (k < y.start) {
            x.start <- adj.start + 0.5
          } else {
            x.start <- adj.start - 0.5
          }
          if (k > y.end) {
            x.finis <- dates.fsubs$woty[nrow(dates.fsubs)] - 0.5
          } else {
            x.finis <- dates.fsubs$woty[nrow(dates.fsubs)] + 0.5
          }
          grid::grid.lines(x = c(x.start, x.finis), y = c(k -0.5, k - 0.5),
                     default.units = "native", gp=grid::gpar(col = "grey", lwd = 1))
        }
        if (adj.start <  2) {
          grid::grid.lines(x = c( 0.5,  0.5), y = c(6.5, y.start-0.5),
                     default.units = "native", gp=grid::gpar(col = "grey", lwd = 1))
          grid::grid.lines(x = c(1.5, 1.5), y = c(6.5, -0.5), default.units = "native",
                     gp=grid::gpar(col = "grey", lwd = 1))
          grid::grid.lines(x = c(x.finis, x.finis),
                     y = c(dates.fsubs$dotw[dates.len] -0.5, -0.5), default.units = "native",
                     gp=grid::gpar(col = "grey", lwd = 1))
          if (dates.fsubs$dotw[dates.len] != 6) {
            grid::grid.lines(x = c(x.finis + 1, x.finis + 1),
                       y = c(dates.fsubs$dotw[dates.len] -0.5, -0.5), default.units = "native",
                       gp=grid::gpar(col = "grey", lwd = 1))
          }
          grid::grid.lines(x = c(x.finis, x.finis),
                     y = c(dates.fsubs$dotw[dates.len] -0.5, -0.5), default.units = "native",
                     gp=grid::gpar(col = "grey", lwd = 1))
        }
        for (n in 1:51) {
          grid::grid.lines(x = c(n + 1.5, n + 1.5),
                     y = c(-0.5, 6.5), default.units = "native", gp=grid::gpar(col = "grey", lwd = 1))
        }
        x.start <- adj.start - 0.5

        if (y.start > 0) {
          grid::grid.lines(x = c(x.start, x.start + 1),
                     y = c(y.start - 0.5, y.start -  0.5), default.units = "native",
                     gp=grid::gpar(col = "black", lwd = 1.75))
          grid::grid.lines(x = c(x.start + 1, x.start + 1),
                     y = c(y.start - 0.5 , -0.5), default.units = "native",
                     gp=grid::gpar(col = "black", lwd = 1.75))
          grid::grid.lines(x = c(x.start, x.start),
                     y = c(y.start - 0.5, 6.5), default.units = "native",
                     gp=grid::gpar(col = "black", lwd = 1.75))
          if (y.end < 6  ) {
            grid::grid.lines(x = c(x.start + 1, x.finis + 1),
                       y = c(-0.5, -0.5), default.units = "native",
                       gp=grid::gpar(col = "black", lwd = 1.75))
            grid::grid.lines(x = c(x.start, x.finis),
                       y = c(6.5, 6.5), default.units = "native",
                       gp=grid::gpar(col = "black", lwd = 1.75))
          } else {
            grid::grid.lines(x = c(x.start + 1, x.finis),
                       y = c(-0.5, -0.5), default.units = "native",
                       gp=grid::gpar(col = "black", lwd = 1.75))
            grid::grid.lines(x = c(x.start, x.finis),
                       y = c(6.5, 6.5), default.units = "native",
                       gp=grid::gpar(col = "black", lwd = 1.75))
          }
        } else {
          grid::grid.lines(x = c(x.start, x.start),
                     y = c( - 0.5, 6.5), default.units = "native",
                     gp=grid::gpar(col = "black", lwd = 1.75))
        }

        if (y.start == 0 ) {
          if (y.end < 6  ) {
            grid::grid.lines(x = c(x.start, x.finis + 1),
                       y = c(-0.5, -0.5), default.units = "native",
                       gp=grid::gpar(col = "black", lwd = 1.75))
            grid::grid.lines(x = c(x.start, x.finis),
                       y = c(6.5, 6.5), default.units = "native",
                       gp=grid::gpar(col = "black", lwd = 1.75))
          } else {
            grid::grid.lines(x = c(x.start + 1, x.finis),
                       y = c(-0.5, -0.5), default.units = "native",
                       gp=grid::gpar(col = "black", lwd = 1.75))
            grid::grid.lines(x = c(x.start, x.finis),
                       y = c(6.5, 6.5), default.units = "native",
                       gp=grid::gpar(col = "black", lwd = 1.75))
          }
        }
        for (j in 1:12)  {
          last.month <- max(dates.fsubs$seq[dates.fsubs$month == j])
          x.last.m <- dates.fsubs$woty[last.month] + 0.5
          y.last.m <- dates.fsubs$dotw[last.month] + 0.5
          grid::grid.lines(x = c(x.last.m, x.last.m), y = c(-0.5, y.last.m),
                     default.units = "native", gp=grid::gpar(col = "black", lwd = 1.75))
          if ((y.last.m) < 6) {
            grid::grid.lines(x = c(x.last.m, x.last.m - 1), y = c(y.last.m, y.last.m),
                       default.units = "native", gp=grid::gpar(col = "black", lwd = 1.75))
            grid::grid.lines(x = c(x.last.m - 1, x.last.m - 1), y = c(y.last.m, 6.5),
                       default.units = "native", gp=grid::gpar(col = "black", lwd = 1.75))
          } else {
            grid::grid.lines(x = c(x.last.m, x.last.m), y = c(- 0.5, 6.5),
                       default.units = "native", gp=grid::gpar(col = "black", lwd = 1.75))
          }
        }
      }
    }
    trellis.unfocus()
  }
  lattice.options(default.theme = def.theme)
}



cutAndStack <- function(x, number, overlap = 0.1, type = 'l',xlab = "Time", ylab = deparse(substitute(x))) {
  date.year=substring(rownames(x),1,4)
  date.month=substring(rownames(x),6,7)
  date.day=substring(rownames(x),9,10)
  myTime=as.matrix(paste(date.year,paste(date.month,date.day,sep=""),sep="."))
  time <- if (is.ts(x)) time(x) else seq_along(x)
  Time <- equal.count(as.numeric(time),number = number, overlap = overlap)
  xyplot(as.numeric(x) ~ time | Time,
         type = type, xlab = xlab, ylab = ylab,
         default.scales = list(x = list(relation = "free"),
                               y = list(relation = "free")))
}
