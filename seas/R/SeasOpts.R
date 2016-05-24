"setSeasOpts" <-
function() {
  tstdate <- as.Date("2001-08-01")  # try and left-justify the day of month
  monthday <- if (format(tstdate, "%-d") == "1") {
    "%b %-d"  # unix-like
  } else if (format(tstdate, "%#d") == "1") {
    "%b %#d"  # Windows
  } else {  # other; two spaces if the day of month is < 10)
    "%b %d"
  }
  options(seas.main=list(fmt="%s\n%s", rngsep="-", show.id=TRUE, show.fun=TRUE),
          seas.label=list(fmt="%s (%s)", monthday=monthday, month="%B",
            ann="annual"),
          seas.month.grid=list(abb=TRUE, len=NULL, force=TRUE, label=TRUE,
            col="lightgrey", lwd=1, lty=1),
          seas.bxp=list(boxcol="lightgrey", outcex=1),
          seas.temp=list(col=c("lightgrey", "red"), lwd=3),
          seas.precip=list(col="grey", density=NULL, angle=45, lwd=1),
          seas.rain=list(col="#66CCFF", density=NULL, angle=45, lwd=1),
          seas.snow=list(col="#666666", density=NULL, angle=-45, lwd=1),
          seas.interarrival=list(col=c("lightblue", "orange")),
          seas.median=list(col="red", lwd=1, lty=1),
          seas.mean=list(col="orange", lwd=1, lty=1),
          seas.na=list(col="red", pch="x")
          )
}
