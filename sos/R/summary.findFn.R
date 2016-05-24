summary.findFn <- function(object, minPackages = 12,
                           minCount=NA, ...) {
##
## 1.  Retrieve the summary attribute
##
  Sum <- attr(object, 'PackageSummary')
##
## 2.  Limit it
##
  nrows <- nrow(Sum)
  minPackages <- min(nrows, minPackages, na.rm=TRUE)
  minCount <- {
    if(minPackages<1) 0 else
    min(minCount, Sum$Count[minPackages], na.rm=TRUE)
  }
#
  sel <- (Sum[, 'Count'] >= minCount)
  sumTh <- Sum[sel,, drop=FALSE]
  structure(list(PackageSummary = sumTh,
                 minPackages = minPackages,
                 minCount = minCount,
                 matches = attr(object, "matches"),
                 nrow = nrow(object),
                 nPackages = length(sel),
                 string = attr(object, 'string'),
                 call = attr(object, "call")),
            class = c("summary.findFn", "list"))
}

print.summary.findFn <- function(x, ...) {
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Total number of matches: ", sum( x$matches) , "\n", sep = "")
  cat('Downloaded ', x$nrow, ' links in ', x$nPackages,
      " package", c('.', 's.')[1 + (x$nPackages > 1)], "\n\n", sep='')
  string <- x$string
  cat("Packages with at least ", x$minCount, " match",
      if(x$minCount == 1) "" else "es",
      " using pattern\n  '", string, "'\n", sep = "")
  packSum <- x$PackageSummary
  ## strip of %H:%M:%S time stamp which won't make any difference
  packSum$Date <-
    substr(packSum$Date, 1, regexpr(" ", packSum$Date) - 1)
  row.names(packSum) <- NULL
  print(packSum, ...)
  invisible()
}
