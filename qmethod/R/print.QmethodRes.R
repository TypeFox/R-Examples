print.QmethodRes <- function(x, length=10, digits=2, ...) {
  old.dig <- getOption("digits")
  options(digits=digits)
  nn <- c("Summary", "Original data", "Q-sort factor loadings", "Flagged Q-sorts", "Statement z-scores", "Statement factor scores", "Factor characteristics", "Distinguishing and consensus statements")
  names(nn) <- c("brief", "dataset", "loa", "flagged", "zsc", "zsc_n", "f_char", "qdc")
  ll <- length(x)
  nl <- nn[1:ll]
  dimsorts <- min(length, x$brief$nqsorts)
  dimstats <- min(length, x$brief$nstat)
  cat(x$brief$info, sep="\n")
  cat("\n")
  cat(nl[2], ":\n")
  print(x$dataset[1:dimstats, 1:dimsorts])
  if (dimstats < x$brief$nstat) cat(" (...) See item '...$dataset' for the full data.\n")
  nxt <- c("loa", "flagged")
  for (i in nxt) {
    cat("\n")
    cat(nl[i], ":\n")
    print(format(x[[i]][1:dimsorts, ], scientific=F))
    if (dimsorts < x$brief$nqsorts) cat(" (...) See item '...$", i, "' for the full data.\n", sep="")
  }
  nxt <- c("zsc", "zsc_n")
  for (i in nxt) {
    cat("\n")
    cat(nl[i], ":\n")
    print(format(x[[i]][1:dimstats, ], scientific=F))
    if (dimstats < x$brief$nstat) cat(" (...) See item '...$", i, "' for the full data.\n", sep="")
  }
  cat("\n", nl[7], ":\n", sep="")
  fcl <- c("   General factor characteristics:", "   Correlation between factor z-scores:", "   Standard error of differences between factors:")
  for (i in 1:length(x$f_char)) {
    cat(fcl[[i]], "\n")
    print(x$f_char[[i]])
    cat("\n")
  }
  if (ll == 8) {
    cat(nl[8], ":\n")
    print(x$qdc[1:dimstats, ])
    if (dimstats < x$brief$nstat)  cat(" (...) See item '...$qdc' for the full data.\n")
  }
  options(digits=old.dig)
  invisible(x)
}