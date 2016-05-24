Merge <-
function(data1, data2, by=NULL, quiet=getOption("quiet"), ...) {

  if (missing(data1)) {
    cat("\n"); stop(call.=FALSE, "\n","------\n",
      "Specify first data frame (table) to merge with:  data1\n\n")
  }

  if (missing(data2)) {
    cat("\n"); stop(call.=FALSE, "\n","------\n",
      "Specify second data frame (table) to merge with:  data2\n\n")
  }

  dname1 <- deparse(substitute(data1))
  dname2 <- deparse(substitute(data2))

  if (!quiet) {
    cat("\n")
    .dash(17)
    cat("Before the merge\n")
    .dash(17)
    cat("\n")
    cat("First five rows of data for first data frame:", dname1, "\n")
    .dash(68)
    print(head(data1, n=5))
    cat("\n")
    cat("First five rows of data for second data frame:", dname2, "\n")
    .dash(68)
    print(head(data2, n=5))
    cat("\n")
  }


  # do the merge
  if (missing(by)) {
    if (!identical(names(data1), names(data2))) { 
      cat("\n"); stop(call.=FALSE, "\n","------\n",
        "To do a vertical merge, both data sets must have the same variables.\n\n")
    }
    mylabels <- attr(data1, which="variable.labels") # save variable labels
    myunits <- attr(data1, which="variable.units") # save variable units
    type <- "vertical"
    data <- rbind(data1, data2)
  }

  else {
    type <- "horizontal"
    mylabels1 <- attr(data1, which="variable.labels") # save variable labels
    mylabels2 <- attr(data2, which="variable.labels") # save variable labels
    myunits1 <- attr(data1, which="variable.units") # save variable units
    myunits2 <- attr(data2, which="variable.units") # save variable units
    data <- merge(data1, data2, by=by, ...)
    mylabels <- c(mylabels1, mylabels2)
    myunits <- c(myunits1, myunits2)
  }


  if (!quiet) {
    cat("\n")
    .dash(16+nchar(type))
    cat("After the", type, "merge\n")
    .dash(16+nchar(type))
    cat("\n")
    cat("First five rows of data ")
    cat( "\n")
    .dash(68)
    print(head(data, n=5))
  }

  # restore any variable labels, units
  if (!is.null(mylabels)) attr(data, which="variable.labels") <- mylabels
  if (!is.null(myunits)) attr(data, which="variable.units") <- myunits

  return(data)

}
