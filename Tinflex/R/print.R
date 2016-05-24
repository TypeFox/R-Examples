#############################################################################
##
##  Print data about generator.
##
#############################################################################

print.Tinflex <- function(x, debug=FALSE, ...) {
  ## ------------------------------------------------------------------------
  ## S3 method for printing class 'Tinflex'.
  ## ------------------------------------------------------------------------
  ##   x     ... S3 object of class 'Tinflex'
  ##   debug ... if TRUE print intervals
  ## ------------------------------------------------------------------------

  ## Prepare log-density.
  lpdf.args <- paste(names(formals(x$lpdf)), collapse=",")
  lpdf.body <- gsub("\\s", "", paste(deparse(x$lpdf), collapse=""), perl=TRUE)
  lpdf.body <- gsub("^function\\(\\w+\\)", "", lpdf.body, perl=TRUE)

  ## Compute total area below squeeze.
  A.sq.tot <- sum(x$ivs["A.sq",], na.rm=TRUE)

  ## Print on console.
  cat("Object of class 'Tinflex':\n\n")
  cat("       log-density(",lpdf.args,") = ",lpdf.body,"\n\n", sep="")

  cat("       area below hat =",x$A.ht,"\n")
  cat("   area below squeeze =",A.sq.tot,"\n")
  cat("    ratio hat/squeeze =",x$A.ht/A.sq.tot,"   [= upper bound for rejection constant ]\n")
  cat("          # intervals =",ncol(x$ivs)-1,"\n\n")

  ## Print boundaries and 'c' values for initial intervals.
  cat(x$iniv,"\n")
  
  if (is.TRUE(debug)) {
    ## Print data about for intervals.
    cat("Interval data:\n")
    print(x$ivs)
    cat("\nCumulated areas:\n")
    print(x$Acum)
    cat("\nGuide table:\n")
    print(x$gt)
    cat("\n")
  }
}

## --------------------------------------------------------------------------

