Subset <-
function(rows, columns, data=mydata, holdout=FALSE,
    random=0, quiet=getOption("quiet"), ...) {

  # save variable labels, units (NULL if no labels, units) 
  mylabels <- attr(data, which="variable.labels")
  myunits <- attr(data, which="variable.units")

  dname <- deparse(substitute(data))

  n.vars <- ncol(data)
  n.obs <- nrow(data)

  if (!quiet) {
    cat("\n")
    .dash(17)
    cat("Before the subset\n")
    .dash(17)
    cat("\n")
    cat("Number of variables in ", dname, ": ", n.vars, "\n", sep="")
    cat("Number of cases (rows) in ", dname, ": ", n.obs, "\n\n", sep="")
    cat("First five rows of data for data frame:", dname, "\n")
    .dash(68)
    print(head(data, n=5))
    cat("\n")
  }

  # set rows parameter by a random selection
  if (random > 0) {
    if (!missing(rows)) {
      cat("\n"); stop(call.=FALSE, "\n","------\n",
        "Cannot specify both rows, to retain or exclude, and random.\n\n")
    }
    if (random > 1)
      n.obs.new <- round(random,0)
    else
      n.obs.new <- round(random*n.obs,0)
    if (!quiet) {
      cat("\n")
      cat("Rows of data randomly extracted\n")
      .dash(42)
      if (random <= 1) cat("Proportion of randomly retained rows: ", random, "\n")
      cat("Number of randomly retained rows: ", n.obs.new, "\n")
    }
    rand.rows <- sample(1:n.obs, size=n.obs.new, replace=FALSE)
    rand.rows <- sort(rand.rows)
    rows <- logical(length=n.obs)  # initial default is FALSE
    j <- 1
    for (i in 1:n.obs) {
      if (i == rand.rows[j]) {
        rows[i] <- TRUE 
        if (j < length(rand.rows)) j <- j + 1
      }
    }
  }

  
  # r is the logical vector or row to retain
  # set r by default to all or by specification
  if (missing(rows))
    r <- TRUE  # retain all rows
  else {  # if numeric, read directly from rows
    r <- eval(substitute(rows), envir=data, enclos=parent.frame())

  if (all(!r)) {
    cat("\n"); stop(call.=FALSE, "\n","------\n",
      "No rows of data would be retained.\n",
      "The specified value does not occur in these data.\n\n")
  }

    if (is.logical(r))
      r <- r & !is.na(r)  # set missing for a row to FALSE
  }

  if (is.numeric(r)) {
    if (length(rows) == 1) {
      if ((rows > 0)  &&  (rows < 1)) {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
          "Behavior has changed from previous documentation.\n",
          "Now use this parameter to draw a random subset:  random\n\n")
      }
      else if (rows > 1) {
        cat("\n"); warning(call.=FALSE, "\n","------\n",
          "Behavior has changed from previous documentation.\n",
          "Now use this parameter to draw a random subset:  random\n\n",
          "Now specifying a single value for rows yields only that\n",
          "  row in the subset, which is OK if that is what is desired.\n",
          "  A negative value deletes the row.\n\n")
      }
    }
  }

  # evaluate columns
  if (missing(columns)) 
    vars <- TRUE  # retain all variables
  else {
    all.vars <- as.list(seq_along(data))
    names(all.vars) <- names(data)
    vars <- eval(substitute(columns), envir=all.vars, parent.frame())
 }

  # do the subset, including any variable labels to be dropped
  data.sub <- data[r, vars, drop=FALSE]
  # if no vars dropped, vars is TRUE, otherwise a char vector of names to retain
  if (!is.logical(vars) && !is.null(mylabels)) {  
    # warning unless length(mydata) a multiple of length(vars)
    keep.index <- which(suppressWarnings(names(mydata) == vars))
    mylabels <- mylabels[keep.index]
  }

  if (!quiet) {
    cat("\n\n")
    .dash(16)
    cat("After the subset\n")
    .dash(16)
    cat("\n")
    cat("Number of variables ", ": ", ncol(data.sub), "\n", sep="")
    cat("Number of cases (rows) ", ": ", nrow(data.sub), "\n\n", sep="")
    cat( "\n")
    cat("First five rows of data ")
    cat( "\n")
    .dash(68)
    print(head(data.sub, n=5))
    cat( "\n")
  }

  if (holdout)  {
    data.hold <- data[!r, vars, drop=FALSE]
    cat("\n\n")
    cat("Holdout Sample\n")
    .dash(70)
    cat("Deleted Rows:", nrow(data.hold), "\n")
    cat("Copy and paste this code into R to create the holdout sample\n",
        "from the original, unmodified data frame:", dname, "\n\n")
    cat(paste(dname,".hold", sep=""), "<- Subset(\n")
    for (i in 1:nrow(data.hold)) {
      if (i > 1) cat ("| ") else cat("  ")
      cat("row.names(", dname, ")==\"", row.names(data.hold)[i], "\"\n", sep="")
    }
    cat(")\n")
  }
 
  # restore any variable units, labels
  # note:  variable labels for all variables of original data frame
  #        even if variables were deleted
  if (!is.null(mylabels)) attr(data.sub, which="variable.labels") <- mylabels
  if (!is.null(myunits)) attr(data.sub, which="variable.units") <- myunits

  return(data.sub)

}
