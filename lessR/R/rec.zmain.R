.rec.main <-
function(x, x.name, new.var, old, new, ivar, n.obs, dname, quiet) {

  n.values <- length(old)

  miss.old <- FALSE
  if (!is.null(old)) if (old[1] == "missing") miss.old <- TRUE

  miss.new <- FALSE
  if (!is.null(new)) {
    if (new[1] != "missing") {
      if (n.values != length(new)) { 
      cat("\n"); stop(call.=FALSE, "\n","------\n",
        "The same number of values must be specified for both\n",
        "old and new values.\n\n")
      }
    }
    else {
      for (i in 1:n.values) new[i] <- "missing"
      miss.new <- TRUE
    }
  }

  # text output
  if (!quiet) {

    if (ivar == 1) {
      cat("\nRecoding Specification\n")
      .dash(22)
      for (i in 1:n.values) cat("  ", old[i], "-->", new[i], "\n")
      cat("\n")
      if (miss.new)
        cat("\nR represents missing data with a NA for 'not assigned'.\n\n")
      cat("Number of cases (rows) to recode:", n.obs, "\n")
      if (is.null(new.var))
        cat("\nReplace existing values of each specified variable",
            ", no value for option: new.var\n", sep="")
    }
    cat("\n")

    old.unique <- sort(unique(x))
    n.unique <- length(old.unique)
    cat("---  Recode:", x.name, "---------------------------------\n")

    if (class(x) == "numeric")
      cat("Unique values of", x.name, "in the data:", old.unique, "\n")
    else if (class(x) == "factor")
      cat("Unique values of", x.name, "in the data:", levels(x), "\n")
    cat("Number of unique values of", x.name, "in the data:", n.unique, "\n")

    # check to ensure that all values to recode exist in the data
    for (i in 1:n.values) {
      is.in <- FALSE
      for (j in 1:n.unique) 
        if (old[i] == old.unique[j]) is.in <- TRUE 
      if (!is.in) {
      cat(">>> Note: A value specified to recode, ", old[i], 
        ", is not in the data.\n\n", sep="")
      }
    }

    cat("Number of values of", x.name, "to recode:", n.values, "\n")

    if (!is.null(new.var)) cat("Recode to variable:", new.var, "\n")
  }  # end text 

  if (class(x) == "factor") x <- as.character(x)

  new.x <- x

  # the recode
  for (i in 1:n.values) {
    for (j in 1:n.obs) {
      if (!miss.old) {
        if (!is.na(new.x[j])) 
          if (x[j] == old[i]) new.x[j] <- ifelse (!miss.new, new[i], NA)
      } 
      else  # miss.old
        if (is.na(new.x[j])) new.x[j] <- new[1]
    } 
  }

  return(new.x)
 
}
