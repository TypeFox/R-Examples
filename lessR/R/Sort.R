Sort <-
function(by, direction=NULL, data=mydata, quiet=getOption("quiet"), ...) {


  if (missing(by)) {
    cat("\n"); stop(call.=FALSE, "\n","------\n",
      "Specify the variables to sort by listing them first or with:  by\n\n")
  }

  dname <- deparse(substitute(data))
  n.obs <- nrow(data)

  all.vars <- as.list(seq_along(data))
  names(all.vars) <- names(data)

  if (!quiet) {
    cat("\nSort Specification\n")
    .dash(31)
  }

  if (deparse(substitute(by)) == "row.names") {
    ord <- "order(row.names(data)"
    cat(" ", "row.names", "-->")
    if (!is.null(direction)) {
      if (direction[1] == "+") txt <- "ascending"
      else if (direction[1] == "-") txt <- "descending"
      else {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
        "Value of direction, the sort direction specification: ", direction[1], "\n\n",
        "Only permissible values are + for ascending and - for descending.\n\n")
      }
    }
    else 
      txt <- "ascending"
    ord.txt <- "decreasing=FALSE"
    if (txt =="descending") ord.txt <- "decreasing=TRUE"
    ord <- paste(ord, ",", ord.txt, ",...)", sep="")
    cat(" ", txt, "\n") 
  }

  else if (deparse(substitute(by)) == "random") {
    cat(" ", "random\n")
    rand.rows <- sample(1:n.obs, size=n.obs, replace=FALSE)
    ord <- paste("order(", "rand.rows", ", ...)", sep="")
  }

  else {  # sort variable(s)

    by.col <- eval(substitute(by), envir=all.vars, enclos=parent.frame())
    n.sort <- length(by.col)

    if (!is.null(direction)) {
      if (n.sort != length(direction)) { 
      cat("\n"); stop(call.=FALSE, "\n","------\n",
        "Number of specified variables to sort: ", n.sort, "\n",
        "Number of + and - signs to indicate direction of sort: ", 
          length(direction), "\n\n",
        "The same number of values must be specified for both\n",
        "the list of values and the list of the sort direction.\n\n")
      }
    }
    else for (i in 1:n.sort) direction[i] <- "+"

    # console output
    if (!quiet) {
      for (i in 1:n.sort) {
        cat(" ", names(data)[by.col[i]], "-->")
        if (direction[i] == "+") txt <- "ascending"
        else if (direction[i] == "-") txt <- "descending"
        else {
          cat("\n"); stop(call.=FALSE, "\n","------\n",
          "Value of direction, the sort direction specification: ", direction[i], "\n\n",
          "Only permissible values are + for ascending and - for descending.\n\n")
        }
        cat(" ", txt, "\n") 
      }
    }

    fvar <- matrix(nrow=n.obs, ncol=n.sort)  # max num of factors is n.sort
    i.fvar <- 0

    # construct the call to the order function
    ord <- ""
    for (i in 1:n.sort) { 
      if ("factor" %in% class(data[, by.col[i]])) {  # 2 attributes if ordered
        i.fvar <- i.fvar + 1  # another factor variable
        fvar[, i.fvar] <- as.character(data[, by.col[i]])
        for (i.row in 1:n.obs)  # replace value with factor integer prefixed 
          fvar[i.row, i.fvar] <- 
             paste(toString(as.numeric(data[, by.col[i]])[i.row]), 
                   as.character(fvar[i.row, i.fvar]), sep="")
        ord <- paste(ord, direction[i], "xtfrm(fvar[, ", toString(i.fvar), "]), ",
                     sep="")
      }
      else   # variable not a factor
        ord <- paste(ord, direction[i], 
                     "data[, by.col[", toString(i), "]], ", sep="")
    }
    ord <- paste("order(", ord, "...)", sep="")
  }


  # finish the console output
  if (!quiet) .dash(31)
  cat("\n")


  # do the sort
  o <- eval(parse(text=ord))
  mydata <- data[o, ]

  if (!quiet) {
    cat("\n")
    .dash(68)
    cat("After the Sort, first four rows of data ")
    cat( "\n")
    .dash(68)
    print(head(mydata, n=4))
    cat("\n")
  }

  return(mydata)

}
