Transform <-
function(data=mydata, quiet=getOption("quiet"), ...) {

  # save variable labels (NULL if no labels) 
  mylabels <- attr(data, which="variable.labels")

  # transformations done here into vector e, e[1] the first trans, etc.
  e <- eval(substitute(list(...)), data, parent.frame())

  n.tvars <- length(e)
  n.obs <- nrow(data)

  if (!quiet) {
    dname <- deparse(substitute(data))

    cat("\n")
    cat("Number of variables of", dname, "to transform:", n.tvars, "\n")
    cat("Number of cases (rows) of ", dname, ": ", n.obs, "\n", sep="")

    cat("\n")
    .dash(69)
    cat("Before Transformation, First four rows of data for",
        "data frame:", dname, "\n")
    .dash(69)
    print(head(get(dname, pos=.GlobalEnv), n=4))
    cat("\n")

    trs.all <- deparse(substitute(list(...)))[[1]]
    trs.all <- substr(trs.all, 6, nchar(trs.all)-1)
    no.break <- FALSE
    if (grepl("factor", trs.all)) no.break <- TRUE
    if (grepl("cut", trs.all)) no.break <- TRUE
    if (!no.break) {
      trs <- as.vector(strsplit(trs.all, ",")[[1]])
      if (n.tvars > 1) for (i in 2:n.tvars) trs[i] <- substr(trs[i], 2, nchar(trs[i]))
    }
    else
      trs <- trs.all  # displayed later if only one transformation
  }

  # vector of positions of existing vars in data, NA if not existing
  inx <- match(names(e), names(data))

  if (!quiet) {
    cat("\n")
    .dash(22)
    cat("Transformation Summary\n")
    .dash(22)
    #if (!no.break || n.tvars==1) {
      #for (i in 1:n.tvars) {
        #if (!is.na(inx[i])) txt <- "rewrite existing" else txt <- "      create new"
        #cat(txt, "variable: ", trs[i], "\n")
      #}
    #}
    #else {  # parsing too difficult, so just include variable names
       for (i in 1:n.tvars) {
        if (!is.na(inx[i])) txt <- "rewrite existing" else txt <- "      create new"
        cat(txt, "variable: ", names(e)[i], "\n")
      }
    #}
  }

  # logical vector, TRUE if transformed var already exists, otherwise FALSE
  existing <- !is.na(inx)

  # inx[existing], e[existing] gives only positions of existing vars in data
  # corresponding FALSE values of existing are deleted
  # add existing var transformations to data
  if (any(existing))   # at least one existing var transformed
    data[inx[existing]] <- e[existing]

  # add new var transformations to data
  if (!all(existing))  # at least one new var created
    data <- do.call("data.frame", c(list(data), e[!existing]))

  if (!quiet) {
    cat("\n\n")
    .dash(65)
    cat("After, First four rows of transformed data ")
    cat( "\n")
    .dash(65)
    print(head(data[, names(e), drop=FALSE], n=4))
    cat("\n")
  }
 
  # restore any variable labels
  if (!is.null(mylabels)) attr(data, which="variable.labels") <- mylabels

  return(data)

}
