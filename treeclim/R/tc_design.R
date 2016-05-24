##' Deparse list structure from month specification and return parameter set
##'
##' Deparse list structure from month specification into single calls, call an
##' aggregating function to collect the data, and return a parameter set for
##' calibration. In general, the months must be specified as a list or nested
##' list structure. Possible cases are:
##' list()
##' list(...)
##' list(list(...))
##' list(list(...), list(...), ...)
##' @param selection the list structure to specify the month selection
##' @param climate the climate data as returned by make_pmat
##' @param check_2 should there be a check for less than 2 variables
##' in the design matrix?
##' @return a data.frame
tc_design <- function(selection, climate, check_2 = TRUE) {
  ## is it a list?
  if (!is.list(selection)) {
    stop("Please supply information about independent variables as list.")
  }
  ## check for nested list structure
  n <- length(as.list(selection))
  if (n == 0) {
    ## empty list, default parameters are taken
    out <- eval_selection(climate, selection)
  }
  if (n == 1) {
    ## could be list("mean") or list(list(...))
    if (is.list(selection[[1]])) {
      ## list(list(...))
      out <- eval_selection(climate, selection[[1]])
    } else {
      ## list("mean")
      out <- eval_selection(climate, selection)
    }
  }
  if (n > 1) {
    ## list(list(...), list(...), ...) or list("mean", 1:10)
    if (is.list(selection[[1]])) {
      ## everything has to be a list
      if (all(sapply(selection, is.list))) {
        ## list(list(...), list(...), ...)
        OUT <- list()
        for (i in 1:n) {
          OUT[[i]] <- eval_selection(climate, selection[[i]])
        }
        out <- list()
        out$month <- list()
        out$month$match <- sapply(OUT, function(x) {
          x$month$match })
        out$month$names <- sapply(OUT, function(x) {
          x$month$names})
        out$method <- sapply(OUT, function(x) {
          x$method })
        out$param <- sapply(OUT, function(x) {
          x$param })
        out$aggregate <- as.data.frame(lapply(OUT, function(x) {
          x$aggregate }))
        out$names <- unlist(as.vector(sapply(OUT, function(x) {
          x$names })))
      } else {
        stop("You may not mix modifiers/lists and other objects for parameter specifications.")
      }
    } else {
      ## everything has to be not a list
      if (all(sapply(selection, function(x) { ifelse(is.list(x), FALSE,
                                                     TRUE)}))) {
        ## list("mean", 1:10)
        out <- eval_selection(climate, selection)
      } else {
        stop("You may not mix lists and other objects for parameter specifications.")
      }
    }
  }

  ## avoid duplicate parameters
  dupes <- duplicated(out$names)
  if (any(dupes)) {
    out$aggregate <- out$aggregate[,!dupes]
    out$names <- out$names[!dupes]
  }

  ## throw error when we have only one variable left; point the user to using
  ## lm() instead
  if (dim(out$aggregate)[2] < 2 & check_2)
    stop("You supplied only one climate variable for calibration. treeclim needs at least two. Consider using lm() in this case. Thanks.")

  ## reorder and add pretty names for plotting required: months as numeric values
  ## from 1:25 (25 is for month aggregations); we do this _after_ removing the
  ## potential duplicates, which is safer, but also results in a somewhat lengthy
  ## code...

  get_month_numeric <- function(m) {
    switch(m,
           "jan" = 1,
           "feb" = 2,
           "mar" = 3,
           "apr" = 4,
           "may" = 5,
           "jun" = 6,
           "jul" = 7,
           "aug" = 8,
           "sep" = 9,
           "oct" = 10,
           "nov" = 11,
           "dec" = 12)
  }
  
  labels <- vars <- out$names
  .months <- numeric(length(labels))
  n <- length(out$names)
  for (i in 1:n) {
    split <- strsplit(out$names[i], "\\.")[[1]]
    mmatches <- !is.na((match(split, c("curr", "prev"))))
    fmatch <- which(mmatches)[1] - 1
    vars[i] <- paste(split[1:fmatch], collapse = ".")
    fmonth_char <- split[which(mmatches)[1] + 1]
    if (sum(mmatches) > 1) {
      .months[i] <- 25                  # 25: unique number for
                                        # aggregates!
      ## long variable name
      season <- split[mmatches]
      smatches <- c(FALSE, mmatches[-length(mmatches)])
      months <- split[smatches]
      ## get first and last month and season
      fseason <- season[1]
      fmonth <- months[1]
      if (fseason == "prev") {
        fmonth <- tolower(fmonth)
      } else {
        fmonth <- toupper(fmonth)
      }
      lseason <- tail(season, 1)
      lmonth <- tail(months, 1)
      if (lseason == "prev") {
        lmonth <- tolower(lmonth)
      } else {
        lmonth <- toupper(lmonth)
      }
      labels[i] <- paste(fmonth, "...", lmonth, sep = "")
    } else {
      .months[i] <- get_month_numeric(fmonth_char)
      ## short variable name
      tseason <- split[mmatches]
      tmonth <- split[which(mmatches)[1] + 1]
      if (tseason == "prev") {
        labels[i] <- tolower(tmonth)
      } else {
        labels[i] <- toupper(tmonth)
        .months[i] <- .months[i] + 12
      }
    }
  }

  pretty_names <- data.frame(
    month_no = .months,
    varname = vars,
    month_label = labels
    )

  pretty_order <- with(pretty_names, order(varname, month_no))
  
  if (dim(out$aggregate)[2] > 1) {
    out$aggregate <- out$aggregate[,pretty_order]
  }
  names(out$aggregate) <- paste("X", 1:dim(out$aggregate)[2],
                                sep = "")
  
  out$names <- out$names[pretty_order]

  out$pretty_names <- pretty_names[pretty_order,]
  out$pretty_names$id <- rownames(out$pretty_names) <- 1:dim(out$pretty_names)[1]
  
  class(out) <- list("tc_design", "list")
  out
}
