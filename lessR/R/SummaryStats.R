SummaryStats <-
function(x=NULL, by=NULL, data=mydata, n.cat=getOption("n.cat"), 
    digits.d=NULL, brief=getOption("brief"), ...)  {


  # get variable name before potential call of data$x
  x.name <- deparse(substitute(x))  # could be a vars list
  options(xname = x.name)

  df.name <- deparse(substitute(data))
  options(dname = df.name)

  is.df <- FALSE  # is data frame


  # evaluate by, i.e., get y.call
  if (!missing(by)) {

    # get actual variable name before potential call of data$x
    y.name <- deparse(substitute(by)) 
    options(yname = y.name)

    # see if y exists from a function call
    # indicate a function call with sys.frame returns larger than 1 
    if (exists(y.name, where=parent.frame(n=1)) && sys.nframe() > 1) 
      in.call <- TRUE else in.call <- FALSE

    # get conditions and check for data existing
    if (!in.call) {
      xs <- .xstatus(y.name, df.name)
      in.global <- xs$ig 
    }
    else in.global <- FALSE

    # see if var exists in data frame, if x not in Global Env or function call 
    if (!in.global && !in.call) .xcheck(y.name, df.name, data)

    if (!in.global) y.call <- eval(substitute(data$by))
    else {  # vars that are function names get assigned to global
      y.call <- by
      if (is.function(y.call)) y.call <- eval(substitute(data$by))
    }
  }
  else
    y.call <- NULL


# -----------------------------------------------------------------
# establish if a data frame, if not then identify variable(s)

  if (!missing(x)) {

    if (!exists(x.name, where=.GlobalEnv)) {  # x not in global env, in df
      .nodf(df.name)  # check to see if data frame container exists 
      .xcheck(x.name, df.name, data)  # var in df?, vars lists not checked
      all.vars <- as.list(seq_along(data))  # even if only a single var
      names(all.vars) <- names(data)  # all data in data frame
      x.col <- eval(substitute(x), envir=all.vars)  # col num selected vars
      if (!("list" %in% class(data))) {
        data <- data[, x.col]
        if (length(x.col) == 1) {  # x is 1 var
          data <- data.frame(data)
          names(data) <- x.name
         }
      }
      else {  # class of data is "list"
        data <- data.frame(data[[x.col]])
        names(data) <- x.name
      }
    }

    else { # x is in the global environment (vector or data frame)
      if (is.data.frame(x))  # x a data frame
        data <- x
      else {  # x a vector in global
        if (!is.function(x))
          data <- data.frame(x)  # x is 1 var
        else
          data <- data.frame(eval(substitute(data$x)))  # x is 1 var
        names(data) <- x.name
      }
    }

  }


# -----------------------------------------------------
# data is now set
# do the analysis for a single variable or a data frame 

  if (ncol(data) == 1) {
    x.call <- data[,1]

    num.cat <- .is.num.cat(x.call, n.cat)
    nu <- length(unique(na.omit(x.call)))

    if (!num.cat && is.integer(x.call) && nu <= n.cat) {
      cat("\n")
      cat(">>> ", x.name, " has only only ", nu, " unique ",
          "integer values, but not equally spaced,\n",
          "      so treat as numerical in this analysis\n",
          "    Maybe convert to an R factor to treat as categorical\n\n",
          sep="")
    }

    if (num.cat || is.character(x.call)) {
      x.call <- as.factor(x.call)
      .ncat("summary statistics", x.name, nu, n.cat)
    }
  }

  if (ncol(data) > 1)
    .ss.data.frame(data, n.cat, brief, ...) 

  else if (!is.factor(x.call)) {
    sk <- NA; kt <- NA; q1 <- NA; q3 <- NA;  qr <- NA;
    stuff <- .ss.numeric(x.call, y.call, digits.d, brief, ...)
    txsts <- stuff$tx
    txotl <- .outliers(x.call)
  }

  # ordered factors have two attributes, "ordered" and "factor"
  else if (is.factor(x.call)) {

    gl <- .getlabels(xlab=NULL, ylab=NULL, main=NULL, cex.lab=NULL)
    x.name <- gl$xn; x.lab <- gl$xb; x.lbl <- gl$xl
    y.name <- gl$yn; y.lab <- gl$yb; y.lbl <- gl$yl

    if (is.factor(x.call)) 
      stuff <- .ss.factor(x.call, y.call, brief, digits.d,
                        x.name, y.name, x.lbl, y.lbl, ...) 
    else if (is.character(x.call))
      if (nlevels(factor(x.call)) < length(x.call)) { 
        stuff <- .ss.factor(factor(x.call), by, brief, digits.d=NULL,
                        x.name, y.name, x.lbl, y.lbl, ...) 
      }
      else cat("\n Appears to contain unique Names or IDs", "\n")

    n.dim <- stuff$n.dim
    if (n.dim == 1) {
      txttl <- stuff$title
      txsts <- stuff$counts
      txchi <- stuff$chi
      frq <- stuff$freq
      prp <- stuff$prop
    }
    else if (n.dim == 2) {
      txttl <- stuff$txttl
      txXV <- stuff$txXV
      txfrq <- stuff$txfrq
      txprp <- stuff$txprp
      txcol <- stuff$txcol
      txrow <- stuff$txrow
    }

  }

  else {
    cat("\n"); stop(call.=FALSE, "\n","------\n",
        "The variable to be analyzed must be numeric or a factor, or have\n",
        "character values that can be converted to a factor.\n")
  }


  if (ncol(data) == 1) { 

    if (!is.factor(x.call)) {
      class(txsts) <- "out_piece"
      class(txotl) <- "out_piece"

      output <- list(out_stats=txsts, out_outliers=txotl, n=stuff$n,
           n.miss=stuff$n.miss, mean=stuff$m, sd=stuff$s, skew=stuff$sk,
           kurtosis=stuff$kt, min=stuff$mn, quartile1=stuff$q1,
           median=stuff$md, quartile3=stuff$q3, max=stuff$mx, IQR=stuff$qr)
    }

    else if (is.factor(x.call)) {
      if (n.dim == 1) {
        class(txttl) <- "out_piece"
        class(txsts) <- "out_piece"
        class(txchi) <- "out_piece"
        class(frq) <- "out_piece"
        class(prp) <- "out_piece"
        output <- list(out_title=txttl, out_stats=txsts, out_chi=txchi, freq=frq, prop=prp)
      }
      else if (n.dim == 2) {
        class(txttl) <- "out_piece"
        class(txfrq) <- "out_piece"
        if (brief) {
          output <- list(out_title=txttl, out_freq=txfrq, out_chi=txXV)
        }
        else {
          class(txXV) <- "out_piece"
          class(txprp) <- "out_piece"
          class(txcol) <- "out_piece"
          class(txrow) <- "out_piece"
          output <- list(out_title=txttl, out_freq=txfrq, out_XV=txXV,
             out_prop=txprp, out_colsum=txcol, out_rowsum=txrow)
        }
      }
    }

    class(output) <- "out_all"

    return(output)

  }

}
