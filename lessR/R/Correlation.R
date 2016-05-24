Correlation <-
function(x, y, data=mydata, # x can be a data frame, or variables in a data frame
         miss=c("pairwise", "listwise", "everything"),
         show.n=NULL, brief=FALSE, 
         digits.d=NULL, graphics=FALSE,
         main=NULL, bottom=3, right=3,
         pdf=FALSE, pdf.width=5, pdf.height=5, ...) {


  miss <- match.arg(miss)
  
  x.name <- deparse(substitute(x))
  options(xname = x.name)

  is.df <- FALSE  # is data frame

  if (missing(y)) {  # is a data frame or a list of variables

    if (missing(x)) {
      x.name <- ""  # in case x is missing, i.e., data frame mydata
      is.df <- TRUE
      if (missing(data)) data <- eval(substitute(mydata))
    }

    else if ( (!grepl(":", x.name) && !grepl(",", x.name)) ) {
      if (is.data.frame(x)) {  # specified data name
        if (exists(x.name, where=.GlobalEnv)) {
          data <- x
          is.df <- TRUE
        }
      }
    } 

    else {

      all.vars <- as.list(seq_along(data))
      names(all.vars) <- names(data)

      # proceed here only if x.name is in data or is a list
      if ( (x.name %in% names(all.vars)) || grepl(":", x.name) || grepl(",", x.name) ) {
        x.col <- eval(substitute(x), envir=all.vars, enclos=parent.frame())

        # if x is a variable list, create subset data frame
        if (length(x.col) > 1) {
          data <- data[, x.col]
          is.df <- TRUE
        }
      }
    }
  }  # end missing y


  if (!is.df) {

    dname <- deparse(substitute(data))
    options(dname = dname)

    # get conditions and check for data existing
    xs <- .xstatus(x.name, dname)
    is.frml <- xs$ifr
    in.global <- xs$ig 

    # see if the variable exists in data frame, if x not in Global Env 
    if (!in.global) .xcheck(x.name, dname, data)

    if (in.global) x.call <- x else x.call <- eval(substitute(data$x))

    # evaluate y
    if (!missing(y)) {

      # get actual variable name before potential call of data$x
      y.name <- deparse(substitute(y)) 
      options(yname = y.name)

      # see if y exists from a function call
      # indicate a function call with sys.frame returns larger than 1 
      #if (exists(y.name, where=parent.frame(n=1)) && sys.nframe() > 1) 
        #in.call <- TRUE else in.call <- FALSE
      in.call <- FALSE

      # get conditions and check for data existing
      if (!in.call) {
        xs <- .xstatus(y.name, dname)
        in.global <- xs$ig 
      }
      else in.global <- FALSE

      # see if var exists in data frame, if x not in Global Env or function call 
      if (!in.global && !in.call) .xcheck(y.name, dname, data)

      if (in.global) y.call <- y 
      else y.call <- eval(substitute(data$y))
    }

    else
      y.call <- NULL

  }  # x not data frame


  if (is.df) { 
    if (is.null(show.n))
      if (nrow(data) <= 15) show.n <- TRUE else show.n <- FALSE
    stuff <- .cr.data.frame(data, miss, show.n, digits.d,
                   graphics, main, bottom, right, 
                   pdf, pdf.width, pdf.height, ...) 

    txbck <- stuff$txb
    txmis <- stuff$txm
    txcor <- stuff$txc

    class(txbck) <- "out_piece"
    class(txmis) <- "out_piece"
    class(txcor) <- "out_piece"

    output <- list(
      out_background=txbck, out_missing=txmis, out_cor=txcor,
      cors=stuff$cors)
  }

  else {
  
    if (pdf) {
      cat("\n");   warning(call.=FALSE, "\n","------\n",
        " To produce a scatter plot, use the function:  ScatterPlot\n",
        " No graphics produced here and yet pdf is specified.\n\n",
         sep="")
     }

    stuff <- .cr.main(x.call, y.call, brief, ...) 
    txbck <- stuff$txb
    txdsc <- stuff$txd
    txinf <- stuff$txi

    class(txbck) <- "out_piece"
    class(txdsc) <- "out_piece"
    class(txinf) <- "out_piece"

    output <- list(out_background=txbck, out_describe=txdsc, out_inference=txinf,
      r=stuff$r, tvalue=stuff$tvalue, df=stuff$df, pvalue=stuff$pvalue,
      lb=stuff$lb, ub=stuff$ub)
  }

    class(output) <- "out_all"
    return(output)
}
