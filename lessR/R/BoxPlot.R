BoxPlot <-
function(x=NULL, data=mydata, n.cat=getOption("n.cat"),
    Rmd=NULL,

    color.fill=getOption("color.fill.bar"),
    color.stroke=getOption("color.stroke.bar"), 
    color.bg=getOption("color.bg"),
    color.grid=getOption("color.grid"),
    color.box=getOption("color.box"),

    cex.axis=0.75, color.axis="gray30",
    xlab=NULL, main=NULL, sub=NULL, digits.d=NULL,

    rotate.values=0, offset=0.5,

    horiz=TRUE, add.points=FALSE,

    quiet=getOption("quiet"),
    pdf.file=NULL, pdf.width=5, pdf.height=5, 
    fun.call=NULL, ...) {


  if (is.null(fun.call)) fun.call <- match.call()

  for (i in 1:length(color.fill))
    if (color.fill[i] == "off") color.fill[i] <- "transparent"
  for (i in 1:length(color.stroke))
    if (color.stroke[i] == "off") color.stroke[i] <- "transparent"
  if (color.bg == "off") color.bg <- "transparent"
  if (color.grid == "off" ) color.grid <- "transparent"
  if (color.box == "off") color.box <- "transparent"

  if (getOption("colors") == "gray") color.stroke <- "black"
  if (getOption("colors") == "gray.black") color.stroke <- getOption("color.stroke.pt")

  if (!is.null(pdf.file))
    if (!grepl(".pdf", pdf.file)) pdf.file <- paste(pdf.file, ".pdf", sep="")

  dots <- list(...)  # check for deprecated parameters
  if (length(dots) > 0) {
    for (i in 1:length(dots)) {
      if (names(dots)[i] == "knitr.file") {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
          "knitr.file  no longer used\n",
          "Instead use  Rmd  for R Markdown file\n\n")
      }
    }
    if (substr(names(dots)[i], 1, 4) == "col.") {
      cat("\n"); stop(call.=FALSE, "\n","------\n",
        "options that began with the abbreviation  col  now begin with  ",
        "color \n\n")
    }
  }

  # get actual variable name before potential call of data$x
  x.name <- deparse(substitute(x))
  options(xname = x.name)

  df.name <- deparse(substitute(data))
  options(dname = df.name)

  pdf.nm <- FALSE
  if (!missing(pdf.file)) pdf.nm <- TRUE

# -----------------------------------------------------------
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


# ---------------
# do the analysis

  go.pdf <- FALSE
  if (pdf.nm || ncol(data) > 1) go.pdf <- TRUE

  for (i in 1:ncol(data)) {

    nu <- length(unique(na.omit(data[,i])))

    x.name <- names(data)[i]
    options(xname = x.name)

    if (is.numeric(data[,i])) {
      # let 1 variable go through, even if num.cat
      if (ncol(data) == 1  ||  !.is.num.cat(data[,i], n.cat)) {

      pdf.fnm <- .pdfname("BoxPlot", x.name, go.pdf, pdf.nm, pdf.file)
     .opendev(pdf.fnm, pdf.width, pdf.height)

      stuff <- .bx.main(data[,i], color.fill, color.stroke, color.bg, color.grid,
         color.box, cex.axis, color.axis, rotate.values, offset, 
         horiz, add.points, xlab, main, sub, digits.d, quiet, ...)
      txsts <- stuff$tx
      if (length(txsts)==0) txsts <- ""

      txotl <- ""
      if (!quiet) {
        txotl <- .outliers(data[,i])
        if (length(txotl)==0) txotl <- "No outliers"
      }

      if (ncol(data) > 1) {  # for a variable range, just the text output
        class(txsts) <- "out_piece"
        class(txotl) <- "out_piece"
        output <- list(out_stats=txsts, out_outliers=txotl)
        class(output) <- "out_all"
        print(output)
      }

      if (go.pdf) {
        dev.off()
        if (!quiet) .showfile(pdf.fnm, "Box Plot")
      }

    }  # nu > n.cat
    else
      .ncat("Box Plot", x.name, nu, n.cat)

    }  # is.numeric(data[,i])
  }  # for


  dev.set(which=2)  # reset graphics window for standard R functions

  if (ncol(data)==1) {

    # R Markdown
    txkfl <- ""
    if (!is.null(Rmd)) {
      if (!grepl(".Rmd", Rmd)) Rmd <- paste(Rmd, ".Rmd", sep="")
      txknt <- .dist.Rmd(x.name, df.name, fun.call, digits.d)
      cat(txknt, file=Rmd, sep="\n")
      txkfl <- .showfile2(Rmd, "R Markdown instructions")
    }
 
    class(txsts) <- "out_piece"
    class(txotl) <- "out_piece"
    class(txkfl) <- "out_piece"

    output <- list(type="BoxPlot",
      call=fun.call,
      out_stats=txsts, out_outliers=txotl, out_file=txkfl,
      n=stuff$n, n.miss=stuff$n.miss, min=stuff$mn, lower_whisker=stuff$lw,
      lower_hinge=stuff$lh, median=stuff$md, upper_hinge=stuff$uh,
      upper_whisker=stuff$uw, max=stuff$mx, IQR=stuff$IQR)

    class(output) <- "out_all"

    return(output)

  }

}


