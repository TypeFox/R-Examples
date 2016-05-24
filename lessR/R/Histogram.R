Histogram <-
function(x=NULL, data=mydata, n.cat=getOption("n.cat"), Rmd=NULL,

    color.fill=getOption("color.fill.bar"), 
    color.stroke=getOption("color.stroke.bar"),
    color.bg=getOption("color.bg"),
    color.grid=getOption("color.grid"),
    color.box=getOption("color.box"),

    color.reg="snow2", over.grid=FALSE,
    cex.axis=0.75, color.axis="gray30", 

    rotate.values=0, offset=0.5,

    breaks="Sturges", bin.start=NULL, bin.width=NULL, bin.end=NULL,

    prop=FALSE, cumul=c("off", "on", "both"), hist.counts=FALSE,
    digits.d=NULL, xlab=NULL, ylab=NULL, main=NULL, sub=NULL,

    quiet=getOption("quiet"),
    pdf.file=NULL, pdf.width=5, pdf.height=5,
    fun.call=NULL, ...) {


  if (is.null(fun.call)) fun.call <- match.call()

  # limit actual argument to alternatives, perhaps abbreviated
  cumul <- match.arg(cumul)

  for (i in 1:length(color.fill))
    if (color.fill[i] == "off") color.fill[i] <- "transparent"
  for (i in 1:length(color.stroke))
    if (color.stroke[i] == "off") color.stroke[i] <- "transparent"
  if (color.bg == "off") color.bg <- "transparent"
  if (color.grid == "off" ) color.grid <- "transparent"
  if (color.box == "off") color.box <- "transparent"

  if (!is.null(pdf.file))
    if (!grepl(".pdf", pdf.file)) pdf.file <- paste(pdf.file, ".pdf", sep="")

  kf <- FALSE
  lbls <- FALSE
  dots <- list(...)  # check for deprecated parameters
  if (length(dots) > 0) {
    old.nm <- c("col.fill", "col.stroke", "col.bg", "col.grid", "col.box",
                "col.reg", "col.axis") 
    for (i in 1:length(dots)) {
      if (names(dots)[i] %in% old.nm) {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
          "options that began with the abbreviation  col  now begin with  ",
          "color \n\n")
      }
      if (names(dots)[i] == "knitr.file") kf <- TRUE 
      if (names(dots)[i] == "labels") lbls <- TRUE 
    }
  }

  if (lbls) {
    cat("\n"); stop(call.=FALSE, "\n","------\n",
      "labels  has multiple definitions in R\n",
      "Instead use  hist.counts  to get the bar labels displayed\n\n")
  }

  if (kf) {
    cat("\n"); stop(call.=FALSE, "\n","------\n",
      "knitr.file  no longer used\n",
      "Instead use  Rmd  for R Markdown file\n\n")
  }
        
  if (is.numeric(breaks) && !is.null(bin.start)) { 
    cat("\n"); stop(call.=FALSE, "\n","------\n",
      "Choose only one option to specify a start value.\n",
      "Either choose the option  breaks  or the option  bin.start.\n\n")
  }

  # get actual variable name before potential call of data$x
  x.name <- deparse(substitute(x))  # could be a list of var names
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
          if (!is.numeric(data)) { 
            cat("\n"); stop(call.=FALSE, "\n","------\n",
              "A histogram is only computed from a numeric variable\n",
              "For the frequencies of a categorical variable:\n\n",
              "  Plot(", x.name, ", topic=\"count\")\n",
              "or\n",
              "  BarChart(", x.name, ")\n\n", sep="")
          }
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


  if (ncol(data) > 1) {
    sug <- getOption("suggest")
    options(suggest = FALSE)
  }

  for (i in 1:ncol(data)) {

    nu <- length(unique(na.omit(data[,i])))

    x.name <- names(data)[i]
    options(xname = x.name)

    if (is.numeric(data[,i])) {
      # let 1 variable go through, even if num.cat
      if (ncol(data) == 1  ||  !.is.num.cat(data[,i], n.cat)) {

      pdf.fnm <- .pdfname("Hist", x.name, go.pdf, pdf.nm, pdf.file)
     .opendev(pdf.fnm, pdf.width, pdf.height)

      txss <- ""
      if (!quiet) {
        ssstuff <- .ss.numeric(data[,i], digits.d=digits.d, brief=TRUE)
        txss <- ssstuff$tx
      }

      # nothing returned if quiet=TRUE
      stuff <- .hst.main(data[,i], color.fill, color.stroke, color.bg,
          color.grid, color.box, color.reg,
          over.grid, cex.axis, color.axis, rotate.values, offset,
          breaks, bin.start, bin.width,
          bin.end, prop, hist.counts, cumul, xlab, ylab, main, sub,
          quiet, fun.call=fun.call, ...)
      txsug <- stuff$txsug
      if (is.null(txsug)) txsug <- ""
      txdst <- stuff$ttx
      if (is.null(txdst)) txdst <- ""

      txotl <- ""
      if (!quiet) {
        txotl <- .outliers(data[,i])
        if (length(txotl)==0) txotl <- "No outliers"
      }

      if (ncol(data) > 1) {  # for a variable range, print the text output
        class(txss) <- "out_piece"
        class(txdst) <- "out_piece"
        class(txotl) <- "out_piece"
        output <- list(out_ss=txss, out_freq=txdst, out_outliers=txotl)
        class(output) <- "out_all"
        print(output)
      }

      if (go.pdf) {
        dev.off()
        if (!quiet) .showfile(pdf.fnm, "Histogram")
      }

    }  # nu > n.cat
    else
      .ncat("Histogram", x.name, nu, n.cat)

    }  # is.numeric(data[,i])
  }  # for

  if (ncol(data) > 1) options(suggest = sug)

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
 
    class(txsug) <- "out_piece"
    class(txss) <- "out_piece"
    class(txdst) <- "out_piece"
    class(txotl) <- "out_piece"
    class(txkfl) <- "out_piece"

    output <- list(type="Histogram",
      call=fun.call, 
      out_suggest=txsug, out_ss=txss, out_outliers=txotl, out_freq=txdst,
      out_file=txkfl,
      bin_width=stuff$bin.width, n_bins=stuff$n.bins,
      breaks=stuff$breaks,
      mids=stuff$mids, counts=stuff$counts, prop=stuff$prop,
      counts_cumul=stuff$counts_cum, prop_cumul=stuff$prop_cum)

    class(output) <- "out_all"

    return(output)

  }

}
