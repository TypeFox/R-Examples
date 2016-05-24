Density <-
function(x, data=mydata, n.cat=getOption("n.cat"), 

    bw="nrd0", type=c("both", "general", "normal"),
    bin.start=NULL, bin.width=NULL,

    Rmd=NULL, digits.d=NULL,

    color.fill=getOption("color.fill.pt"),
    color.bg=getOption("color.bg"),
    color.grid=getOption("color.grid"),
    color.box=getOption("color.box"),

    color.nrm="black", color.gen="black",
    color.fill.nrm=NULL, color.fill.gen=NULL,

    cex.axis=0.75, color.axis="gray30",

    rotate.values=0, offset=0.5,

    x.pt=NULL, xlab=NULL, main=NULL, sub=NULL, y.axis=FALSE, 
    x.min=NULL, x.max=NULL, band=FALSE, 

    quiet=getOption("quiet"),
    pdf.file=NULL, pdf.width=5, pdf.height=5,
    fun.call=NULL, ...) {


  if (is.null(fun.call)) fun.call <- match.call()

  for (i in 1:length(color.fill))
    if (color.fill[i] == "off") color.fill[i] <- "transparent"
  if (color.bg == "off") color.bg <- "transparent"
  if (color.grid == "off" ) color.grid <- "transparent"
  if (color.box == "off") color.box <- "transparent"

  if (!is.null(pdf.file))
    if (!grepl(".pdf", pdf.file)) pdf.file <- paste(pdf.file, ".pdf", sep="")

  dots <- list(...)  # check for deprecated parameters
  if (length(dots) > 0) {
    for (i in 1:length(dots)) {
      old.nm <- c("col.fill", "col.bg", "col.grid", "col.box", "col.nrm",
                  "col.gen", "col.fill.nrm", "col.fill.gen", "col.axis")
      if (names(dots)[i] %in% old.nm) {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
          "options that began with the abbreviation  col  now begin with  ",
          "color \n\n")
      }
      if (names(dots)[i] == "knitr.file") {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
          "knitr.file  no longer used\n",
          "Instead use  Rmd  for R Markdown file\n\n")
      }
    }
  }

  clr <- getOption("colors")  # color theme not used except for monochrome 

  if (missing(color.fill))
    if (.Platform$OS == "windows")
      color.fill <- "gray80"
    else
      color.fill <- "gray86"

  if (missing(color.bg)) color.bg <- "ghostwhite"

  if (missing(color.fill.nrm))
      color.fill.nrm <- rgb(80,150,200, alpha=70, maxColorValue=255)

  if (missing(color.fill.gen))
      color.fill.gen <- rgb(250,210,230, alpha=70, maxColorValue=255)

  if (clr == "gray" || clr == "gray.black") {
    color.fill.nrm <- "transparent"
    color.fill.gen <- "transparent"
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

      pdf.fnm <- .pdfname("Density", x.name, go.pdf, pdf.nm, pdf.file)
     .opendev(pdf.fnm, pdf.width, pdf.height)

      gl <- .getlabels()
      x.name <- gl$xn; x.lbl <- gl$xl;
      y.name <- gl$yn; y.lbl <- gl$yl
      if (!quiet)
        ttlns <- .title2(x.name, y.name, x.lbl, y.lbl, TRUE)
      else
        ttlns <- ""
 
      stuff <- .dn.main(data[,i], bw, type, bin.start, bin.width, 
            color.fill, color.bg, color.grid, color.box, color.nrm, color.gen,
            color.fill.nrm, color.fill.gen, 
            cex.axis, color.axis, rotate.values, offset, 
            x.pt, xlab, main, sub, y.axis, x.min, x.max, band, quiet, ...)
      txdst <- stuff$tx
      if (length(txdst)==0) txdst <- ""

      if (ncol(data) > 1) {  # for a variable range, just print text output
        class(ttlns) <- "out_piece"
        class(txdst) <- "out_piece"
        output <- list(out_title=ttlns, out_stats=txdst)
        class(output) <- "out_all"
        print(output)
      }

      if (go.pdf) {
        dev.off()
        if (!quiet) .showfile(pdf.fnm, "density curve")
      }

    }  # nu > n.cat
    else
      .ncat("Density curve", x.name, nu, n.cat)

    }  # is.numeric(data[,i])
  }  # for


  dev.set(which=2)  # reset graphics window for standard R functions

  if (ncol(data)==1  &&  nu>n.cat) {

    # R Markdown
    if (is.null(digits.d)) {
      dig.dec <- .max.dd(data[,i]) + 1
      if (dig.dec == 1) dig.dec <- 2
    }
    else dig.dec <- digits.d
    options(digits.d=dig.dec)

    txkfl <- ""
    if (!is.null(Rmd)) {
      if (!grepl(".Rmd", Rmd)) Rmd <- paste(Rmd, ".Rmd", sep="")
      txknt <- .dist.Rmd(x.name, df.name, fun.call, digits.d)
      cat(txknt, file=Rmd, sep="\n")
      txkfl <- .showfile2(Rmd, "R Markdown instructions")
    }
 
    class(txdst) <- "out_piece"
    class(txkfl) <- "out_piece"

    output <- list(type="Density",
      out_title=ttlns, out_stats=txdst, out_file=txkfl,
      bw=stuff$bw, n=stuff$n, n.miss=stuff$n.miss, W=stuff$W,
         pvalue=stuff$pvalue)

    class(output) <- "out_all"

    return(output)

  }

}
