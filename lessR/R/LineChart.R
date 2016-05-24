LineChart <-
function(x, data=mydata, n.cat=getOption("n.cat"), type=NULL, 

         color.fill=getOption("color.fill.bar"), 
         color.stroke=getOption("color.stroke.pt"),
         color.bg=getOption("color.bg"),
         color.grid=getOption("color.grid"),
         color.box=getOption("color.box"),
         color.line=getOption("color.stroke.pt"),

         color.area=NULL, 

         shape.pts=21, cex.axis=0.75, color.axis="gray30",

         rotate.values=0, offset=.5,

         xy.ticks=TRUE, line.width=1,
         xlab=NULL, ylab=NULL, main=NULL, sub=NULL, cex=NULL,

         time.start=NULL, time.by=NULL, time.reverse=FALSE,

         center.line=c("default", "mean", "median", "zero", "off"),

         show.runs=FALSE, quiet=getOption("quiet"),
         pdf.file=NULL, pdf.width=5, pdf.height=5, ...) {


  center.line <- match.arg(center.line)

  for (i in 1:length(color.fill))
    if (color.fill[i] == "off") color.fill[i] <- "transparent"
  for (i in 1:length(color.stroke))
    if (color.stroke[i] == "off") color.stroke[i] <- "transparent"
  if (color.bg == "off") color.bg <- "transparent"
  if (color.grid == "off" ) color.grid <- "transparent"
  if (color.box == "off") color.box <- "transparent"
  if (color.line == "off") color.box <- "transparent"
  if (!is.null(color.area)) if (color.area == "off") color.area <- "transparent"

  if (!is.null(pdf.file))
    if (!grepl(".pdf", pdf.file)) pdf.file <- paste(pdf.file, ".pdf", sep="")

  dots <- list(...)  # check for deprecated parameters
  if (length(dots) > 0) {
    for (i in 1:length(dots)) {
      old.nm <- c("col.fill", "col.stroke", "col.bg", "col.grid", "col.box",
                  "col.line", "col.axis", "col.area")
      if (names(dots)[i] %in% old.nm) {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
          "options that began with the abbreviation  col  now begin with  ",
          "color \n\n")
      }
    }
  }

  # get actual variable name before potential call of data$x
  x.name <- deparse(substitute(x))
  options(xname = x.name)
  options(yname = x.name)  # for .lc.main, which uses y as the var

  df.name <- deparse(substitute(data))
  options(dname = df.name)

  pdf.nm <- FALSE
  if (!missing(pdf.file)) pdf.nm <- TRUE

  if (!quiet && !show.runs)
    cat("[To view the individual runs: show.runs=TRUE]\n")

# -----------------------------------------------------------
# establish if a data frame, if not then identify variable(s)

  if (!missing(x)) {
    if (!exists(x.name, where=.GlobalEnv)) {  # x not in global env, in df
      .nodf(df.name)  # check to see if data frame container exists 
      .xcheck(x.name, df.name, data)  # see if var in df, vars lists not checked
      vars.list <- as.list(seq_along(data))
      names(vars.list) <- names(data)
      x.col <- eval(substitute(x), envir=vars.list)  # col num of each var
      if (!("list" %in% class(data))) {
        data <- data[, x.col]
        if (length(x.col) == 1) {
          data <- data.frame(data)  # x is 1 var
          names(data) <- x.name
         }
      }
      else {
        data <- data.frame(data[[x.col]])
        names(data) <- x.name
      }
    }
    else { # x is in the global environment (vector or data frame)
      if (is.data.frame(x))  # x a data frame
        data <- x
      else {  # x a vector in global
        data <- data.frame(x)  # x is 1 var
        names(data) <- x.name
      }
    }
  }


# ---------------
# do the analysis

  go.pdf <- FALSE
  if (pdf.nm || ncol(data) > 1) go.pdf <- TRUE

  for (i in 1:ncol(data)) {
    cat("\n")

    nu <- length(unique(na.omit(data[,i])))

    x.name <- names(data)[i]
    options(xname = x.name)

    if (is.numeric(data[,i])) {
      # let 1 variable go through, even if num.cat
      if (ncol(data) == 1  ||  !.is.num.cat(data[,i], n.cat)) {

      pdf.fnm <- .pdfname("LC", x.name, go.pdf, pdf.nm, pdf.file)
     .opendev(pdf.fnm, pdf.width, pdf.height)

      .lc.main(data[,i], type,
         color.line, color.area, color.stroke, color.fill, shape.pts,
         color.grid, color.box, color.bg, cex.axis, color.axis,
         rotate.values, offset, xy.ticks,
         line.width, xlab, ylab, main, sub, cex,
         time.start, time.by, time.reverse, 
         center.line, show.runs, quiet, ...)

      if (go.pdf) {
        dev.off()
        if (!quiet) .showfile(pdf.fnm, "Line Chart")
      }

    }  # nu > n.cat
    else
      .ncat("Line Chart", x.name, nu, n.cat)

    }  # is.numeric(data[,i])
  }  # for

}
