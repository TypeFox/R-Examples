BarChart <-
function(x=NULL, by=NULL, data=mydata, n.cat=getOption("n.cat"), 

         color.fill=getOption("color.fill.bar"),
         color.stroke=getOption("color.stroke.bar"),
         color.bg=getOption("color.bg"),
         color.grid=getOption("color.grid"),
         color.box=getOption("color.box"),

         colors=c("rainbow", "terrain", "heat"),

         horiz=FALSE, over.grid=FALSE, addtop=0.05,
         gap=NULL, proportion=FALSE,
         
         xlab=NULL, ylab=NULL, main=NULL,
         cex.axis=0.75, color.axis="gray30",
         value.labels=NULL, rotate.values=0, offset=0.5,

         beside=FALSE, color.low=NULL, color.hi=NULL, count.labels=NULL,

         legend.title=NULL, legend.loc="right.margin", legend.labels=NULL,
         legend.horiz=FALSE, 

         quiet=getOption("quiet"),
         pdf.file=NULL, pdf.width=5, pdf.height=5, ...)  {


  if (missing(colors)) 
    colors <- getOption("colors")
  else
    colors <- match.arg(colors)

  if (missing(color.stroke))  # default black border unless dark bg
    if (sum(col2rgb(color.bg))/3 > 80) color.stroke <- "black"

  if (!is.null(color.fill)) {
    for (i in 1:length(color.fill))
      if (color.fill[i] == "off") color.fill[i] <- "transparent"
  }
  for (i in 1:length(color.stroke))
    if (color.stroke[i] == "off") color.stroke[i] <- "transparent"
  if (color.bg == "off") color.bg <- "transparent"
  if (color.grid == "off" ) color.grid <- "transparent"
  if (color.box == "off") color.box <- "transparent"

  if (!is.null(pdf.file))
    if (!grepl(".pdf", pdf.file)) pdf.file <- paste(pdf.file, ".pdf", sep="")

  dots <- list(...)  # check for deprecated/changed parameters
  if (length(dots) > 0) {
    for (i in 1:length(dots)) {
      old.nm <- c("col.fill", "col.stroke", "col.bg", "col.grid", "col.box",
                  "col.reg", "col.axis", "col.trans", "col.low", "col.hi")
      if (names(dots)[i] %in% old.nm) {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
          "options that began with the abbreviation  col  now begin with  ",
          "color \n\n")
        }
      if (names(dots)[i] == "addtop") 
        cat("\naddtop  is now a multiplicative factor instead of additive\n\n")
      if (names(dots)[i] == "count.levels") {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
          "Now use  count.labels  instead of count.levels\n\n")
      }
    }
  }

  x.name <- deparse(substitute(x))
  options(xname = x.name)

  df.name <- deparse(substitute(data))
  options(dname = df.name)


# -----------------------------------------------------------
# establish if a data frame, if not then identify variable(s)

  x.call <- NULL

  if (!missing(x)) {
    if (!exists(x.name, where=.GlobalEnv)) {  # x not in global env, in df
      .nodf(df.name)  # check to see if data frame container exists 
      .xcheck(x.name, df.name, data)  # see if var in df, vars lists not checked
      vars.list <- as.list(seq_along(data))
      names(vars.list) <- names(data)
      x.col <- eval(substitute(x), envir=vars.list)  # col num of each var
      if (length(x.col) > 1) data <- data[, x.col]  # x is a vars list
      if (length(x.col) == 1) x.call <- eval(substitute(data$x))  # x is 1 var
    }
    else {  # x is in the global environment (vector, matrix or data frame)
      if (is.data.frame(x))  # x a data frame
        data <- x
      else {  # x a vector or matrix in global
        if (exists(x.name, where=.GlobalEnv)) if (is.matrix(x)) { 
          x.name <- xlab
          xlab <- NULL
          y.name <- legend.title
          options(xname = x.name)
          options(yname = y.name)
        }
        x.call <- x
        if (is.function(x.call)) x.call <- eval(substitute(data$x))
      }
    }
  }


  if (!is.null(x.call)) {  # x is a single var, not a data frame

    # evaluate by
    if (!missing(by)) {

      # get actual variable name before potential call of data$x
      y.name <- deparse(substitute(by)) 
      options(yname = y.name)

      # see if y exists from a function call
      # indicate a function call with sys.nframe returns larger than 1 
      #if (exists(y.name, where=parent.frame(n=1)) && sys.nframe() > 1) 
        #in.call <- TRUE else in.call <- FALSE

      # get conditions and check for data existing
      #if (!in.call) {
        xs <- .xstatus(y.name, df.name, quiet)
        in.global <- xs$ig 
      #}
      #else in.global <- FALSE
      # if y is in global, sys.nframe() returns two, in.call is TRUE,
      #   which leads to in.global FALSE
      #if (exists(y.name, where=.GlobalEnv)) in.global <- TRUE

      # see if var exists in data frame, if x not in Global Env or function call 
      if (!in.global) .xcheck(y.name, df.name, data)
      #if (!in.global && !in.call) .xcheck(y.name, df.name, data)
      if (!in.global)
        y.call <- eval(substitute(data$by))
      else {  # vars that are function names get assigned to global
        y.call <- by
        if (is.function(y.call)) y.call <- eval(substitute(data$by))
      }

    }
    else
      y.call <- NULL


    # evaluate count.labels
    #---------------------
    if (!missing(count.labels)) {

      # get actual variable name before potential call of data$x
      x.name <- deparse(substitute(count.labels)) 
      options(xname = x.name)

      # get conditions and check for data existing
      xs <- .xstatus(x.name, df.name, quiet)
      in.global <- xs$ig 

      # see if var exists in data frame, if x not in Global Env or function call 
      if (!missing(x) && !in.global)
        .xcheck(x.name, df.name, data)

      if (!in.global) count.labels.call <- eval(substitute(data$count.labels))
      else {  # vars that are function names get assigned to global
        count.labels.call <- count.labels
        if (is.function(count.labels.call)) 
          count.labels.call <- eval(substitute(data$count.labels))
      }
 
      #if (!.is.integer(x.call)) { 
      #cat("\n"); stop(call.=FALSE, "\n","------\n",
        #"Values to be analyzed must be counts, that is, integers\n",
        #"First two values to analyze: ", x.call[1], " ", x.call[2], "\n\n")
      #}
    }
    else
      count.labels.call <- NULL


    if (length(unique(na.omit(x.call))) == 1) { 
      cat("\n"); stop(call.=FALSE, "\n","------\n",
        "There is only one unique value for the values of ", x.name,
        ": ", na.omit(x.call)[1], "\n",
        "The bar chart is only computed if there is more than one",
        " unique value\n\n")
    }

# ---------------
# do the analysis

    .opendev(pdf.file, pdf.width, pdf.height)

    orig.params <- par(no.readonly=TRUE)
    on.exit(par(orig.params))

    bc <- .bc.main(x.call, y.call,
         color.fill, color.stroke, color.bg, color.grid, color.box, colors,
         horiz, over.grid, addtop, gap, proportion, xlab, ylab, main, value.labels,
         cex.axis, color.axis, rotate.values, offset, beside, color.low, color.hi,
         count.labels.call,
         legend.title, legend.loc, legend.labels, legend.horiz, quiet, ...)

    if (!is.null(pdf.file)) {
      dev.off()
      if (!quiet) .showfile(pdf.file, "barchart")
    }

    invisible(bc)
  }
  

  else
    bc.data.frame(data, n.cat,
      color.fill, color.stroke, color.bg, color.grid, color.box, colors,
      horiz, over.grid, addtop, gap, proportion, xlab, ylab, main, value.labels,
      cex.axis, color.axis, rotate.values, offset, beside, color.low, color.hi,
      count.labels,
      legend.title, legend.loc, legend.labels, legend.horiz, quiet,
      pdf.width, pdf.height, ...)

}

