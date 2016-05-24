Plot <- 
function(x, y=NULL, by=NULL, data=mydata, n.cat=getOption("n.cat"),

         topic=c("data", "count", "prop", "mean", "sd", "min", "median", "max",
                 "diff"),
         object=c("point", "line", "both", "bubble", "sunflower", "bar", "off"),

         color.fill=getOption("color.fill.pt"),
         color.stroke=getOption("color.stroke.pt"),
         color.bg=getOption("color.bg"),
         color.grid=getOption("color.grid"),
         color.box=getOption("color.box"),

         color=NULL, color.trans=NULL, color.area=NULL,

         cex.axis=0.76, color.axis="gray30", xy.ticks=TRUE,
         xlab=NULL, ylab=NULL, main=NULL, sub=NULL,
         value.labels=NULL, rotate.values=0, offset=0.5,

         size=NULL, shape="circle", means=TRUE, 
         sort.yx=FALSE, segments.y=FALSE, segments.x=FALSE,

         bubble.size=0.25, bubble.power=0.6, bubble.counts=TRUE,
         color.low=NULL, color.hi=NULL,

         fit.line=NULL, color.fit.line="gray55",

         ellipse=FALSE, color.ellipse="lightslategray",
         color.fill.ellipse="off", 

         method="overplot", pt.reg="circle", pt.out="circle", 
         color.out30="firebrick2", color.out15="firebrick4", new=TRUE,

         breaks="Sturges", bin.start=NULL, bin.width=NULL, bin.end=NULL,
         prop=FALSE, cumul=c("off", "on", "both"), hist.counts=FALSE,
         color.reg="snow2",

         beside=FALSE, horiz=FALSE, 
         over.grid=FALSE, addtop=0.05, gap=NULL, count.labels=NULL,
         legend.title=NULL, legend.loc="right.margin", legend.labels=NULL,
         legend.horiz=FALSE, 

         digits.d=NULL, quiet=getOption("quiet"),
         pdf.file=NULL, pdf.width=NULL, pdf.height=NULL,
         fun.call=NULL, ...) {


  if (is.null(fun.call)) fun.call <- match.call()

  object.miss <- ifelse (missing(object), TRUE, FALSE)
  topic.miss <- ifelse (missing(topic), TRUE, FALSE)

  object <- match.arg(object)
  topic <- match.arg(topic)
  cumul <- match.arg(cumul)

  if (object.miss) object <- "default"

  # any bubble parameter actives a bubble plot
  if (!missing(bubble.size) || !missing(bubble.size) || !missing(bubble.counts))
    object <- "bubble"

  # any ellipse parameter actives an ellipse
  if (!missing(color.ellipse) || !missing(color.fill.ellipse))
    ellipse <- TRUE

  # any bin parameter activates bins
  if (!missing(breaks)  ||  !missing(bin.start)  ||  !missing(bin.width)  ||
      !missing(bin.end))
    object <- "bar"

  if (object == "bar") {
    if (topic.miss) topic <- "count"  # default topic for bars
    if (topic == "data") {
      cat("\n"); stop(call.=FALSE, "\n","------\n",
      "Bars do not apply to data, only to count, means, etc\n\n")
    }
    if (missing(color.fill)) color.fill <- getOption("color.fill.bar")
  }

  for (i in 1:length(color.fill))
    if (color.fill[i] == "off") color.fill[i] <- "transparent"
  for (i in 1:length(color.stroke))
    if (color.stroke[i] == "off") color.stroke[i] <- "transparent"
  if (color.bg == "off") color.bg <- "transparent"
  if (color.grid == "off" ) color.grid <- "transparent"
  if (color.box == "off") color.box <- "transparent"
  if (color.fill.ellipse == "off") color.fill.ellipse <- "transparent"

  if (topic == "diff"  &&  missing(y)) {
    cat("\n"); stop(call.=FALSE, "\n","------\n",
        "Need both an x-variable and a y-variable to compute \n",
        "  their difference\n\n")
  }

  if (!is.null(color)) {
    color.stroke <- color
    color.fill <- color
  }

  if (is.null(fit.line)) fit.ln <- "none"
  if (is.logical(fit.line))
    fit.ln <- ifelse (fit.line, "loess", "none")
  if (is.character(fit.line)) {
    if (!(fit.line %in% c("loess", "ls", "none"))) { 
      cat("\n"); stop(call.=FALSE, "\n","------\n",
        "fit.line applies only for  loess  or  ls  (least squares)\n\n")
    }
    fit.ln <- fit.line
  }

  if (is.logical(ellipse)) if (ellipse) ellipse <- 0.95
  if (as.logical(ellipse[1])) {
    txt <- "[Ellipse with Murdoch and Chow's function ellipse"
    cat(txt, "from the ellipse package]\n") 
  }

  if (!is.null(pdf.file))
    if (!grepl(".pdf", pdf.file)) pdf.file <- paste(pdf.file, ".pdf", sep="")


  dots <- list(...)  # check for deprecated parameters
  if (length(dots) > 0) {
    for (i in 1:length(dots)) {
      old.nm <- c("col.fill", "col.stroke", "col.bg", "col.grid", "col.box",
                  "col.reg", "col.axis", "col.trans", "col.low", "col.hi",
                  "col.ellipse", "col.fill.ellipse", "col.fit.line", "col.out30",
                  "col.out15") 
      if (names(dots)[i] %in% old.nm) {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
          "options that began with the abbreviation  col  now begin with  ",
          "color \n\n")
      }
      if (names(dots)[i] %in% c("x.start","x.end","y.start","y.end")) {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
          "x.start, x.end, y.start, and y.end no longer used\n\n",
          "Instead use the standard R xlim and ylim parameters,\n",
          "such as xlim=c(0,40) to specify from 0 to 40. Same for ylim.\n\n")
      }
      if (names(dots)[i] == "kind") {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
          "option  kind  is renamed  object\n\n")
      }
      if (names(dots)[i] == "knitr.file") {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
          "knitr.file  no longer used\n",
          "Instead use  Rmd  for R Markdown file\n\n")
      }
      if (names(dots)[i] == "type") {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
          "type  option replaced with  object\n\n")
      }
      if (names(dots)[i] == "diag") {
        cat("\n"); stop(call.=FALSE, "\n","------\n",
          "diag  option no longer used\n\n")
      }
    }
  }

  if (method == "stack") { 
    cat("\n"); stop(call.=FALSE, "\n","------\n",
      "Option  stack  does not work\n\n")
  }

  if (missing(x)) { 
    cat("\n"); stop(call.=FALSE, "\n","------\n",
      "Must specify at least one variable to analyze\n\n")
  }

  if (topic %in% c("mean", "sd", "min", "max") && missing(y)) {
      cat("\n"); stop(call.=FALSE, "\n","------\n",
      "Must specify a numeric y-variable from which to compute the\n ",
      " ", topic, " for each level of ", deparse(substitute(x)), "\n\n")
  }

  # conflict between graphical and cor parameters, so avoid
  if (method %in% c("spearman", "kendall")) {
    cat("\n"); stop(call.=FALSE, "\n","------\n",
      "The  method  parameter has another meaning for ScatterPlot\n\n",
      "Compute a Spearman or Kendall correlation\n",
      "  with the Correlation function\n\n")
  }
        
  if (is.numeric(breaks) && !is.null(bin.start)) { 
    cat("\n"); stop(call.=FALSE, "\n","------\n",
      "Choose only one option to specify a start value.\n",
      "Either choose the option  breaks  or the option  bin.start.\n\n")
  }
        
  if (topic != "data"  &&  object == "sunflower") { 
    cat("\n"); stop(call.=FALSE, "\n","------\n",
      "Sunflowers are only plotted for data\n\n")
  }


  # process shapes
  bad.shape <- NULL
  shapes <- c("circle", "square", "diamond", "triup", "tridown")
  shapes.all <- c(shapes, c(21:25), letters, LETTERS, 0:9, "+", "*", "#",
                  "%", "!", "=", "-", "&", "$", "?", "|", ">", "<", "@")

  num.flag <- FALSE
  for (i in 1:length(shape)) {
    if (!(shape[i] %in% shapes.all)) 
      bad.shape <- shape[i]
    else
      if (shape[i] %in% shapes) {
        shape[i] <- which(shape[i]==shapes)+20
        num.flag <- TRUE
      }
  }
  if (num.flag) shape <- as.numeric(shape)

  if (pt.reg %in% shapes) 
    pt.reg <- which(pt.reg==shapes) + 20
  else
    if (!(pt.reg %in% c(21:25))) bad.shape <- pt.reg

  if (pt.out %in% shapes) 
    pt.out <- which(pt.out==shapes) + 20
  else
    if (!(pt.out %in% c(21:25))) bad.shape <- pt.out

  if (!is.null(bad.shape)) {
      message("\nValid shapes") 
      message("------------") 
      for (j in 1:length(shapes)) message(shapes[j])
      message("all uppercase and lowercase letters")
      message("all digits")
      message("+ * # % ! = - & $ ? | < > @")
      cat("\n")
      stop(call.=FALSE, "\n","------\n",
      "Not a valid shape: ", bad.shape, "\n\n")
  }


  # ------
  # get actual variable name before potential call of data$x
  x.name <- deparse(substitute(x)) 
  options(xname = x.name)

  # get data frame name
  df.name <- deparse(substitute(data))
  options(dname = df.name)

  x.only <- FALSE  # flag for two-variable Cleveland dot plot 

  if (deparse(substitute(x)) == "row.names") {
    data.x <- data.frame(factor(row.names(data)))
    if (is.null(xlab)) xlab <- ""  # unless specified, drop the axis label
  }

  else {
    if (!exists(x.name, where=.GlobalEnv)) {  # x not in global env, in df
      .nodf(df.name)  # check to see if data frame container exists 
      .xcheck(x.name, df.name, data)  # var in df?, vars lists not checked
      all.vars <- as.list(seq_along(data))  # even if only a single var
      names(all.vars) <- names(data)  # all data in data frame
      x.col <- eval(substitute(x), envir=all.vars)  # col num selected vars
      if (!("list" %in% class(data))) {
        data.x <- data[, x.col]
        if (length(x.col) == 1) {  # x is 1 var
          data.x <- data.frame(data.x)
          names(data.x) <- x.name
        }
      }
      else {  # class of data is "list"
        data.x <- data.frame(data[[x.col]])
        names(data.x) <- x.name
      }
    }

    else { # x is in the global environment (vector or data frame)
      if (is.data.frame(x))  # x a data frame
        data.x <- x
      else {  # x a vector in global
        if (!is.function(x)) {
          data.x <- data.frame(x)  # x is 1 var
          x.only <- TRUE  # to flag later for analyzing differences
        }
        else 
          data.x <- data.frame(eval(substitute(data$x)))  # x is 1 var
        names(data.x) <- x.name
      }
    }
  }

  # just one x variable for now, a vector of cat or num values
  if (ncol(data.x) == 1) { 
    x.call <- data.x[,1]

    nu <- length(unique(na.omit(x.call)))
    num.cat.x <- .is.num.cat(x.call, n.cat)

    if (!num.cat.x && is.integer(x.call) && nu <= n.cat) {
      cat("\n")
      cat(">>> ", x.name, " has only only ", nu, " unique ",
          "integer values, but not equally spaced,\n",
          "      so treat as numerical in this analysis\n",
          "    Maybe convert to an R factor to treat as categorical\n",
          sep="")
    }

    if (num.cat.x && !quiet) .ncat("ScatterPlot", x.name, nu, n.cat)
    cat.x <- ifelse (num.cat.x || is.factor(x.call), TRUE, FALSE)

    if (cat.x  &&  object %in% c("line", "both")) {
      cat("\n"); stop(call.=FALSE, "\n","------\n",
        "A line applies only to continuous variables\n",
        x.name, " is a categorical variable\n\n")
    }
  }

  else {  # more than one x-variable
    if (missing(y)  &&  object == "point") { # no y, so all selected vars are categorical
      for (i in 1:length(x.col)) {
        nu <- length(unique(na.omit(data.x[,i])))
        num.cat <- .is.num.cat(data.x[,i], n.cat)
        cat.all <- ifelse (num.cat || is.factor(data.x[,i]), TRUE, FALSE)
        if (!cat.all && object == "point") { 
          cat("\n"); stop(call.=FALSE, "\n","------\n",
            "All variables for multiple x-variables with no y-variable\n",
            "  must be categorical for a Bubble Plot Frequency Matrix\n\n",
            "A categorical variable is either an R factor variable,\n",
            "  or is numeric with not more than  n.cat  unique values\n\n",
            "Can set  n.cat  locally when calling  ScatterPlot,\n",
            "  or globally with the function:  theme\n\n")
        }
      }
    }
    else {  # y is present, all x must be numerical

      if (as.logical(ellipse[1])) { 
        cat("\n"); stop(call.=FALSE, "\n","------\n",
          "An ellipse only applies to an analysis of a single x-variable\n",
          "  with a single y-variable\n\n")
      }

      # everything goes to numeric, NA if character, internal code if factor
      x.call <- data.matrix(data.x, rownames.force=FALSE)
      cat.x <- FALSE
    }
  }  # end more than 1 x-variable


  # evaluate y
  #-----------
  if (!missing(y)) {

    # get actual variable name before potential call of data$x
    y.name <- deparse(substitute(y)) 
    options(yname = y.name)

    if (deparse(substitute(y)) == "row.names") {
       y.call <- factor(row.names(data))
       if (is.null(ylab)) ylab <- ""  # unless specified, drop the axis label
    }

    else {
      if (!exists(y.name, where=.GlobalEnv)) {  # y not in global env, in df
        .nodf(df.name)  # check to see if data frame container exists 
        .xcheck(y.name, df.name, data)  # var in df?, vars lists not checked
        all.vars <- as.list(seq_along(data))  # even if only a single var
        names(all.vars) <- names(data)  # all data in data frame
        y.col <- eval(substitute(y), envir=all.vars)  # col num selected vars
        if (!("list" %in% class(data))) {
          data.y <- data[, y.col]
          if (length(y.col) == 1) {  # y is 1 var
            data.y <- data.frame(data.y)
            names(data.y) <- y.name
          }
        }
        else {  # class of data is "list"
          data.y <- data.frame(data[[y.col]])
          names(data.y) <- y.name
        }
      }

      else { # y is in the global environment (vector or data frame)
        if (is.data.frame(y))  # y a data frame
          data.y <- y
        else {  # y a vector in global
          if (!is.function(y))
            data.y <- data.frame(y)  # y is 1 var
          else
            data.y <- data.frame(eval(substitute(data$y)))  # y is 1 var
          names(data.y) <- y.name
        }
      }

      if (ncol(data.y) == 1) { 
        y.call <- data.y[,1]

        nu <- length(unique(na.omit(y.call)))
        num.cat.y <- .is.num.cat(y.call, n.cat)

        if (!num.cat.y && is.integer(y.call) && nu <= n.cat) {
          cat("\n")
          cat(">>> ", y.name, " has only only ", nu, " unique ",
              "integer values, but not equally spaced,\n",
              "      so treat as numerical in this analysis\n",
              "    Maybe convert to an R factor to treat as categorical\n",
              sep="")
        }

        if (num.cat.y && !quiet) .ncat("ScatterPlot", y.name, nu, n.cat)
        cat.y <- ifelse (num.cat.y || is.factor(y.call), TRUE, FALSE)
      }
      else {
        cat.y <- FALSE  #  multiple y-vars must be numerical
        if (missing(ylab)) ylab <- ""  # use legend instead
        y.call <- data.matrix(data.y, rownames.force=FALSE)
      }
    }
  }
  else
    y.call <- NULL


  # evaluate by
  #-----------
  if (!missing(by)) {

    # get actual variable name before potential call of data$x
    by.name <- deparse(substitute(by)) 
    options(byname = by.name)

    # get conditions and check for data existing
    xs <- .xstatus(by.name, df.name, quiet)
    in.global <- xs$ig 

    # see if var exists in data frame, if x not in Global Env or function call 
    if (!missing(x) && !in.global)
      .xcheck(by.name, df.name, data)

    if (!in.global)
      by.call <- eval(substitute(data$by))
    else {  # vars that are function names get assigned to global
      by.call <- by
      if (is.function(by.call)) by.call <- eval(substitute(data$by))
    }

    if (!is.factor(by.call)) by.call <- factor(by.call)
  }
  else
   by.call <- NULL 

  # graphics 
  # --------
  if (is.null(y.call)  &&  ncol(data.x) == 1  &&  method!="stack"
      &&  !(object %in% c("line", "both"))  && object != "bar")
    plt.h <- ifelse(is.null(main), 2.5, 3.1)  # narrow for 1-D plot
  else
    plt.h <- 4.5
  # for BPFM with more than 7 variables, make extra long
  if (ncol(data.x) > 7) plt.h <- plt.h + ((ncol(data.x) - 7) * 0.5)

  if (is.null(pdf.file)) {  
    if (options("device") != "RStudioGD"  &&  is.null(options()$knitr.in.progress)) {
      if (missing(by)) {  # set up graphics system to manage
        if (is.null(y.call)) {  # works for 1-D scatter plots and BPFM
          .graphwin(1, d.h=plt.h) 
        }
        else
          .graphwin(1) 
      }
      else
        .graphwin(d.w=pdf.width)  # add width to default of 4.5 for legend

      }
    }
  else  {  # pdf file
    if (is.null(pdf.width)) pdf.width <- 4.5
    if (!missing(by.call)) pdf.width <- pdf.width + 0.6
    if (is.null(pdf.height)) pdf.height <- plt.h
    pdf(file=pdf.file, width=pdf.width, height=pdf.height)
  }


  # ------------------------------------------------
  # set object and topic where needed

  if (is.null(y.call) && cat.x  &&  topic == "data") {
    y.call <- rep(0, length(x.call))
    if (object == "default") object <- "bubble"
  }

  # if numeric x is sorted with equal intervals, line chart
  if (is.numeric(x.call)) {
    if (object == "default") {
      if (sum(is.na(x.call)) > 0)
        eq.int <- FALSE  # missing data in x present
      else {
        eq.int <- TRUE
        if (is.matrix(x.call))
          d.x <- diff(x.call[,1])  # only look at first x-variable
        else
          d.x <- diff(x.call)  # only look at first x-variable
        for (i in 2:(length(d.x))) 
          if ((abs(d.x[i-1] - d.x[i]) > 0.0000000001)) eq.int <- FALSE
        rm(d.x)
      }  # also no y missing

      if(!is.unsorted(x.call) && eq.int && sum(is.na(y))==0) object <- "line"
    }
  }

  # if numeric y is sorted with equal intervals, line chart
  if (is.numeric(y.call)) {
    if (object == "default") {
      if (sum(is.na(y.call)) > 0)
        eq.int <- FALSE  # missing data in x present
      else {
        eq.int <- TRUE
        if (is.matrix(y.call))
          d.y <- diff(y.call[,1])  # only look at first x-variable
        else
          d.y <- diff(y.call)  # only look at first x-variable
        for (i in 2:(length(d.y))) 
          if ((abs(d.y[i-1] - d.y[i]) > 0.0000000001)) eq.int <- FALSE
        rm(d.y)
      }  # also no y missing

      if(!is.unsorted(y.call) && eq.int && sum(is.na(x))==0) object <- "line"
    }
  }

  # set bubble plot for small number of unique values of x and y
  if (object == "default") {  # set default
    if (!missing(y)) {
      object <- "point"
      if (topic == "data") if (cat.x && cat.y) object <- "bubble"
    }
    else if (topic %in% c("count", "prop"))
      object <- "point"
    else {
      if (topic == "data") if (cat.x) object <- "bubble"
      if (is.matrix(x.call)  &&  missing(y)) object <- "bubble"
    }
  }

  if (!quiet) {
    cat("\n")
    cat(">>> geometric object to plot: object = \"", object, "\"\n", sep="")
    cat(">>> subject of the analysis:  topic = \"", topic, "\"\n", sep="")
    cat("\n")
  }


  # ------------------------------------------------
  # analysis

  if (is.data.frame(data.x) == "data.frame") {  # correlation matrix
    .cr.data.frame(x, miss="pairwise", show.n=FALSE, n.cat, digits.d=2,
                   heat.map=FALSE, colors=NULL,
                   main=NULL, bottom=NULL, right=NULL, ...)
  }

  else if (object == "point"  &&  topic == "diff") {

    # construct x.call
    if (!missing(y)) {
      x.call <- data.frame(data.x, data.y)
      names(x.call) <- c(x.name, y.name)
    }
    else {
      x.call <- data.frame(data.x)
      names(x.call) <- x.name
    }
    x.call <- data.matrix(x.call, rownames.force=FALSE)

    # construct y.call, need unique values for Cleveland dot plot
    if (x.only) data <- data.frame(x.call)  # input data are vectors
    is.unique <- logical(ncol(data))  # initialed to FALSE
    for(i in 1:ncol(data))
      if ((length(data[,i])==length(unique(data[,i]))) && !is.numeric(data[,i]))
        is.unique[i] <- TRUE
    unq <- which(is.unique)[1]  # choose first non-num variable with unique values
    if (!is.na(unq))
      y.call <- factor(data[, unq[1]])
    else { # no non-num var with unique values, so go to row.names
      if (row.names(data)[1] != "1")
        y.call <- factor(row.names(data))
      else
        y.call <- factor(row.names(data), levels=1:nrow(data))  # proper sort
    }

    # Cleveland two-variable dot plot
    .plt.main(x.call, y.call, by=NULL, n.cat=getOption("n.cat"),
       col.area=NULL, col.box="black", col.grid="transparent",
       col.trans=getOption("col.trans.pt"), shape=21, size=.8, ylab="",
       segments.y=TRUE, quiet=TRUE)

  }

  else if (object == "bar"  &&  topic == "count") {

    if (!cat.x) {  # histogram

      h <- .hst.main(x.call, color.fill, color.stroke, color.bg, color.grid,
         color.box, color.reg,
         over.grid=FALSE, cex.axis, color.axis, rotate.values, offset,
         breaks, bin.start, bin.width, bin.end, prop, hist.counts, cumul,
         xlab, ylab, main, sub, quiet, fun.call=NULL, ...)

      if (!quiet) {

        stats <- .hst.stats(h, length(x.call), fun.call)

        txsug <- stats$txsug
        txdst <- stats$tx

        bin.width <- stats$bin.width
        n.bins <- stats$n.bins
        prop <- stats$prop
        cum.c <- stats$counts_cum
        cum.p <- stats$prop_cum

        txotl <- .outliers(data[,i])
        if (length(txotl)==0) txotl <- "No outliers"

        if (ncol(data) > 1) {  # for a variable range, print the text output
          class(txsug) <- "out_piece"
          class(txdst) <- "out_piece"
          class(txotl) <- "out_piece"
          output <- list(out_suggest=txsug, out_freq=txdst, out_outliers=txotl)
          class(output) <- "out_all"
          print(output)
        }
      }
    }

    else {  # bar chart
      .bc.main(x.call, by=y.call,
        color.fill, color.stroke, color.bg, color.grid, color.box,
        colors=getOption("colors"),
        horiz, over.grid, addtop, gap,
        prop, xlab, ylab, main, value.labels,
        cex.axis, color.axis, rotate.values, offset, beside,
        color.low, color.hi, count.labels,
        legend.title, legend.loc, legend.labels, legend.horiz,
        quiet=quiet, ...)
    }
  }
 
  # bubble plot frequency matrix
  else if (ncol(data.x) > 1 && missing(y) && object == "bubble") { 

    # get labels just for subset data matrix
    mylabels <- attr(data, which="variable.labels")
    nm <- names(data.x)
    if (!is.null(mylabels)) {
      mylabs <- character(length=0) 
      for (i in 1:length(nm))
        mylabs[i] <- mylabels[which(names(mylabels) == nm[i])]
    }
    else
      mylabs <- NULL
    if (length(mylabs) == 0) mylabs <- NULL  # when labels, but not relevant

    if (is.null(xlab)) xlab <- ""  # suppress x-axis label if not specified
print(bubble.power)

    .dpmat.main(data[,x.col], mylabs, nm,
      color.fill, color.stroke, color.bg, color.grid, color.trans,
      shape, color.area, color.box, 
      cex.axis, color.axis, color.low, color.hi,
      xy.ticks, xlab, ylab, main, sub, size,
      bubble.size, bubble.counts, bubble.power,
      value.labels, rotate.values, offset, quiet, ...)
  }

  else {

    if (topic != "data") {  # do stats output before reducing data

      n.cat <- 0
      means <- FALSE
      if (missing(size)) size <- 1.3
      if (missing(segments.x)) segments.x <- TRUE
      if (is.null(color.trans)) color.trans <- 0

      if (!quiet) {
        if (!missing(y)) {
          if (cat.y) {
            cat("\n"); stop(call.=FALSE, "\n","------\n",
            y.name, " is not numerical, so cannot compute its mean\n\n")
          }
          options(yname = x.name)  # reverse order of x and y for .ss.numeric
          options(xname = y.name)
          stats <- .ss.numeric(y.call, by=x.call, digits.d=digits.d, brief=TRUE)
          txout <- stats$tx
          options(xname = x.name)  # reverse back
          options(yname = y.name)
        }
        else  {
          stats <- .ss.factor(x.call, digits.d=digits.d, x.name=x.name, brief=TRUE)
          txout <- stats$counts
        }

        class(txout) <- "out_piece"

        output <- list(out_txt=txout)
        class(output) <- "out_all"
        print(output)
      }
    }
 
    # only 1 variable, so get y.call to plot points
    if (is.null(y.call) && cat.x) {
      if (topic %in% c("count", "prop")) {
        ylab <- ifelse(topic=="count", "Count of", "Proportion of")
        ylab <- paste(ylab, x.name)
        frq <- table(x.call)
        if (topic == "prop") frq <- frq / sum(frq)
        if (is.factor(x.call))  # preserve ordering, will lose order attribute
          x.call <- factor(names(frq), levels=levels(x.call))
        else
          x.call <- factor(names(frq))
        y.call <- as.vector(frq)
      }
    }

    # line chart
    if (is.null(y.call) &&  !cat.x  && object %in% c("line", "both")
        && topic == "data") {
      y.call <- x.call
      options(yname = x.name)
      options(xname = "Index")
      #xlab <- "Index"
      if (!is.matrix(x.call)) 
        x.call <- 1:length(x.call) 
      else {
        x.call <- 1:nrow(x.call)
        #ylab <- ""
      }
    }

    # set up new x.call and y.call for stats
    if (topic %in% c("mean", "sd", "min", "median", "max")) {
      if (topic == "mean") {
        ylab <- paste("Mean of", y.name) 
        out <- tapply(y.call, x.call, mean, na.rm=TRUE)
      }
      if (topic == "sd") {
        ylab <- paste("Standard Deviation of", y.name)
        out <- tapply(y.call, x.call, sd, na.rm=TRUE)
      }
      if (topic == "min") {
        ylab <- paste("Minimum of", y.name)
        out <- tapply(y.call, x.call, min, na.rm=TRUE)
      }
      if (topic == "median") {
        ylab <- paste("Median of", y.name)
        out <- tapply(y.call, x.call, median, na.rm=TRUE)
      }
      if (topic == "max") {
        ylab <- paste("Maximum of", y.name)
        out <- tapply(y.call, x.call, max, na.rm=TRUE)
      }

      if (is.factor(x.call))  # preserve ordering, will lose order attribute
        x.call <- factor(names(out), levels=levels(x.call))
      else {
        if (is.numeric(x.call)) {
          m1 <- min(sort(unique(x.call)))
          m2 <- max(sort(unique(x.call)))
          x.call <- factor(names(out), levels=m1:m2)  # get entire numeric range
        }
        else
          x.call <- factor(names(out))
      }
      y.call <- as.vector(out)

      if (object == "bar") {
        names(y.call) <- x.call
        .bc.main(y.call, by=NULL,
          color.fill, color.stroke, color.bg, color.grid, color.box,
          colors=getOption("colors"),
          horiz, over.grid, addtop, gap,
          prop, xlab, ylab, main, value.labels,
          cex.axis, color.axis, rotate.values, offset, beside,
          color.low, color.hi, count.labels,
          legend.title, legend.loc, legend.labels, legend.horiz,
          quiet=quiet, ...)
      }
    }  # mean, sd, min, median, max

    # 2-variable scatter plot
    if (!is.null(y.call)  &&  object != "bar") {  
      .plt.main(x.call, y.call, by.call, n.cat,
         object, topic,
         color.fill, color.stroke, color.bg, color.grid, color.box,
         color.trans, color.area,
         cex.axis, color.axis, xy.ticks,
         xlab, ylab, main, sub, value.labels, rotate.values, offset,
         size, shape, means, sort.yx, segments.y, segments.x,
         bubble.size, bubble.power, bubble.counts, color.low, color.hi,
         fit.ln, color.fit.line,
         ellipse, color.ellipse, color.fill.ellipse, 
         method, pt.reg, pt.out, color.out30, color.out15, new,
         quiet, fun.call, ...)
    }

    # 1-D traditional scatter plot (not bubble plot)
    else if (!cat.x  &&  is.null(y.call)) {  
       
        .dp.main(x.call, by.call, size,
           color.fill, color.stroke, color.bg, color.grid, color.trans,
           shape, cex.axis, color.axis, xlab, main, sub, size, 
           rotate.values, offset, method, pt.reg, pt.out, 
           color.out30, color.out15, quiet, new, ...)

        # R Markdown
        #txkfl <- ""
        #if (!is.null(Rmd)) {
          #if (!grepl(".Rmd", Rmd)) Rmd <- paste(Rmd, ".Rmd", sep="")
          #txknt <- .dist.Rmd(x.name, df.name, fun.call, digits.d)
          #cat(txknt, file=Rmd, sep="\n")
          #txkfl <- .showfile2(Rmd, "R Markdown instructions")
        #}

        #class(txkfl) <- "out_piece"

        #output <- list(
          #call=fun.call,
          #type="1D_ScatterPlot", out_file=txkfl)

        #class(output) <- "out_all"

        #return(output)
 
    }  # end 1-D plot

  }


  # terminate pdf graphics system if used
  if (!is.null(pdf.file)) {
    dev.off()
    .showfile(pdf.file, "plot")
  }

}
