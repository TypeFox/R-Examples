.pc.main <- 
function(x,
         random.col, col.fill, col.low, col.hi,
         colors, cex, cex.main, quiet, main, 
         pdf.file, pdf.width, pdf.height, ...)  {

  # get values for ... parameter values
  stuff <- .getdots(...)
  col.main <- stuff$col.main

  # set the labels
  # use variable label for main if it exists and main not specified
  gl <- .getlabels(main=main)
  x.name <- gl$xn; x.lbl <- gl$xl
  cex.lab <- gl$cex.lab

  if (!is.null(main)) main.lbl <- main
  else if (length(x.lbl) == 0) main.lbl <- x.name else main.lbl <- x.lbl

  cex.lab <- 0.85
  if (strwidth(main.lbl, units="figure", cex=cex.lab) > .85) {
    brk <- nchar(main.lbl)
    while (strwidth(substr(main.lbl,1,brk), units="figure", cex=cex.lab) > .85)
      brk <- brk-1 
    while (substr(main.lbl,brk,brk) != " ") brk <- brk-1
    main.lbl <- paste(substr(main.lbl,1,brk), "\n",
                      substr(main.lbl,brk+1,nchar(main.lbl)))
    while (strwidth(main.lbl, units="figure", cex=cex.lab) > .85)
      cex.lab <- cex.lab-0.05
  }

  # entered counts typically integers as entered but stored as type double
  # if names(x) is null, likely data from sample and c functions
  if (!is.integer(x) && is.double(x) && !is.null(names(x)))  x <- as.table(x)
  if (!is.factor(x) && !is.table(x)) x <- factor(x)
  ncolors <- ifelse (!is.table(x), nlevels(x), length(x))

  # color palette
  # set some default colors in case not assigned below
  if (is.ordered(x)) {
    lowhi <- .ordcolors(colors, col.low, col.hi) 
    col.low <- lowhi$col.low
    col.hi <- lowhi$col.hi

    color.palette <- colorRampPalette(c(col.low, col.hi))
    clr <- color.palette(ncolors)
  }  # is ordered

  else if (!is.null(col.low) && !is.null(col.hi)) {
      color.palette <- colorRampPalette(c(col.low, col.hi))
      clr <- color.palette(ncolors)
  }
  else {
    if (colors == "gray") {
      color.palette <- colorRampPalette(c("gray28","gray92"))
      clr <- color.palette(ncolors)
    }
    else if (colors == "rainbow") clr <- rainbow(ncolors)
    else if (colors == "terrain") clr <- terrain.colors(ncolors)
    else if (colors == "heat") clr <- heat.colors(ncolors)
    else  {  # ordered color range does not make sense here 
      clr <- .col.discrete()
    }
    if (random.col) clr <- clr[sample(length(clr))]
  }

  if (!is.null(col.fill)) {
    for (i in 1:(min(length(col.fill),length(clr)))) clr[i] <- col.fill[i]
    ncolors <- min(length(col.fill),length(clr))
  }

  palette(clr)
  col <- 1:ncolors 

  # plot the pie chart
  if (!is.table(x)) x <- table(x)
  pie(x, col=col, main=main.lbl, cex=cex, cex.main=cex.main, col.main=col.main, ...)

# legend("bottom", legend=unique(na.omit(x)), horiz=TRUE, cex=0.8, fill=col)

  # text output
  if (length(dim(x)) == 1  && !quiet) {  # one variable

    stats <- .ss.factor(x, brief=TRUE, x.name=x.name)

    txttl <- stats$title
    counts <- stats$counts
    chi <- stats$chi
    class(txttl) <- "out_piece"
    class(counts) <- "out_piece"
    class(chi) <- "out_piece"
    output <- list(out_title=txttl, out_counts=counts, out_chi=chi)
    class(output) <- "out_all"
    print(output)      
  }

 
  cat("\n")

}  #  end pc.main
