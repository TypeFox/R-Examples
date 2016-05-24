# History Mar 11 2008  Initial coding
#         Mar 13 2008  Add option to transform the locations 
#         Mar 18 2008  Check for errors in the beginning
#         Mar 27 2008  Add option to split the screen in 2.
#         Apr 08 2008  Add option for the maximum upper limit and
#                      a dashed line.
#         May 13 2008  Add code for UNIX.
#         Jul 25 2008  Add QQ.plot
#         Aug 21 2008  In QQ.plot, plot on -log10 scale
#         Sep 03 2008  Add plotting symbols and colors 
#         Oct 15 2008  Add options for QQ.plot
#         Nov 17 2008  Generalize chromosome.plot
#         Nov 20 2008  Remove missing values in chromosome, for the
#                      chromosome plot function.
#         Nov 25 2008  Add title and legend options in chromosome plot
#         Nov 25 2008  Fix bug with matching snp names in chromosome.plot
#         Mar 17 2009  Add forest.plot function
#                      Change to setDevice
#         Apr 11 2009  Add alternating colors to chromosome plot
#         Apr 27 2009  Allow alternating colors with multiple p-values
#                      in chromsome plot.
#         Apr 29 2009  Allow for character variables in qq, chromosome plot
#         Jul 15 2009  Add optional vec to QQ plot
#         Jul 23 2009  Remove library(rmeta) call in forestPlot
#         Jul 28 2009  Let alt.colors = 1 be default value in chromosome.plot
#                      Move forest.plot to wga_plot2.R
#         Oct 13 2009  Include add option in chromosome.plot
#                      Dashed line option
#                      Add options for padj and las in chromosome.plot
#         Oct 14 2009  Add getColors function
#         Nov 19 2009  Fix bug in chrm.plot.ylim
#         Mar 09 2010  Add functions set.plot and save.plot
#         Mar 15 2010  Add new option to chromosome plot for removing white space
#                      on edges of plot.
#         Mar 23 2010  Add option to not print x-axis label in chromosome.plot
#         Mar 24 2010  Add gene.plot, remove myplot.scan.gwaa and myplot
#         Mar 26 2010  Add OR.plot.main function
#         Apr 26 2010  Add confidence intervals to OR.plot
#         May 03 2010  Add title option and list of sublists of options for OR.plot
#         Sep 23 2010  Add QQ.plot2 function for qq-plots using -2log(pval)
#         Jan 31 2011  Add wrapper function for gene.plot for perm package
#         Nov 04 2011  Fix OR.plot for 1 method, length(levels2)=1
#         Nov 15 2011  Change name of OR.plot to snp.effects.plot

# Function to set a graphics device
setDevice <- function(file, op=NULL) {

  # file     NULL or file to save plot
  ###############################################################
  # op
  #  type    "jpeg", "ps" or "pdf"
  #          The default is "jpeg"
  #  which   0 or 1  0=before plot, 1 = after plot
  #          The default is 0

  if (is.null(file)) return(NULL)

  op   <- default.list(op, c("type", "which"), list("jpeg", 0))
  type <- tolower(op$type)

  winFlag <- (.Platform$OS.type == "windows")

  if (op$which == 0) {
    if (winFlag) {
      return(0)
    }
    if (type == "ps") {
      postscript(file)
    } else if (type == "pdf") {
      pdf(file)
    } else {
      jpeg(file)
    } 
  } else {
    if (!winFlag) {
      graphics.off()
    } else {
      savePlot(filename=file, type=type)
    }
  }

  0

} # END: setDevice

# Function to create a QQ plot on a -log10 scale.
QQ.plot <- function(pvals, op=NULL) {

  # pvals      Vector of p-values
  ###############################################
  # op         List with the optional fields:
  #   title    Title for the plot
  #            The default is "QQ Plot"
  #   ylim     NULL or vector of length 2.
  #            ylim=c(0, 10) will cause the y-axis range to be between
  #            10^{-0} and 10^{-10}
  #            The default is NULL
  #   outfile  File to save the plot or NULL.
  #            The default is NULL
  #   type     "jpeg", "ps", or "pdf"
  #            The default is "jpeg"
  ###############################################

  op <- default.list(op, c("title", "type", "color"), 
                     list("QQ Plot", "jpeg", "blue"))
  op$type <- tolower(op$type)

  pvals <- as.numeric(pvals)

  # Remove non-finite values
  pvals <- pvals[is.finite(pvals)]

  pvals <- sort(pvals)
  n     <- length(pvals)
  ex    <- (1:n)/(n+1)

  # Plot on -log10 scale
  pvals[pvals < 1e-16] <- 1e-16
  pvals <- -log10(pvals)
  ex    <- -log10(ex)

  xlabel <- expression(-log[10]("Expected P-values"))
  ylabel <- expression(-log[10]("Observed P-values"))

  temp <- setDevice(op$outfile, op=list(type=op$type, which=0)) 

  plot(ex, pvals, axes=FALSE, type="p", xlab=xlabel, ylab=ylabel,
        lwd=1, col=op$color[1], ylim=op$ylim)

  maxy <- ceiling(max(pvals))
  maxx <- ceiling(max(ex))
  m    <- max(maxy, maxx)
  my   <- m
  if (!is.null(op$ylim)) my <- max(m, op$ylim)
  ypos <- 0:my
  xpos <- 0:m

  axis(1, at=xpos)
  axis(2, at=ypos)
  box()
  title(main=op$title, col.main="blue", cex.main=1)

  # Add a diagonal line
  abline(a=0, b=1, col="black", lwd=2)

  # See if other p-values need to be plotted
  pvals <- getListName(op, "pvalues")
  if (!is.null(pvals)) {
    pvals <- as.numeric(pvals)
    pvals <- sort(pvals)
    n     <- length(pvals)
    ex    <- (1:n)/(n+1)
    pvals[pvals < 1e-16] <- 1e-16
    pvals <- -log10(pvals)
    ex    <- -log10(ex)
    points(ex, pvals, type="p", lwd=1, col=op$color[2])
  }

  temp <- setDevice(op$outfile, op=list(type=op$type, which=1)) 

  0 

} # END: QQ.plot

# Function to create a QQ plot.
QQ.plot2 <- function(pvals, op=NULL) {

  # pvals      Vector of p-values
  ###############################################
  # op         List with the optional fields:
  #   title    Title for the plot
  #            The default is "QQ Plot"
  #   ylim     NULL or vector of length 2.
  #            ylim=c(0, 10) will cause the y-axis range to be between
  #            10^{-0} and 10^{-10}
  #            The default is NULL
  #   outfile  File to save the plot or NULL.
  #            The default is NULL
  #   type     "jpeg", "ps", or "pdf"
  #            The default is "jpeg"
  ###############################################

  op <- default.list(op, c("title", "type", "color"), 
                     list("QQ Plot", "jpeg", "blue"))
  op$type <- tolower(op$type)

  pvals <- as.numeric(pvals)

  # Remove non-finite values
  pvals <- pvals[is.finite(pvals)]

  pvals <- sort(pvals)
  n     <- length(pvals)
  ex    <- (1:n)/(n)

  # Plot on -log10 scale
  pvals[pvals < 1e-16] <- 1e-16
  pvals <- -2*log(pvals)
  ex    <- qchisq(ex, df=2, lower.tail=FALSE)

  xlabel <- expression("Expected P-values from chi-squared(2) distribution")
  ylabel <- expression("-2log(Observed P-values)")

  temp <- setDevice(op$outfile, op=list(type=op$type, which=0)) 

  plot(ex, pvals, axes=FALSE, type="p", xlab=xlabel, ylab=ylabel,
        lwd=1, col=op$color[1], ylim=op$ylim)

  maxy <- ceiling(max(pvals))
  maxx <- ceiling(max(ex))
  m    <- max(maxy, maxx)
  my   <- m
  if (!is.null(op$ylim)) my <- max(m, op$ylim)
  ypos <- 0:my
  xpos <- 0:m

  axis(1, at=xpos)
  axis(2, at=ypos)
  box()
  title(main=op$title, col.main="blue", cex.main=1)

  # Add a diagonal line
  abline(a=0, b=1, col="black", lwd=2)

  # See if other p-values need to be plotted
  #pvals <- getListName(op, "pvalues")
  #if (!is.null(pvals)) {
  #  pvals <- as.numeric(pvals)
  #  pvals <- sort(pvals)
  #  n     <- length(pvals)
  #  ex    <- (1:n)/(n+1)
  #  pvals[pvals < 1e-16] <- 1e-16
  #  pvals <- -log10(pvals)
  #  ex    <- -log10(ex)
  #  points(ex, pvals, type="p", lwd=1, col=op$color[2])
  #}

  temp <- setDevice(op$outfile, op=list(type=op$type, which=1)) 

  0 

} # END: QQ.plot2

# Function to transform locations for a chromosome plot
transform.loc <- function(chromosome, chrms, map, addToEnd=0) {

  nchrm <- length(chrms)

  # Get the maximum value and length of each chrm
  clen <- rep(NA, times=nchrm)
  cmax <- clen
  for (i in 1:nchrm) {
    temp <- map[(chromosome == chrms[i])]

    # Get the max and min values
    cmax[i] <- max(temp, na.rm=TRUE)
    clen[i] <- cmax[i] - min(temp, na.rm=TRUE) 
  }

  # Remove problem values
  xx <- !is.finite(cmax)
  if (any(xx)) {
    for (i in 1:nchrm) {
      if (xx[i]) {
        temp <- !(chromosome == chrms[i])
        chromosome <- chromosome[temp]
        map        <- map[temp]
      }
    }
    chrms <- chrms[!xx]
  }

  # Scale the lengths
  clen  <- clen/max(clen)
  scale <- clen/cmax 

  # Scale for each chromosome
  chpos <- rep(NA, nchrm)
  add   <- 0
  cmin  <- rep(NA, nchrm)

  for (i in 1:nchrm) {
    temp <- (chromosome == chrms[i])

    # Scale between 0 and 1 and translate
    map[temp] <- scale[i]*map[temp] + add

    # Get min and max
    cmin[i] <- min(map[temp], na.rm=TRUE)
    cmax[i] <- max(map[temp], na.rm=TRUE)

    # Get the new mean
    chpos[i] <- (cmin[i] + cmax[i])/2

    # Update add
    add <- cmax[i] + addToEnd
  }

  list(map=map, chpos=chpos, cmin=cmin, cmax=cmax)

} # END: transform.loc

# Function to create a chromosome plot
chromosome.plot <- function(infile, plot.vars, locusMap.list, op=NULL) {

  # infile           Input data set containing the p-values.
  #                  infile can be the path to the file or a data
  #                  frame containing the pvalues and the locus map
  #                  variables.
  #                  No default.
  # plot.vars        Variables in infile to plot
  # locusMap.list    List containing info about the file containing, 
  #                  snp, chromosome, and location.
  ###################################################################
  # op               List of options:
  #  splitScreen     0 or 1 to split the plot into 2 parts
  #                  The default is 1.
  #  snp.var         The snp variable name in infile
  #                  The default is "SNP".
  #  yaxis.range     Range for the y-axis. Should be on the original scale.
  #                  Ex: c(1e-12, 1e-3)
  #                  The default is NULL.
  #  title           NULL or title of plot
  #                  The default is NULL
  #  legend          0 or 1 for a legend
  #                  The default is 1
  #  legend.names    The default is plot.vars
  #  legend.horiz    TRUE or FALSE for a horizontal legend
  #                  The default is 1
  #  legend.where    "bottomright", "bottom", "bottomleft", 
  #                         "left", "topleft", "top", "topright", 
  #                         "right" and "center". 
  #                  The default is "top"
  #  cex.axis        X-axis label size
  #                  The default is 0.75
  #  colors          Vector of colors to use
  #                  The default is NULL
  #  alt.colors      0 or 1 to alternate colors in plot
  #                  The default is 1
  #  pch             Vector of plotting symbols
  #                  The default is 21 ("circles")
  #                  19 = solid circle, 20 = bullet
  #  add             A number to add spacing between the chromsomes
  #                  The default is 0.
  #  x.padj          X-axis padj option
  #                  The default depends on x.las
  #  x.las           0-3 for x-axis labels 
  #                  0=parallel, 1=horizontal, 2=perpendicular, 3=vertical
  #                  The default is 0
  #  xlim.add        Vector of length 1 or 2 for adding(subtracting) a 
  #                  value to xlim
  #                  The default is c(0, 0)
  #  xlab            The default is "Chromosome" or "Map Position"
  #####################################################################
  #  hline           NULL or list specifying a horizontal line
  #    h             y value
  #    lty
  #                  The default is NULL
  #####################################################################

  op      <- default.list(op, 
            c("splitScreen", "snp.var", "legend", "cex.axis", "alt.colors",
              "legend.horiz", "legend.where", "add", "x.las", "xlim.add"), 
             list(1, "SNP", 1, 0.75, 1, "TRUE", "top", 0, 0, c(0, 0)))
  subset  <- getListName(op, "subset")
  subFlag <- !is.null(subset)

  plot.vars <- unique(plot.vars)
  nvars     <- length(plot.vars)

  dfFlag  <- is.data.frame(infile) | is.matrix(infile)

  if (!dfFlag) {
    # Check the locusMap list
    locusMap.list <- check.locusMap.list(locusMap.list)

    # Read in the locus map data
    snp  <- NULL
    chrm <- NULL
    map  <- NULL
    for (file in locusMap.list$file) {
      temp <- paste(locusMap.list$dir, file, sep="")
      data <- getLocusMap(temp, locusMap.list)
      snp  <- c(snp, data$snp)
      chrm <- c(chrm, data$chrm)
      map  <- c(map, data$loc)
    }
    rm(data)
    temp <- gc(verbose=FALSE)

  } else {
    infile <- unfactor.all(infile)
    snp    <- infile[, locusMap.list$snp.var]
    chrm   <- as.character(infile[, locusMap.list$chrm.var])
    map    <- as.numeric(infile[, locusMap.list$loc.var])
  }

  if (subFlag) {
    temp <- chrm %in% subset
    chrm <- chrm[temp]
    snp  <- snp[temp]
    map  <- map[temp]
  }

  if (!dfFlag) {
    # Read in the p-values
    tlist <- list(file=infile, file.type=3, header=1, delimiter="\t")
    data  <- getColumns(tlist, c(op$snp.var, plot.vars), temp.list=NULL)

    # Match the snp names
    snp2 <- data[[op$snp.var]]
    data[[op$snp.var]] <- NULL
    rm(tlist)
  } else {
    snp2 <- unfactor(infile[, op$snp.var])
    data <- list()
    for (var in plot.vars) data[[var]] <- as.numeric(infile[, var])
    rm(infile)
    temp <- gc(verbose=FALSE)
  }

  # Get the correct subset for data
  temp  <- snp2 %in% snp
  snp2  <- snp2[temp]
  for (var in plot.vars) {
    data[[var]] <- as.numeric(data[[var]][temp])
  }

  # Match chrm and location to the data
  temp  <- match(snp2, snp)
  if (any(is.na(temp))) stop("ERROR: matching SNP names")
  chrm  <- chrm[temp]
  map   <- map[temp]

  maxp <- rep(NA, times=nvars)
  # Transform p-values to a -log10 scale
  for (i in 1:nvars) { 
    temp <- data[[i]] < 1e-30
    data[[i]][temp] <- 1e-30
    data[[i]]       <- -log10(data[[i]])
    maxp[i]         <- max(data[[i]], na.rm=TRUE)
  } 

  # Remove values
  temp <- (is.finite(map)) & (!is.na(chrm))
  for (i in 1:nvars) temp <- (temp & is.finite(data[[i]]))
  map  <- map[temp]
  chrm <- chrm[temp]
  for (i in 1:nvars) data[[i]] <- data[[i]][temp]

  # y-axis limits
  temp <- getListName(op, "yaxis.range")
  ylim <- chrm.plot.ylim(maxp, ylim=temp)

  rm(maxp, var, snp, snp2)
  temp <- gc(verbose=FALSE)
  uchrm <- getUniqueChrm(chrm)
  nchrm <- length(uchrm)

  if (nchrm <= 1) {
    op$splitScreen <- 0
    transLoc       <- 0
  } else {
    transLoc       <- 1
  }

  # axis labels
  xlab <- op[["xlab", exact=TRUE]]
  if (is.null(xlab)) {
    if (nchrm > 1) {
      xlab <- "Chromosome"
    } else {
      xlab <- "Map Position"
    }
  }
  ylab <- expression(-log[10](P-value))

  # Legend
  nclrs <- nvars
  alt.colors <- op$alt.colors
  if (alt.colors) nclrs <- nchrm
  colors <- getListName(op, "colors")
  colors <- chrm.plot.colors(nclrs, colors=colors)
  pch    <- chrm.plot.pch(nvars, pch=op$pch)
  if (op$legend) {
    names <- plot.vars
    temp <- getListName(op, "legend.names")
    if (!is.null(temp)) names <- temp
    legend <- chrm.plot.legend(names, colors, pch, alt.colors,
               horiz=op$legend.horiz, where=op$legend.where)
  } else {
    legend <- NULL
  }

  hline  <- getListName(op, "hline")
  x.padj <- getListName(op, "x.padj") 

  if (op$splitScreen) {

    half  <- floor(nchrm/2)
    chrm1 <- uchrm[1:half]

    # Get the ids
    temp  <- chrm %in% chrm1
    
    # Top plot
    scr <- split.screen(c(2,1))
    screen(1)

    ret <- chrm.plot.main(data, map, chrm, ids=temp, ylim=ylim[1:2],
           xlab=xlab, ylab=ylab, transLoc=transLoc, legend=legend,
           colors=colors, pch=op$pch, title=op$title, cex.axis=op$cex.axis,
           alt.colors=alt.colors, add=op$add, hline=hline,
           x.las=op$x.las, x.padj=x.padj, xlim.add=op$xlim.add)
    #close.screen(scr[1])

    # Get the ids for the bottom plot
    temp <- as.logical(1-temp)

    screen(2)
    ret <- chrm.plot.main(data, map, chrm, ids=temp, ylim=ylim[3:4],
            xlab=xlab, ylab=ylab, transLoc=transLoc, legend=legend,
            colors=colors, pch=op$pch, title=op$title, cex.axis=op$cex.axis,
            alt.colors=alt.colors, add=op$add, hline=hline,
            x.las=op$x.las, x.padj=x.padj, xlim.add=op$xlim.add)
    #close.screen(scr[2])

  } # END: if (splitScreen)
  else {
    ret <- chrm.plot.main(data, map, chrm, ids=NULL, ylim=ylim[1:2],
              xlab=xlab, ylab=ylab, transLoc=transLoc, legend=legend,
            colors=colors, pch=op$pch, title=op$title, cex.axis=op$cex.axis,
            alt.colors, add=op$add, hline=hline,
            x.las=op$x.las, x.padj=x.padj, xlim.add=op$xlim.add)
  }

} # END: chromosome.plot

# Function to produce a plot
chrm.plot.main <- function(pvals, map, chrm, ids=NULL, ylim=c(0, 8),
                  xlab=NULL, ylab=NULL, transLoc=1, legend=NULL,
                  colors=NULL, pch=NULL, title=NULL, cex.axis=0.75,
                  alt.colors=0, add=0, hline=NULL, x.las=0, x.padj=-1.0,
                  xlim.add=c(0,0)) {

  # pvals   List
 
  if (is.null(ids)) ids <- 1:length(map)
  if (is.null(x.padj)) {
    if (x.las %in% c(0, 1)) {
      x.padj <- -1.0
    } else {
      x.padj <- 0.5
    }
  }

  # Subset the data
  chrm  <- chrm[ids]
  map   <- map[ids]
  vars  <- names(pvals)
  nvars <- length(vars)

  for (i in 1:nvars) pvals[[i]] <- pvals[[i]][ids]

  rm(ids)
  temp <- gc(verbose=FALSE)
  
  uchrms <- getUniqueChrm(chrm)
  nchrms <- length(uchrms)

  if (transLoc) {
    temp   <- transform.loc(chrm, uchrms, map, addToEnd=add)
    map    <- temp$map
    chpos  <- temp$chpos
    cmin   <- temp$cmin
    cmax   <- temp$cmax
  }

  nclrs <- nvars
  if (alt.colors) nclrs <- nchrms

  # Get the colors and plotting characters
  col  <- chrm.plot.colors(nclrs, colors=colors)
  pch  <- chrm.plot.pch(nvars, pch=pch)
  if (length(xlim.add) == 1) xlim.add <- rep(xlim.add, times=2)
  xlim <- c(min(map, na.rm=TRUE)+xlim.add[1], max(map, na.rm=TRUE)+xlim.add[2])

  # Initialize the plot
  if (transLoc) {
    
    plot(map[1], pvals[[1]][1], xlab=xlab,ylab=ylab, axes=FALSE, xlim=xlim,
                ylim=ylim, col=col[1], pch=pch[1])
    mxlog <- floor(ylim[2])
    #if (mxlog == 0) mxlog <- maxp[1]
    if (mxlog < 1) {
      axis(2, at=c(ylim[1], mxlog))
    } else {
      axis(2, at=floor(ylim[1]:mxlog))
    }

    # tcl is for tick mark length
    # padj is for moving the labels close/farther from the axis
    axis(1, at=chpos, labels=uchrms, tcl=-0.5, padj=x.padj, las=x.las,
         cex.axis=cex.axis, font=2)
    box()   

  } else {
    plot(map, pvals[[1]], xlab=xlab,ylab=ylab, axes=TRUE, xlim=xlim,
                ylim=ylim, col=col[1], pch=pch[1])
  }

  # Plot the other vars
  if (nvars > 1) {
    for (i in 2:nvars) {
      points(map, pvals[[i]], col=col[i], pch=pch[i])
    }
  }

  # Alternate colors
  if (alt.colors) {
    for (j in 1:nvars) {
      for (i in 1:nchrms) {
        temp <- (chrm == uchrms[i])
        points(map[temp], pvals[[j]][temp], pch=pch[j], col=col[i])
      }
    }
  }

  # Add line segments (tick marks)
  if ((transLoc) & (!alt.colors)) addLineSegments(cmin, cmax, ylim) 

  # Add a legend
  if (!is.null(legend)) {
    legend(x=legend$x, y=legend$y, legend=legend$legend,
           bty=legend$bty, title=legend$title,
           pch=legend$pch, col=legend$col, cex=legend$cex, 
           horiz=legend$horiz)
  }

  # Add title 
  if (!is.null(title)) title(title)

  if (!is.null(hline)) {
    hline <- default.list(hline, c("h", "lty"), list(7, 2))
    abline(h=hline$h, lty=hline$lty)
  }

} # END: chrm.plot.main

# Function to get the unique chromosomes
getUniqueChrm <- function(chrm) {

  u    <- unique(chrm)
  n    <- as.numeric(u)
  nas  <- is.na(n)
  ch   <- u[nas]
  s    <- sort(n)
  c(s, ch)

} # END: getUniqueChrm

# Function to return a vector of colors
chrm.plot.colors <- function(n, colors=NULL) {

  if (is.null(colors)) {
    col    <- c("blue", "pink", "green", "red", "black", "orange",
                "yellow", "gray", "brown", "turquoise", "gold",
                "violet", "olivedrab", "skyblue", "purple",
                "maroon")
    col <- c(col, col)
    #col <- rainbow(n)
  } else {
    col <- rep(colors, times=n)
  }
  col    <- col[1:n]
  col

} # END: chrm.plot.colors

# Function to return a vector of plotting characters
chrm.plot.pch <- function(n, pch=NULL) {

  if (is.null(pch)) {
    pch <- rep(21, times=n)
  } else {
    
    pch <- rep(pch, times=n)
  }
  pch <- pch[1:n]
  pch

} # END: chrm.plot.pch

# Function to return the legend list
chrm.plot.legend <- function(plot.vars, colors, pch, alt.colors,
                             horiz=TRUE, where="top") {

  nvars  <- length(plot.vars)
  col    <- chrm.plot.colors(nvars, colors=colors)
  if ((alt.colors) && (nvars > 1)) col <- rep("black", times=nvars)
  title  <- NULL
  x      <- where
  y      <- NULL
  legend <- plot.vars
  cex    <- 0.7
  horiz  <- horiz

  list(x=x, y=y, bty="n", pch=pch, col=col, title=title, 
       legend=legend, cex=cex, horiz=horiz)  

} # END: chrm.plot.legend

# Function to add line segments to a graph
addLineSegments <- function(cmin, cmax, ylim) {

  n    <- length(cmin) + 1
  temp <- cmin 
  temp[1] <- cmin[1] - 0.035
  temp[n] <- cmax[n-1] + 0.035
  if (n > 2) {
    for (i in 2:(n-1)) {
      temp[i] <- (cmax[i-1] + cmin[i])/2
    }
  }
  temp0 <- rep.int(-1, times=n)
  temp1 <- rep.int(ylim[1]+0.1, times=n)
  segments(temp, temp0, temp, temp1, lwd=1, font=2)

} # END: addLineSegments

# Function to return the y-axis limits
chrm.plot.ylim <- function(maxp, ylim=NULL) {

  # maxp is on a -log10 scale
  temp <- max(maxp) + 1
  y0   <- c(0, temp, 0, temp)
  if (is.null(ylim)) return(y0)
  n <- length(ylim)
  if (!(n %in% c(2, 4))) return(y0)
  
  for (i in 1:n) {
    if ((ylim[i] > 0) && (ylim[i] <= 1)) ylim[i] <- -log10(ylim[i])
  }

  # check ylim
  if (ylim[1] > ylim[2]) {
    temp    <- ylim[1]
    ylim[1] <- ylim[2]
    ylim[2] <- temp
  }
  if (n == 4) {
    if (ylim[3] > ylim[4]) {
      temp    <- ylim[3]
      ylim[3] <- ylim[4]
      ylim[4] <- temp
    }
  }

  if (n == 2) ylim <- c(ylim, ylim)

  ylim

} # END: chrm.plot.ylim

# Function to return colors
getColors <- function(colVec, op=NULL) {

 # op      List with names
 #  ncolors  (Maximum) number of colors to return
 #           Default is NULL
 #  plot     0 or 1
 #  print    0 or 1
 #  seed     Default is -1
 #  exclude  Character vector of colors to exclude

 op <- default.list(op, c("plot", "print", "seed", "sample"), list(0, 0, -1, 0))

 all <- colors()
 cls <- NULL
 for (cc in colVec) {
   temp <- grep(cc, all, fixed=TRUE)
   if (length(temp)) cls <- c(cls, all[temp])
 }
 
 rm(all)
 gc()
 exclude <- getListName(op, "exclude")
 if (!is.null(exclude)) {
   for (cc in exclude) {
     temp <- grep(cc, cls, fixed=TRUE)
     if (length(temp)) cls <- cls[-temp]
   }
 } 

 seed <- op$seed
 if (seed > 0) set.seed(seed)
 n <- length(cls)
 ncolors <- getListName(op, "ncolors")
 if (!is.null(ncolors)) {
   # If n != ncolors, sample from the colors
   if (n != ncolors) {
     if (n > ncolors) {
       cls <- sample(cls, ncolors, replace=FALSE)
     } else {
       temp <- ncolors - n
       if (temp > n) {
         temp <- sample(cls, temp, replace=TRUE)
       } else {
         temp <- sample(cls, temp, replace=FALSE)
       }
       cls <- c(cls, temp) 
     }
   }
 } else {
   ncolors <- n
 } 

 if (op$sample) cls <- sample(cls, ncolors, replace=FALSE)
 if (op$print) print(cls)
 if (op$plot) pie(rep.int(1, ncolors), col=cls)

 cls

} # END: getColors

# Function to call before plotting
set.plot <- function(op=NULL) {

  winFlag <- (.Platform$OS.type == "windows")
  if (winFlag) return(0)

  op <- default.list(op, c("out", "type", "res"),
           list("ERROR", "jpeg", 100), error=c(1, 0, 0))

  host <- callOS("hostname", intern=TRUE)
  temp <- op[["R_GSCMD", exact=TRUE]]
  if (is.null(temp)) {
    if (host == "biowulf.nih.gov") {
      temp <- "/data/wheelerb/nilanjan/wga/software/ghostscript-8.64/bin/gs" 
    } else {
      temp <- "/data/software/ghostscript-8.64/bin/gs"
    }
  }

  # Set the R_GSCMD environment variable for plots
  Sys.setenv(R_GSCMD=temp)

  bitmap(op$out, type=op$type, res=op$res)

  0 

} # END: set.plot

# Function to call after plotting
save.plot <- function(op=NULL) {

  winFlag <- (.Platform$OS.type == "windows")
  if (!winFlag) {
    #dev.off()
    return(0)
  }

  op <- default.list(op, c("out", "type"),
           list("ERROR", "jpeg"), error=c(1, 0))

  savePlot(filename=op$out, type=op$type)

  0 

} # END: save.plot

# Function for a gene plot
gene.plot <- function(infile, plot.vars, locusMap.list, op=NULL) {

  # infile           Input data set containing the p-values.
  #                  infile can be the path to the file or a data
  #                  frame containing the pvalues and the locus map
  #                  variables.
  #                  No default.
  # plot.vars        Variables in infile to plot
  # locusMap.list    List containing info about the file containing, 
  #                  snp, gene, chromosome, location.
  ###################################################################
  # op               List of options:
  #  splitScreen     0 or 1 to split the plot into 2 parts
  #                  The default is 1.
  #  snp.var         The snp variable name in infile
  #                  The default is "SNP".
  #  yaxis.range     Range for the y-axis. Should be on the original scale.
  #                  Ex: c(1e-12, 1e-3)
  #                  The default is NULL.
  #  title           NULL or title of plot
  #                  The default is NULL
  #  subtitle        NULL or sub-title of plot
  #                  The default is "Gene"
  #  legend          0 or 1 for a legend
  #                  The default is 0
  #  legend.names    The default is plot.vars
  #  legend.horiz    TRUE or FALSE for a horizontal legend
  #                  The default is 1
  #  legend.where    "bottomright", "bottom", "bottomleft", 
  #                         "left", "topleft", "top", "topright", 
  #                         "right" and "center". 
  #                  The default is "top"
  #  cex.axis        X-axis label size
  #                  The default is 0.75
  #  colors          Vector of colors to use
  #                  The default is NULL
  #  pch             Vector of plotting symbols
  #                  The default is 20 ("circles")
  #                  19 = solid circle, 20 = bullet
  #  add             A number to add spacing between the chromsomes
  #                  The default is 0.
  #  x.padj          X-axis padj option
  #                  The default depends on x.las
  #  x.las           0-3 for x-axis labels 
  #                  0=parallel, 1=horizontal, 2=perpendicular, 3=vertical
  #                  The default is 2
  #  xlim.add        Vector of length 1 or 2 for adding(subtracting) a 
  #                  value to xlim
  #                  The default is c(0, 0)
  #  xlab            The default is ""
  #  subset          Character vector of genes to plot
  #                  The default is NULL
  #  chr.text        For the chromosome numbers
  #                  The default is NULL
  #  maxLabelLen     Maximum length of x-axis labels
  #                  The default is NULL
  #  uniform         0 or 1 for uniform spacing of genes in plot
  #                  The default is 1
  #####################################################################
  #  hline           NA or list specifying a horizontal line
  #    h             y value The default is 10e-7
  #    lty
  #                  The default is NULL
  #####################################################################

  op      <- default.list(op, 
            c("splitScreen", "snp.var", "legend", "cex.axis", "alt.colors",
              "legend.horiz", "legend.where", "add", "x.las", "xlim.add", "pch",
              "xlab", "subtitle", "uniform"), 
             list(0, "SNP", 0, 0.75, 1, "TRUE", "top", 0.1, 2, c(0, 0), 20, 
                  "", "Gene", 1))
  subset  <- op[["subset", exact=TRUE]]
  subFlag <- !is.null(subset)

  plot.vars <- unique(plot.vars)
  nvars     <- length(plot.vars)

  dfFlag  <- is.data.frame(infile) | is.matrix(infile)

  locusMap.list <- default.list(locusMap.list, c("gene.var"), list("ERROR"), error=c(1))

  if (!dfFlag) {
    # Check the locusMap list
    locusMap.list <- check.locusMap.list(locusMap.list)

    # Read in the locus map data
    snp  <- NULL
    chrm <- NULL
    map  <- NULL
    gene <- NULL
    for (file in locusMap.list$file) {
      temp <- paste(locusMap.list$dir, file, sep="")
      data <- getLocusMap(temp, locusMap.list)
      snp  <- c(snp, data$snp)
      chrm <- c(chrm, data$chrm)
      map  <- c(map, data$loc)
      gene <- c(gene, data$gene)
    }
    rm(data)
    temp <- gc(verbose=FALSE)

  } else {
    infile <- unfactor.all(infile)
    snp    <- infile[, locusMap.list$snp.var]
    chrm   <- as.character(infile[, locusMap.list$chrm.var])
    map    <- as.numeric(infile[, locusMap.list$loc.var])
    gene   <- infile[, locusMap.list$gene.var]
  }

  if (subFlag) {
    temp <- gene %in% subset
    gene <- gene[temp]
    chrm <- chrm[temp]
    snp  <- snp[temp]
    map  <- map[temp]
  }

  if (!dfFlag) {
    # Read in the p-values
    tlist <- list(file=infile, file.type=3, header=1, delimiter="\t")
    data  <- getColumns(tlist, c(op$snp.var, plot.vars), temp.list=NULL)

    # Match the snp names
    snp2 <- data[[op$snp.var]]
    data[[op$snp.var]] <- NULL
    rm(tlist)
  } else {
    snp2 <- unfactor(infile[, op$snp.var])
    data <- list()
    for (var in plot.vars) data[[var]] <- as.numeric(infile[, var])
    rm(infile)
    temp <- gc(verbose=FALSE)
  }

  # Get the correct subset for data
  temp  <- snp2 %in% snp
  snp2  <- snp2[temp]
  for (var in plot.vars) {
    data[[var]] <- as.numeric(data[[var]][temp])
  }

  # Match chrm and location to the data
  temp  <- match(snp2, snp)
  if (any(is.na(temp))) stop("ERROR: matching SNP names")
  chrm  <- chrm[temp]
  map   <- map[temp]
  gene  <- gene[temp]

  maxp <- rep(NA, times=nvars)
  # Transform p-values to a -log10 scale
  for (i in 1:nvars) { 
    temp <- data[[i]] < 1e-30
    data[[i]][temp] <- 1e-30
    data[[i]]       <- -log10(data[[i]])
    maxp[i]         <- max(data[[i]], na.rm=TRUE)
  } 

  # Remove values
  temp <- (is.finite(map)) & (!is.na(chrm)) & (!is.na(gene))
  for (i in 1:nvars) temp <- (temp & is.finite(data[[i]]))
  map  <- map[temp]
  chrm <- chrm[temp]
  gene <- gene[temp]
  for (i in 1:nvars) data[[i]] <- data[[i]][temp]

  # y-axis limits
  temp <- getListName(op, "yaxis.range")
  ylim <- chrm.plot.ylim(maxp, ylim=temp)

  rm(maxp, var, snp, snp2)
  temp <- gc(verbose=FALSE)
  uchrm <- getUniqueChrm(chrm)

  nchrm <- length(uchrm)
  ngene <- length(unique(gene))

  # Get the genes for each chromosome and order them
  glist <- gene.getOrder(chrm, gene, map)

  if (nchrm <= 1) {
    op$splitScreen <- 0
    transLoc       <- 1
  } else {
    transLoc       <- 1
  }

  # axis labels
  xlab <- op[["xlab", exact=TRUE]]
  if (is.null(xlab)) {
    if (nchrm > 0) {
      xlab <- "Gene"
    } else {
      xlab <- "Map Position"
    }
  }
  xlab <- ""
  ylab <- expression(-log[10](P-value))

  # Legend
  colors <- op[["colors", exact=TRUE]]
  colors <- gene.plot.colors(glist, colors=colors)
  pch    <- chrm.plot.pch(nvars, pch=op$pch)
  if (op$legend) {
    names <- plot.vars
    temp <- getListName(op, "legend.names")
    if (!is.null(temp)) names <- temp
    legend <- chrm.plot.legend(names, colors, pch, NULL,
               horiz=op$legend.horiz, where=op$legend.where)
  } else {
    legend <- NULL
  }

  hline  <- getListName(op, "hline")
  x.padj <- getListName(op, "x.padj") 

  if (op$splitScreen) {

    half  <- floor(nchrm/2)
    chrm1 <- uchrm[1:half]
    chrm2 <- uchrm[(half+1):nchrm]

    # Get the ids
    temp   <- chrm %in% chrm1
    glist1 <- glist[chrm1]
    glist2 <- glist[chrm2]
    
    rm(chrm, chrm1)
    gc()

    # Top plot
    split.screen(c(2,1))
    screen(1)

    ret <- gene.plot.main(data, map, gene, glist1, ids=temp, ylim=ylim[1:2],
           xlab=xlab, ylab=ylab, transLoc=transLoc, legend=legend,
           colors=colors, pch=op$pch, title=op$title, cex.axis=op$cex.axis,
           alt.colors=op$alt.colors, add=op$add, hline=hline,
           x.las=op$x.las, x.padj=x.padj, xlim.add=op$xlim.add,
           maxLabelLen=op$maxLabelLen, subtitle=op$subtitle, chr.text=op$chr.text,
           uniform=op$uniform)

    # Get the ids for the bottom plot
    temp <- as.logical(1-temp)

    screen(2)
    ret <- gene.plot.main(data, map, gene, glist2, ids=temp, ylim=ylim[3:4],
            xlab=xlab, ylab=ylab, transLoc=transLoc, legend=legend,
            colors=colors, pch=op$pch, title=op$title, cex.axis=op$cex.axis,
            alt.colors=op$alt.colors, add=op$add, hline=hline,
            x.las=op$x.las, x.padj=x.padj, xlim.add=op$xlim.add,
            maxLabelLen=op$maxLabelLen, subtitle=op$subtitle, chr.text=op$chr.text,
            uniform=op$uniform)

  } # END: if (splitScreen)
  else {
    rm(chrm)
    gc()

    ret <- gene.plot.main(data, map, gene, glist, ids=NULL, ylim=ylim[1:2],
              xlab=xlab, ylab=ylab, transLoc=transLoc, legend=legend,
            colors=colors, pch=op$pch, title=op$title, cex.axis=op$cex.axis,
            alt.colors=op$alt.colors, add=op$add, hline=hline,
            x.las=op$x.las, x.padj=x.padj, xlim.add=op$xlim.add,
            maxLabelLen=op$maxLabelLen, subtitle=op$subtitle, chr.text=op$chr.text,
            uniform=op$uniform)
  }

} # END: gene.plot

# Function to produce a plot
gene.plot.main <- function(pvals, map, gene, glist, ids=NULL, ylim=c(0, 8),
                  xlab="", ylab=NULL, transLoc=1, legend=NULL,
                  colors=NULL, pch=NULL, title=NULL, cex.axis=0.75,
                  alt.colors=0, add=0, hline=NULL, x.las=0, x.padj=-1.0,
                  xlim.add=c(0,0), maxLabelLen=NULL, subtitle="Gene", chr.text=NULL,
                  uniform=0) {

  # pvals   List
 
  if (is.null(ids)) ids <- 1:length(map)
  if (is.null(x.padj)) {
    if (x.las %in% c(0, 1)) {
      x.padj <- -1.0
    } else {
      x.padj <- 0.5
    }
  }

  # Subset the data
  map   <- map[ids]
  gene  <- gene[ids]
  vars  <- names(pvals)
  nvars <- length(vars)

  for (i in 1:nvars) pvals[[i]] <- pvals[[i]][ids]

  rm(ids)
  temp <- gc(verbose=FALSE)

  nchrms <- length(glist)
  if (transLoc) {
    temp   <- gene.transform.loc(gene, glist, map, addToEnd=add, uniform=uniform)
    map    <- temp$map
    chpos  <- temp$chpos
    cmin   <- temp$cmin
    cmax   <- temp$cmax
    midMap <- temp$midChrMap
  }

  nclrs <- nchrms

  # Get the colors and plotting characters
  col  <- gene.plot.colors(glist, colors=colors)
  pch  <- chrm.plot.pch(nvars, pch=pch)
  if (length(xlim.add) == 1) xlim.add <- rep(xlim.add, times=2)
  xlim <- c(min(map, na.rm=TRUE)+xlim.add[1], max(map, na.rm=TRUE)+xlim.add[2])

  # Initialize the plot
  if (transLoc) {
    
    plot(map[1], pvals[[1]][1], xlab=xlab,ylab=ylab, axes=FALSE, xlim=xlim,
                ylim=ylim, col=col[1], pch=pch[1])
    mxlog <- floor(ylim[2])
    #if (mxlog == 0) mxlog <- maxp[1]
    if (mxlog < 1) {
      axis(2, at=c(ylim[1], mxlog))
    } else {
      axis(2, at=floor(ylim[1]:mxlog))
    }

    # Get the labels (gene names)
    labels <- NULL
    for (i in 1:nchrms) labels <- c(labels, glist[[i]])
    if (!is.null(maxLabelLen)) labels <- substr(labels, 1, maxLabelLen)

    # tcl is for tick mark length
    # padj is for moving the labels close/farther from the axis
    axis(1, at=chpos, labels=labels, tcl=-0.5, padj=x.padj, las=x.las,
         cex.axis=cex.axis, font=2)
    box()   

  } else {
    plot(map, pvals[[1]], xlab=xlab,ylab=ylab, axes=TRUE, xlim=xlim,
                ylim=ylim, col=col[1], pch=pch[1])
  }

  # Plot the other vars
  if (nvars > 1) {
    for (i in 2:nvars) {
      points(map, pvals[[i]], col=col[i], pch=pch[i])
    }
  }

  # Alternate colors
  for (j in 1:nvars) {
    for (i in 1:nchrms) {
      gg  <- glist[[i]]
      ngg <- length(gg)
      for (k in 1:ngg) {
        temp <- (gene == gg[k])
        points(map[temp], pvals[[j]][temp], pch=pch[j], col=col[i])
      }
    }
  }

  # Add line segments (tick marks)
  if ((transLoc) & (!alt.colors)) addLineSegments(cmin, cmax, ylim) 

  # Add a legend
  if (!is.null(legend)) {
    legend(x=legend$x, y=legend$y, legend=legend$legend,
           bty=legend$bty, title=legend$title,
           pch=legend$pch, col=legend$col, cex=legend$cex, 
           horiz=legend$horiz)
  }

  # Add title 
  if (!is.null(title)) title(title, sub=subtitle)

  if (!is.null(hline)) {
    hline <- default.list(hline, c("h", "lty"), list(7, 2))
    abline(h=hline$h, lty=hline$lty)
  }

  # Add Chromosome numbers
  flag <- 1
  if (!is.null(chr.text)) {
    if (!is.list(chr.text)) flag <- 0
  }
  if (flag) {  
    chr.text <- default.list(chr.text, c("cex", "at", "text"), list(0.75, 0, "Chr"))
    mtext(names(glist), side=3, line=0, outer=FALSE, at=midMap, cex=chr.text$cex, font=2)
    mtext(chr.text$text, side=3, line=0, outer=FALSE, at=chr.text$at, cex=chr.text$cex, font=2)
  }

} # END: gene.plot.main

# Function to return a vector of colors
gene.plot.colors <- function(glist, colors=NULL) {

  cls <- colors
  n <- length(glist)
  if (!length(cls)) cls <- NULL
  if (is.null(cls)) {
    cls <- c("blue", "pink", "green", "red", "black", "orange",
                "yellow", "gray", "brown", "turquoise", "gold",
                "violet", "olivedrab", "skyblue", "purple",
                "maroon")
  } 
  lenc <- length(cls)
  if (lenc < n) {
    temp <- ceiling(n/lenc) 
    cls  <- rep(cls, times=temp)
  }
  cls <- cls[1:n]

  cls

} # END: gene.plot.colors

# Function to transform locations for a chromosome plot
gene.transform.loc <- function(gene, glist, map, addToEnd=0, uniform=0) {

  

  nchrm <- length(glist)
  ngene <- length(unique(gene))
  clen <- rep(NA, times=ngene)
  cmax <- clen

  # Get the maximum value and length of each gene
  if (!uniform) {
    ii   <- 1
    for (i in 1:nchrm) {
      gg  <- glist[[i]]
      ngg <- length(gg)
      for (j in 1:ngg) {
        temp <- map[(gene == gg[j])]

        # Get the max and min values
        cmax[ii] <- max(temp, na.rm=TRUE)
        clen[ii] <- cmax[ii] - min(temp, na.rm=TRUE)
        ii       <- ii + 1
      } 
    }
  } else {
    map[]  <- 1
    cmax[] <- 1
    clen[] <- 1
  }

  # Remove problem values
  #xx <- !is.finite(cmax)
  #if (any(xx)) {
  #  for (i in 1:nchrm) {
  #    if (xx[i]) {
  #      temp <- !(chromosome == chrms[i])
  #      chromosome <- chromosome[temp]
  #      map        <- map[temp]
  #    }
  #  }
  #  chrms <- chrms[!xx]
  #}

  # Scale the lengths
  clen  <- clen/max(clen)
  scale <- clen/cmax 

  # Scale for each chromosome
  chpos <- rep(NA, ngene)
  add   <- 0
  cmin  <- rep(NA, ngene)
  midChrMap <- rep(NA, nchrm)

  ii <- 1
  for (i in 1:nchrm) {
    gg   <- glist[[i]]
    ngg  <- length(gg)
    mmin <- 1e100
    mmax <- -999999
    for (j in 1:ngg) {
      temp <- (gene == gg[j])

      # Scale between 0 and 1 and translate
      map[temp] <- scale[ii]*map[temp] + add

      # Get min and max
      cmin[ii] <- min(map[temp], na.rm=TRUE)
      cmax[ii] <- max(map[temp], na.rm=TRUE)

      # Get the new mean
      chpos[ii] <- (cmin[ii] + cmax[ii])/2

      # Update add
      add <- cmax[ii] + addToEnd

      # Get midpoint for chromosome number
      mmin <- min(cmin[ii], mmin)
      mmax <- max(cmax[ii], mmax)

      ii  <- ii + 1
    }
    midChrMap[i] <- (mmin + mmax)/2
  }

  list(map=map, chpos=chpos, cmin=cmin, cmax=cmax, midChrMap=midChrMap)

} # END: gene.transform.loc

# Function to get the genes for each chromosome and order them
gene.getOrder <- function(chrm, gene, map) {

  glist <- list()
  #uchrm <- unique(chrm)
  uchrm <- getUniqueChrm(chrm)
  nchrm <- length(uchrm)
  for (i in 1:nchrm) {
    temp  <- (chrm == uchrm[i])
    gg    <- unique(gene[temp]) 
    ngg   <- length(gg)
    mloc  <- double(ngg)
    # Get the smallest location
    for (j in 1:ngg) {
      temp <- (gene == gg[j]) 
      mloc[j] <- min(map[temp])
    }
    temp <- sort.int(mloc, index.return=TRUE)
    gg <- gg[temp$ix]
    glist[[uchrm[i]]] <- gg
  }

  glist

} # END: gene.getOrder

# Main function for OR plot
OR.plot.main <- function(or, lower=NULL, upper=NULL, op=NULL) {

  # or        Vector or matrix of odds-ratios
  #           A vector is considered as a matrix of dimension c(1, length(or))
  # lower     Vector or matrix of lower confidence limits
  #           The default is NULL
  # upper     Vector or matrix of upper confidence limits
  #           The default is NULL
  # op        List with names
  #  by       1 or 2 for the grouping of ORs
  #           1 = rows, 2 = columns
  #           The default is 1
  #  xlab
  #  ylab
  #  title    Character string or list with names "main" and "sub"
  #  yticks   y-axis tick labels
  #  colors
  #  digits   Number of significant digits on y-axis to round to (if yticks is NULL)
  #           The default is 2
  #  xticks   x-axis tick labels. The default is the names of the or vector/matrix
  #  ylim     Vector of length 1 or 2 to set the y limits
  #           The default is NULL
  #  add      For spacing in plot
  #           The default is either 0 or 0.05 depending on other options
  #  add2     For spacing in plot
  #           The default is 1
  ##################################################################################
  #  legend   List for the legend. Set to NA for no legend.
  #           The default is NULL so that a legend will be created depending on
  #           other options selected and the input or object.
  #################################################################################

  # Function to define the polygon
  definePoly <- function(x0, y0) {

    if ((y0 <= ylim0[2]) && (y0 >= ylim0[1])) {
      x <- c(x0, x0+1, x0+1, x0)
      y <- c(1,     1,   y0, y0)
      return(list(x=x, y=y))
    }
 
    # Define the sawtooth
    x1 <- x0 + 1/3
    x2 <- x0 + 0.5
    x3 <- x0 + 2/3
    x4 <- x0 + 1
    if (y0 > 1) {
      y0 <- ylim0[2] 
      y1 <- y0 - pstep
      y2 <- y0 - pstep/4
      y3 <- y1
      y4 <- y0 - pstep/2
    } else {
      y0 <- ylim0[1] 
      y1 <- y0 + pstep
      y2 <- y0 + pstep/4
      y3 <- y1
      y4 <- y0 + pstep/2
    } 
    x <- c(x0, x1, x2, x3, x4, x4, x0)
    y <- c(y0, y1, y2, y3, y4, 1,  1)

    list(x=x, y=y)
  } # END: definePoly

  lFlag <- !is.null(lower)
  uFlag <- !is.null(upper)
  op    <- default.list(op, c("by", "add2", "digits"), list(1, 1, 2))
  if (is.null(op[["add", exact=TRUE]])) {
    if ((lFlag) || (uFlag)) {
      op$add <- 0.05
    } else {
      op$add <- 0
    }
  }

  legend <- op[["legend", exact=TRUE]]
  if (is.null(legend)) legend <- list()
  temp <- !is.na(legend[1])
  if (!length(temp)) temp <- FALSE
  legendFlag <- temp
  if (legendFlag) {
    legend <- default.list(legend, 
    c("x", "horiz", "cex", "bty", "inset", "y.intersp", "x.intersp"), 
    list("top", TRUE, 1.0, "n", -0.02, 0.75, 1))
  }

  # Set up the matrix of ORs
  d <- dim(or)
  if (is.null(d)) {
    temp <- names(or)
    len     <- length(or)
    dim(or) <- c(1, len)
    colnames(or) <- temp
    if (lFlag) {
      dim(lower) <- c(1, len)
      colnames(lower) <- temp
    }
    if (uFlag) {
      dim(upper) <- c(1, len)
      colnames(upper) <- temp
    }
  }
  if (op$by == 2) {
    or <- t(or)
    if (lFlag) lower <- t(lower)
    if (uFlag) upper <- t(upper)
  }
  d <- dim(or)
  
  temp  <- is.finite(or)
  if (uFlag) temp <- temp & is.finite(upper)
  minor <- min(or[temp], na.rm=TRUE)
  maxor <- max(or[temp], na.rm=TRUE)
  if (lFlag) minor <- min(minor, lower[temp], na.rm=TRUE)
  if (uFlag) maxor <- max(maxor, upper[temp], na.rm=TRUE)
  ymin  <- min(1, minor)
  #h.legend <- (maxor-minor)/25

  ylim  <- op[["ylim", exact=TRUE]]
  ylim0 <- ylim
  if (is.null(ylim)) {
    ylim     <- c(ymin, maxor)
    ylim0    <- ylim
    ylimFlag <- 0
  } else {
    ylimFlag <- 1
  } 
  h.legend <- (ylim[2]-ylim[1])/25 

  # Get the legend position
  if (legendFlag) {
    temp <- legend$x
    if ((!is.null(temp)) && (is.character(temp))) {
      temp <- toupper(substr(temp, 1, 1))
      if (temp == "T") {
        ylim[2] <- ylim[2] + h.legend
      } else if (temp == "B") {
        ylim[2] <- ylim[2] - h.legend
      } 
    }
  }
  pstep <- (ylim[2] - ylim[1])/50

  # Get the x limits
  add  <- op$add
  add2 <- op$add2
  len  <- length(or)
  xmax <- len + len*add + (d[2]-1)*add2
  xlim <- c(0, xmax)  

  ylab <- op[["ylab", exact=TRUE]]
  if (is.null(ylab)) ylab <- "Odds-ratio"
  xlab <- op[["xlab", exact=TRUE]]
  if (is.null(xlab)) xlab <- ""

  plot(1, 1, xlab=xlab,ylab=ylab, axes=FALSE, xlim=xlim,
                ylim=ylim, type="n")
  ylim <- ylim0

  # Get the colors
  col <- op[["colors", exact=TRUE]]
  if (is.null(col)) {
    col <- c("red", "orange", "blue", "green", "purple",
             "brown", "yellow", "pink", "olivedrab", "powderblue", "gray",
             "aquamarine", "magenta", "gold", "darkblue")
  }

  # Use polygon function for each OR
  x0   <- 0
  xavg <- rep(0, times=d[2])
  for (j in 1:d[2]) {
    for (i in 1:d[1]) {
      if (d[1] == 1) {
        a <- x0
        b <- x0 + 1 + add
      } else if (i == 1) {
        a <- x0
      } else if (i == d[1]) {
        b <- x0 + 1
      }

      if (lFlag) {
        temp <- definePoly(x0, lower[i, j])
        x    <- temp$x
        y    <- temp$y
        polygon(x, y, border=col[i], col="white")
      }
      if (uFlag) {
        temp <- definePoly(x0, upper[i, j])
        x    <- temp$x
        y    <- temp$y
        polygon(x, y, border=col[i], col="white")
      }

      val <- or[i, j]
      if (abs(val - 1) >= 0.001) {
        # Define the rectangle
        temp <- definePoly(x0, val)
        x    <- temp$x
        y    <- temp$y
        polygon(x, y, border=NA, col=col[i])
      } else {
        x <- c(x0, x0+1)
        y <- c(1,     1)
        lines(x, y, col=col[i])
      }

      # Update x0
      x0 <- x0 + 1 + add
    }
    # Update x0
    x0 <- x0 + add2
 
    xavg[j] <- (a + b)/2
  }

  # x-axis
  xlabs  <- op[["xticks", exact=TRUE]]
  if (is.null(xlabs)) {
    # Use or matrix names
    xlabs <- colnames(or)
    if (is.null(xlabs)) xlabs <- FALSE
  }

  axis(1, at=xavg, tick=FALSE, labels=xlabs, font=2)

  # y-axis tick marks
  at <- op[["yticks", exact=TRUE]]
  if (is.null(at)) {
    # Determine step size
    if (ylimFlag) {
      min20 <- abs(ylim[1]-1)
      min21 <- abs(ylim[2]-1)
      h     <- max(min20, min21)/4 
      min20 <- ylim[1]
      min21 <- ylim[2]
    } else {
      min20 <- min(abs(minor-1), abs(ylim[1]-1))
      min21 <- min(abs(maxor-1), abs(ylim[2]-1))
      h     <- max(min20, min21)/4 
      min20 <- min(minor, ylim[1])
      min21 <- min(maxor, ylim[2])
    }
    at    <- 1 
    for (i in 1:100) {
      val <- 1 + i*h
      if (val > min21) {
        break
      } else {
        at <- c(at, val)
      }
    }
    for (i in 1:100) {
      val <- 1 - i*h
      if (val < min20) {
        break
      } else {
        at <- c(at, val)
      }
    }
    at <- round(at, digits=op$digits)
  }

  axis(2, at=at, font=2)
  box()

  title <- op[["title", exact=TRUE]]
  if (!is.null(title)) {
    if (is.character(title)) {
      main <- title
      sub  <- NULL
    } else {
      main <- title[["main", exact=TRUE]]
      sub  <- title[["sub", exact=TRUE]]
    }
    title(main=main, sub=sub, font.sub=2)
  }

  # Legend
  if (!legendFlag) return(NULL)

  # Get legend text
  str <- legend[["legend", exact=TRUE]]
  if (is.null(str)) {
    str <- rownames(or)
    if (is.null(str)) return(NULL)  
  }
  legend(x=legend$x, y=NULL, legend=str, bty=legend$bty, title=legend$title, 
           cex=legend$cex, fill=col[1:d[1]], horiz=legend$horiz, border=col[1:d[1]],
           inset=legend$inset, y.intersp=legend$y.intersp, x.intersp=legend$x.intersp)
  
  NULL

} # END: OR.plot.main

# Function to plot ORs from snp.effects object of a particular method
OR.plot.type.method <- function(obj, type, method, op=NULL) {

  # obj    Object of class snp.effects or snp.effects.method
  ##############################################################
  # op
  #  levels1
  #  levels2

  op <- default.list(op, c("addCI"), list(0))
  addCI <- op$addCI

  method <- toupper(method)
  cls <- class(obj)
  if ("snp.effects" %in% cls) {
    obj <- obj[[method, exact=TRUE]]
  } 
  if (is.null(obj)) stop("ERROR in OR.plot.type.method: with obj") 

  # Determine the size of the OR matrix
  temp <- obj[[type, exact=TRUE]]
  or <- temp$effects
  if (addCI) {
    lower <- temp$lower95
    upper <- temp$upper95
  } else {
    lower <- NULL
    upper <- NULL
  }
  if (is.null(or)) stop("ERROR in OR.plot.type.method: with effects")
  var1 <- attr(temp, "var1", exact=TRUE)
  var2 <- attr(temp, "var2", exact=TRUE)
  lev1 <- attr(temp, "levels1", exact=TRUE)
  lev2 <- attr(temp, "levels2", exact=TRUE)

  # Check the levels
  d <- dim(or)
  if (is.null(d)) {
    len     <- length(or)
    dim(or) <- c(1, len)
    d <- dim(or)    
    if (addCI) {
      dim(lower) <- c(1, len)
      dim(upper) <- c(1, len)
    }
  }
  levels1 <- op[["levels1", exact=TRUE]]
  if (is.null(levels1)) levels1 <- lev1
  levels2 <- op[["levels2", exact=TRUE]]
  if (is.null(levels2)) levels2 <- lev2

  rows <- 1:d[1]
  cols <- 1:d[2]

  temp <- lev1 %in% levels1
  if (any(temp)) {
    lev1 <- lev1[temp]
    rows <- rows[temp]
  }
  temp <- lev2 %in% levels2
  if (any(temp)) {
    lev2 <- lev2[temp]
    cols <- cols[temp]
  }
  or <- or[rows, cols]
  if (addCI) {
    lower <- lower[rows, cols]
    upper <- upper[rows, cols]
  }
  if (is.null(dim(or))) {
    len     <- length(or)
    dim(or) <- c(1, len)
    if (addCI) {
      dim(lower) <- c(1, len)
      dim(upper) <- c(1, len)
    }
  }

  # Get the options
  op <- OR.getTitle(2, op, type, var1, var2, levels1=lev1, method=method, levels2=lev2)

  legend <- op[["legend", exact=TRUE]]
  flag   <- 0 
  if (is.null(legend)) {
    if (nrow(or) == 1) {
      legend <- NA
      flag   <- 1
    } else {
      legend <- list(title=var1, legend=lev1)
    }
  }
  if (!flag) legend <- default.list(legend, c("title", "legend"), list(var1, lev1))
  xticks <- lev2

  # Watch out for 1 level2
  if (length(lev2) == 1) {
    legend <- NA
    xticks <- lev1
  }

  op$legend <- legend
  op$xticks <- xticks

  # Call the main function
  OR.plot.main(or, lower=lower, upper=upper, op=op)

  NULL

} # END: OR.plot.type.method

# Function to plot ORs from snp.effects object
OR.plot.type <- function(obj, type, op=NULL) {

  # op
  #  method
  #  levels1    Single level
  #             The default is 1
  #  levels2    Vector of levels or NULL
  #             The default is NULL 

  allMethods <- c("UML", "CML", "EB", "HCL", "CCL", "CLR")
  methods <- op[["method", exact=TRUE]]
  if (is.null(methods)) {
    methods <- allMethods
  } else {
    methods <- toupper(methods)
    temp <- methods %in% allMethods
    methods <- methods[temp]
  }
  if (!length(methods)) stop("Incorrect methods selected")
 
  temp <- methods %in% names(obj)
  methods <- methods[temp]
  n <- length(methods)
  if (!n) stop("ERROR with op$method and/or obj")

  if (n == 1) return(OR.plot.type.method(obj, type, methods, op=op))

  op <- default.list(op, c("levels1", "addCI"), list(1, 0))
  addCI <- op$addCI

  # Determine the size of the OR matrix
  mat <- obj[[1]][[type]]$effects
  if (is.null(mat)) stop("ERROR in OR.plot.type")

  # Get the levels
  levels1 <- op$levels1
  if (length(levels) > 1) stop("ERROR: op$levels1 must be of length 1")
  levels2 <- op[["levels2", exact=TRUE]]

  # Get the matrix of ORs
  for (i in 1:n) {
    temp <- obj[[methods[i]]][[type]]
    eff  <- temp$effects
    if (addCI) {
      l95 <- temp$lower95
      u95 <- temp$upper95
    }
    if (i == 1) {
      var1 <- attr(temp, "var1", exact=TRUE)
      var2 <- attr(temp, "var2", exact=TRUE)
      lev1 <- attr(temp, "levels1", exact=TRUE)
      lev2 <- attr(temp, "levels2", exact=TRUE)
      nlev1 <- length(lev1)
      nlev2 <- length(lev2)
      if (!(levels1 %in% lev1)) stop("ERROR: with op$levels1")
      temp1 <- lev1 == levels1
      cols  <- (1:nlev2)
      if (!is.null(levels2)) {
        temp2 <- lev2 %in% levels2
        if (any(temp2)) {
          lev2 <- lev2[temp2]
          cols <- cols[temp2]
        }
      }
      nc    <- length(lev2)
      or    <- matrix(data=NA, nrow=n, ncol=nc)
      rows  <- (1:nlev1)[temp1]

      if (addCI) {
        lower <- matrix(data=NA, nrow=n, ncol=nc)
        upper <- matrix(data=NA, nrow=n, ncol=nc)
      } else {
        lower <- NULL
        upper <- NULL
      }
    }
    or[i, ] <- eff[rows, cols]
    if (addCI) {
      lower[i, ] <- l95[rows, cols]
      upper[i, ] <- u95[rows, cols]
    } 
  }

  # Get the options
  op <- OR.getTitle(1, op, type, var1, var2, levels1=levels1)

  legend <- op[["legend", exact=TRUE]]
  flag   <- 0 
  if (is.null(legend)) {
    if (nrow(or) == 1) {
      legend <- NA
      flag   <- 1
    } else {
      legend <- list(legend=methods)
    }
  }
  if (!flag) legend <- default.list(legend, c("legend"), list(methods))
  op$legend <- legend
  op$xticks <- lev2

  # Call the main function
  OR.plot.main(or, lower=lower, upper=upper, op=op)


} # END: OR.plot.type

# Function to return the title
OR.getTitle <- function(which, op, type, var1, var2, levels1=NULL, method=NULL, levels2=NULL) {

  n1 <- length(levels1)
  n2 <- length(levels2)
  title <- op[["title", exact=TRUE]]
  if (is.null(title)) title <- list()
  if (is.character(title)) title <- list(main=title)
  mainFlag <- !is.null(title[["main", exact=TRUE]])
  subFlag  <- !is.null(title[["sub", exact=TRUE]])
  if (!mainFlag) {
    allTypes <- c("JointEffects", "StratEffects", "StratEffects.2")
    temp     <- c("Joint Effects", "Stratified Effects", "Stratified Effects")  
    i        <- match(type, allTypes)
    if (which == 1) {
      str   <- paste(temp[i], " of ", var1, "=", levels1, " and ", var2, sep="")
    } else {
      str <- paste(method, " ", temp[i], " of ", var1, sep="")
      if (n1 == 1) str <- paste(str, "=", levels1, sep="")
      str <- paste(str, " and ", var2, sep="")
      if (n2 == 1) str <- paste(str, "=", levels2, sep="")
    }
    title$main <- str
  }
  if (!subFlag) {
    if ((n1 > 1) && (n2 == 1)) {
      title$sub <- var1
    } else {
      title$sub <- var2
    }
  }
  op$title <- title
  
  op

} # END: OR.getTitle 

# Function for an OR plot 
snp.effects.plot <- function(obj.list, op=NULL) {

  # obj.list     Return object from snp.effects or list of return objects from snp.effects
  #              No default
  ########################################################################################
  # NOTE: if obj.list contains more than 1 object, then the options list can
  #       be a list of sublists of the following options.
  # op           List of options
  #  method      Character vector of "UML", "CML", "EB", "HCL", "CCL", "CLR"
  #              The default is NULL
  #  type        One of "JointEffects", "StratEffects", "StratEffects.2"
  #              The default is "StratEffects"
  #  split.screen NULL or 2-element vector to split the plot screen
  #               The default is NULL
  #  levels1
  #  levels2
  #  legend      List with names
  #              The default is NULL
  #  ylim        Numeric vector of length 2 to specify the y limits to be used
  #              for all plots
  #              The default is NULL
  #  addCI       0 or 1 to add 95% confidence intervals to the plot
  #              The default is 0.
  #  title       Character string or list with names "main" and "sub" for the
  #              main and sub titles.
  #              The default is NULL.
 
  
  allMethods <- c("UML", "CML", "EB", "HCL", "CCL", "CLR")
  allTypes   <- c("JointEffects", "StratEffects", "StratEffects.2")

  if (is.null(op)) op <- list(list())
  temp <- c("method", "type", "split.screen", "levels1", "levels2", "ylim", 
            "addCI", "title")
  if (any(names(op) %in% temp)) op <- list(op)
  nop <- length(op)
  for (i in 1:nop) {
    op[[i]] <- default.list(op[[i]], c("type", "method", "addCI"), list("StratEffects", NULL, 0),
            checkList=list(allTypes, allMethods, 0:1))
    type <- op[[i]]$type
    if (length(type) > 1) stop("ERROR: op$type must have length 1")
  }
  
  if (class(obj.list) == "snp.effects") obj.list <- list(obj.list)
  n <- length(obj.list)
  temp <- op[[1]]
  svec <- temp[["split.screen", exact=TRUE]]
  if (is.null(svec)) {
    nrow <- floor(n/2)
    if (!nrow) nrow <- 1
    if (nrow == 1) {
      ncol <- n
    } else {
      for (i in 2:100) {
        if (i*nrow >= n) {
          ncol <- i
          break
        }
      }
    }
    svec <- c(nrow, ncol) 
  }  
  if (length(svec) != 2) stop("ERROR with op$split.screen")
  
  split.screen(svec)
  j <- 1
  for (i in 1:n) {
    screen(i)
    temp <- try(OR.plot.type(obj.list[[i]], op[[j]]$type, op=op[[j]]), silent=TRUE)
    close.screen(i)
    j <- j + 1
    if (j > nop) j <- 1
  }

  invisible(NULL)

} # END: snp.effects.plot

# Function in ARTP package to call gene.plot
plot_genes <- function(obs.outfile, gene.list, op=NULL) {

  op <- default.list(op, c("cex.axis", "colors", "maxLabelLen", "chr.text", "x.las", "x.padj"),
                        list(0.75, NULL, NULL, NULL, 2, NULL))
  gene.list <- default.list(gene.list, c("snp.var", "gene.var", "chrm.var"), 
                 list("SNP", "Gene", "Chr"))
  temp <- list(add.CI=1, gene.list=gene.list)
  x <- convert.obsfile(obs.outfile, op=temp)

  # Check for the chromosome var 
  chr <- gene.list$chrm.var
  if (chr %in% colnames(x)) {
    flag <- 1
    chr.text <- op[["chr.text", exact=TRUE]]
  } else {
    flag <- 0
    chr.text <- 0
    chr <- "chr"
    x[, chr] <- "1"
  } 
  x[, "loc"] <- 1
  op$chr.text <- chr.text

  lmap.list <- list(snp.var="SNP", gene.var=gene.list$gene, chrm.var=chr, loc.var="loc")
  ret <- gene.plot(x, "Pvalue", lmap.list, op=op)

  vars <- "loc"
  if (!flag) vars <- c(vars, chr)
  x <- removeOrKeepCols(x, vars, which=-1)

  x

} # END: plot_genes


