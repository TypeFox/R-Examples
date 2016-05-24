###########################################################################/**
# @set "class=AffymetrixCelFile"
# @RdocMethod plotDensity
#
# @title "Plots the density of the probe signals on the array"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{subset}{The subset of probes to considered \emph{before} any
#     filtering by probe type is applied.
#     If a @vector of @doubles, the cell indices.
#     If a scalar @double in [0,1], the fraction of cells, which can
#     be used to speed up the plotting if approximate densities are
#     acceptable.
#     if @NULL, all cells are considered.
#   }
#   \item{types}{The type of probes to include, e.g. \code{"all"},
#     \code{"pmmm"}, \code{"pm"}, and \code{"mm"}.}
#   \item{...}{Additional arguments passed to
#              @see "aroma.light::plotDensity.numeric".}
#   \item{xlim}{The range on the x axis.}
#   \item{xlab,ylab}{The labels on the x and the y axes.}
#   \item{log}{If @TRUE, the density of the log (base 2) values are
#      used, otherwise the non-logged values.}
#   \item{verbose}{A @logical or a @see "R.utils::Verbose" object.}
# }
#
# \value{
#  Returns nothing.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("plotDensity", "AffymetrixCelFile", function(this, subset=NULL, types=NULL, ..., xlim=c(0,16), xlab=NULL, ylab="density (integrates to one)", log=TRUE, annotate=TRUE, verbose=FALSE) {
  ## aroma.light::plotDensity()
  requireNamespace("aroma.light") || throw("Package aroma.light not loaded.")


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'subset':

  # Argument 'xlab':
  if (is.null(xlab)) {
    if (log) {
      xlab <- expression(log[2](y));
    } else {
      xlab <- expression(y);
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify the subset of probes to be updated
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);
  verbose && enter(verbose, "Identifying subset of probes");
  suppressWarnings({
    subset <- identifyCells(cdf, indices=subset, types=types,
                                                    verbose=less(verbose));
  })
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Plot density
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Plotting the density");
  verbose && cat(verbose, "Array: ", getName(this));
  suppressWarnings({
    verbose && enter(verbose, "Loading probe intensities");
    y <- getData(this, indices=subset, fields="intensities");
    y <- y$intensities;
    verbose && exit(verbose);
    if (log) {
      verbose && cat(verbose, "Taking the logarithm (base 2)");
      y <- log(y, base=2);
    }
    verbose && cat(verbose, "Plotting");
    plotDensity(y, xlim=xlim, xlab=xlab, ylab=ylab, ...);
  })

  if (annotate) {
    stextChipType(getChipType(this));
    stextLabels(this);
    stextSize(this, size=length(y));
  }

  verbose && exit(verbose);
})



setMethodS3("getAm", "AffymetrixCelFile", function(this, reference, indices=NULL, ..., zeros=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'reference':
  reference <- Arguments$getInstanceOf(reference, "AffymetrixCelFile");

  # Argument 'indices':
  nbrOfCells <- nbrOfCells(this);
  if (is.null(indices)) {
  } else {
    indices <- Arguments$getIndices(indices, max=nbrOfCells);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Further validation
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check if the two CEL files are compatible
  if (nbrOfCells != nbrOfCells(reference)) {
    throw("This and the 'reference' CEL file have different number of cells: ",
                                   nbrOfCells, " != ", nbrOfCells(reference));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the signals for this channel
  y1 <- getData(this, indices=indices, fields="intensities")[,1];


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Offset signals?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  offset <- this$offset;
  if (is.null(offset))
    offset <- 0;
  if (offset != 0)
    cat("Offset: ", offset, "\n", sep="");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Remove signals that are zero?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!zeros) {
    keep <- which(y1 != 0);
    y1 <- y1[keep];
  } else {
    keep <- seq_along(y1);
  }
  y1 <- y1 + offset;
  y1 <- log(y1, base=2);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get reference signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (length(y1) == 0) {
    y2 <- y1;
  } else {
    # Get the signals for the reference channel
    if (is.null(indices)) {
      indices <- keep;
    } else {
      indices <- indices[keep];
    }
    y2 <- getData(reference, indices=indices, fields="intensities")[,1];
    y2 <- y2 + offset;
    y2 <- log(y2, base=2);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return (A,M)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  am <- matrix(c((y1+y2)/2, y1-y2), ncol=2);
  colnames(am) <- c("A", "M");

  am;
})



setMethodS3("annotateMvsA", "AffymetrixCelFile", function(this, reference, ..., what="M") {
  if (identical(what, "M")) {
    abline(h=0, lty=1, col="blue");
  }
  stextChipType(getChipType(this));
  stextLabels(this, others=reference);
}, private=TRUE)




###########################################################################/**
# @RdocMethod plotMvsA
#
# @title "Plots log-ratio versus log-intensity in a scatter plot"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{reference}{The reference channel, i.e. the denominator in the
#     log ratios.}
#   \item{indices}{Indices of the probe signals to be plotted.}
#   \item{pch}{The plot symbol.}
#   \item{xlim,ylim}{The range of the x and the y axes.}
#   \item{xlab,ylab}{The labels on the x and the y axes.}
#   \item{...}{Additional arguments passed to @see "graphics::plot".}
#   \item{annotate}{If @TRUE, the plot is annotated with information about
#     the data plotted, otherwise not.}
# }
#
# \value{
#  Returns (invisibly) a @data.frame with the plotted data in the
#  first two columns.
# }
#
# @author "HB"
#
# \seealso{
#   @seemethod "smoothScatterMvsA".
#   @seemethod "plotMvsX".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("plotMvsA", "AffymetrixCelFile", function(this, reference, indices=NULL, pch=176, xlim=c(0,16), ylim=c(-1,1)*diff(xlim), xlab=expression(A==1/2%*%log[2](y[1]*y[2])), ylab=expression(M==log[2](y[1]/y[2])), ..., annotate=TRUE) {
  ma <- getAm(this, reference, indices=indices);
  plot(ma, pch=pch, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...);
  if (annotate) {
    annotateMvsA(this, reference);
    stextSize(this, size=nrow(ma));
  }
  this$lastPlotData <- ma;
  invisible(ma);
})



###########################################################################/**
# @RdocMethod smoothScatterMvsA
#
# @title "Plots log-ratio versus log-intensity in a smooth scatter plot"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{reference}{The reference channel, i.e. the denominator in the
#     log ratios.}
#   \item{indices}{Indices of the probe signals to be plotted.}
#   \item{pch}{The plot symbol.}
#   \item{xlim,ylim}{The range of the x and the y axes.}
#   \item{xlab,ylab}{The labels on the x and the y axes.}
#   \item{...}{Additional arguments passed to @see "graphics::plot".}
#   \item{annotate}{If @TRUE, the plot is annotated with information about
#     the data plotted, otherwise not.}
# }
#
# \value{
#  Returns (invisibly) a @data.frame with the plotted data in the
#  first two columns.
# }
#
# @author "HB"
#
# \seealso{
#   @seemethod "plotMvsA".
#   Internally @see "graphics::smoothScatter" is used.
#   @seeclass
# }
#*/###########################################################################
setMethodS3("smoothScatterMvsA", "AffymetrixCelFile", function(this, reference, indices=NULL, pch=176, xlim=c(0,16), ylim=c(-1,1)*diff(xlim), xlab=expression(A==1/2%*%log[2](y[1]*y[2])), ylab=expression(M==log[2](y[1]/y[2])), ..., annotate=TRUE) {
  ma <- getAm(this, reference, indices=indices);
  smoothScatter(ma, pch=pch, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...);
  if (annotate) {
    annotateMvsA(this, reference);
    stextSize(this, size=nrow(ma));
  }
  this$lastPlotData <- ma;
  invisible(ma);
})




###########################################################################/**
# @RdocMethod plotMvsX
#
# @title "Plots log-ratio versus another variable in a scatter plot"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{reference}{The reference channel, i.e. the denominator in the
#     log ratios.}
#   \item{x}{The other variable.  A @double @vector.}
#   \item{indices}{Indices of the probe signals to be plotted.}
#   \item{pch}{The plot symbol.}
#   \item{xlim,ylim}{The range of the x and the y axes.}
#   \item{xlab,ylab}{The labels on the x and the y axes.}
#   \item{...}{Additional arguments passed to @see "graphics::plot".}
#   \item{annotate}{If @TRUE, the plot is annotated with information about
#     the data plotted, otherwise not.}
# }
#
# \value{
#  Returns (invisibly) a @data.frame with the plotted data in the
#  first two columns, and remaining data in the following columns.
# }
#
# @author "HB"
#
# \seealso{
#   @seemethod "plotMvsA".
#   @seemethod "smoothScatterMvsA".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("plotMvsX", "AffymetrixCelFile", function(this, reference, x, indices=NULL, pch=176, ylim=c(-1,1)*2, ylab=NULL, ..., what=c("M", "A"), add=FALSE, annotate=!add) {
  # Argument 'what':
  what <- match.arg(what);

  # Get the log-ratios
  ma <- getAm(this, reference, indices=indices, zeros=TRUE);
  nobs <- nrow(ma);
  if (nobs == 0)
    throw("Cannot plot M vs X because there is not non-zero data.");

  if (nobs != length(x)) {
    throw("The number of log-ratios does not match the number of elements in argument 'x': ", length(nobs), " != ", length(x));
  }

  if (what == "M") {
    ylab <- expression(M==log[2](y1/y2))
  } else {
    ma <- ma[,2:1];
    ylab <- expression(A==1/2%*%log[2](y1*y2))
  }

  if (add) {
    points(x, ma[,1], pch=pch, ...);
  } else {
    plot(x, ma[,1], pch=pch, ylim=ylim, ylab=ylab, ...);
    if (annotate) {
      annotateMvsA(this, reference, what=what);
      stextSize(this, size=length(x));
    }
  }

  # The first two columns should always be the data plotted
  ma <- cbind(x=x, ma);

  this$lastPlotData <- ma;
  invisible(ma);
})





setMethodS3("highlight", "AffymetrixCelFile", function(this, indices=NULL, ...) {
  data <- this$lastPlotData;
  if (!is.null(indices))
    data <- data[indices,,drop=FALSE];
  points(data[,1:2], ...);
  invisible(data);
}, protected=TRUE)


###########################################################################/**
# @RdocMethod image270
#
# @title "Displays all or a subset of the data spatially"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{xrange}{A @numeric @vector of length two giving the left and right
#          coordinates of the cells to be returned.}
#   \item{yrange}{A @numeric @vector of length two giving the top and bottom
#          coordinates of the cells to be returned.}
#   \item{...}{Additional arguments passed @see "graphics::image" and [...].}
#   \item{field}{The data field to be displayed.}
#   \item{col}{The color map to be used.}
#   \item{main}{The main title of the plot.}
# }
#
# \value{
#  Returns the (270-degrees rotated) data @matrix.
# }
#
# @author "HB"
#
# \seealso{
#   @seemethod "updateUnits".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("image270", "AffymetrixCelFile", function(this, xrange=c(0,Inf), yrange=c(0,Inf), takeLog=TRUE, interleaved=FALSE, ..., field=c("intensities", "stdvs", "pixels"), col=gray.colors(256), main=getName(this)) {
  rotate270 <- function(x, ...) {
    x <- t(x)
    nc <- ncol(x)
    if (nc < 2) return(x)
    x[,nc:1,drop=FALSE]
  }

  # Argument 'field':
  field <- match.arg(field);

  suppressWarnings({
    y <- readRawDataRectangle(this, xrange=xrange, yrange=yrange,
                                              fields=field, ..., drop=TRUE);
  })

  # if only PM locations have signal, add a fake row

  nr <- nrow(y);
  if (interleaved) {
    idxEven <- which((1:nr) %% 2 == 0);
    y[idxEven-1,] <- y[idxEven,];
  }

  suppressWarnings({
    if (takeLog) {
      image(log2(rotate270(y)), col=col, ..., axes=FALSE, main=main);
    } else {
      image(rotate270(y), col=col, ..., axes=FALSE, main=main);
    }
  })

  if (is.null(xrange) || xrange[2] == Inf)
    xrange <- c(0,ncol(y)-1);
  if (is.null(yrange) || yrange[2] == Inf)
    yrange <- c(0,nrow(y)-1);

  cdf <- getCdf(this);
  dim <- paste(getDimension(cdf), collapse="x");
  label <- sprintf("Chip type: %s [%s]", getChipType(this), dim);
  text(x=0, y=0, labels=label, adj=c(0,1.2), cex=0.8, xpd=TRUE)
  label <- sprintf("(%d,%d)", as.integer(xrange[1]), as.integer(yrange[1]));
  text(x=0, y=1, labels=label, adj=c(0,-0.7), cex=0.8, xpd=TRUE)
  label <- sprintf("(%d,%d)", as.integer(xrange[2]), as.integer(yrange[2]));
  text(x=1, y=0, labels=label, adj=c(1,1.2), cex=0.8, xpd=TRUE)

  # Return the plotted data.
  invisible(y);
})




###########################################################################/**
# @RdocMethod getImage
#
# @title "Creates an RGB Image object from a CEL file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{other}{An optional @see "AffymetrixCelFile" of the same chip type,
#      that is used for calculating the ratio (non-logged).  Note, to get
#      the log-ratios, the \code{log}() function has to be specified as
#      the first transform in the @list of \code{transformations}.}
#   \item{xrange, yrange}{@vectors of length two specifying the
#      (x0,x1) and (y0,y1) regions to be extracted.  If @NULL, the
#      complete regions is used.}
#   \item{field}{One of the CEL file fields, i.e. \code{"intensities"},
#      \code{stdvs}, or \code{pixels}.}
#   \item{zoom}{A @numeric scale factor in (0,+Inf) for resizing the
#     imaging. If \code{1}, no resizing is done.}
#   \item{palette}{An optional @vector of color code.}
#   \item{...}{Additional arguments passed to
#      @seemethod "readRawDataRectangle" and more function.}
#   \item{readRectFcn}{A @function taking arguments 'xrange' and 'yrange',
#     or @NULL for the default read function.}
#   \item{verbose}{A @logical or a @see "R.utils::Verbose" object.}
# }
#
# \value{
#   Returns an Image object as defined by the EBImage package.
#   If \code{palette==NULL}, the color code is \code{Grayscale}, otherwise
#   \code{TrueColor}.
# }
#
# @author "KS, HB"
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getImage", "AffymetrixCelFile", function(this, other=NULL, transforms=list(sqrt), xrange=c(0,Inf), yrange=xrange, zrange=c(0,sqrt(2^16)), field=c("intensities", "stdvs", "pixels"), zoom=1, ..., readRectFcn=NULL, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  readRectangleByField <- function(this, other=NULL, xrange, yrange, ...) {
    suppressWarnings({
      y <- readRawDataRectangle(this, xrange=xrange, yrange=yrange,
                                               fields=field, ..., drop=TRUE);
    });

    if (is.null(other)) {
    } else {
      if (inherits(other, "AffymetrixCelFile")) {
        suppressWarnings({
          yR <- readRawDataRectangle(other, xrange=xrange, yrange=yrange,
                                               fields=field, ..., drop=TRUE);
        });
      } else {
        yR <- other;
      }

      y <- y/yR;
    }

    y;
  } # readRectangleByField()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'other':
  if (!is.null(other)) {
    if (inherits(other, "AffymetrixCelFile")) {
      hdr1 <- getHeader(this);
      hdr2 <- getHeader(other);
      fields <- c("rows", "cols");
      if (!identical(hdr1[fields], hdr2[fields])) {
        throw("Argument 'other' contains an ", class(other)[1], " with a dimension not compatible with the main ", class(this)[1], "");
      }
    } else {
      throw("Argument 'other' is of an unknown class: ", other);
    }
  }

  # Argument 'field':
  field <- match.arg(field);

  # Argument 'zoom':
  zoom <- Arguments$getDouble(zoom, range=c(0,Inf));

  # Argument 'readRectFcn':
  if (is.null(readRectFcn)) {
    readRectFcn <- readRectangleByField;
  } else if (!is.function(readRectFcn)) {
    throw("Argument 'readRectFcn' is not a function: ", readRectFcn);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Getting CEL image");

  verbose && enter(verbose, "Reading CEL image");
  y <- readRectFcn(this, other=other, xrange=xrange, yrange=yrange, ...);
  verbose && str(verbose, y);
  verbose && summary(verbose, as.vector(y[is.finite(y) & (y != 0)]));
  verbose && printf(verbose, "RAM: %.1fMB\n", object.size(y)/1024^2);
  verbose && exit(verbose);

  verbose && enter(verbose, "Creating Image");
  img <- getImage(y, transforms=transforms, scale=zoom, lim=zrange, ...,
                                                   verbose=less(verbose, 1));
  verbose && print(verbose, img);
  verbose && exit(verbose);

  verbose && exit(verbose);

  # Return the 'field'
  attr(img, "field") <- field;

  img;
})



###########################################################################/**
# @RdocMethod plotImage
#
# @title "Displays a spatial plot of a CEL file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @seemethod "getImage".}
# }
#
# \value{
#   Returns (invisibly) the image displayed.
# }
#
# @author "KS, HB"
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("plotImage", "AffymetrixCelFile", function(this, ...) {
  # Get the image
  img <- getImage(this, ...);

  # Display image
  display(img);  # Using display() for Image of aroma.core

  invisible(img);
})


###########################################################################/**
# @RdocMethod writeImage
#
# @title "Writes a spatial image of the signals in the CEL file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{A @character string specifying the filename of
#     the output file.}
#   \item{fullname}{A @character string specifying the full name of
#     the output file.}
#   \item{tags}{A @character @vector of optional tags added to the
#     already existing tags of the CEL file.}
#   \item{imgFormat}{A @character string specifying the filename extension
#     which also defines the image file format.}
#   \item{path}{The path where the image file is stored.}
#   \item{...}{Arguments passed to @seemethod "getImage".}
#   \item{verbose}{A @logical or a @see "R.utils::Verbose" object.}
# }
#
# \value{
#   Returns the pathname to the image file created.
# }
#
# \examples{\dontrun{
#   yellow.colors <- function(n) { hsv(h=0.15, v=0:(n-1)/(n-1)) }
#   df <- ds[[1]]
#   writeImage(df, tags="gray", palette=gray.colors(256), xrange=c(0,200))
#   writeImage(df, tags="yellow", palette=yellow.colors(256), xrange=c(0,200))
#   writeImage(df, tags="heat", palette=heat.colors(256), xrange=c(0,200))
#   writeImage(df, tags="terrain", palette=terrain.colors(256), xrange=c(0,200))
#   writeImage(df, tags="topo", palette=topo.colors(256), xrange=c(0,200))
#   writeImage(df, tags="cm", palette=cm.colors(256), xrange=c(0,200))
#   writeImage(df, tags="rainbow", palette=rainbow(256), xrange=c(0,200))
# }}
#
# @author "KS, HB"
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("writeImage", "AffymetrixCelFile", function(this, filename=NULL, fullname=NULL, tags=c("*", "sqrt", "gray"), imgFormat="png", path=NULL,  field=c("intensities", "stdvs", "pixels"), ..., skip=TRUE, verbose=FALSE) {
  # Argument 'path':
  if (is.null(path)) {
    rootPath <- "reports";

    # Infer the data set name and the tags from the path
    path <- getPath(this);
    parent <- getParent(path); # chip type
    parent <- getParent(path); # data set
    parts <- unlist(strsplit(basename(parent), split=","));
    dataSet <- parts[1];
    dataSetTags <- parts[-1];
    if (length(dataSetTags) == 0) {
      dataSetTags <- "raw";
    } else {
      dataSetTags <- paste(dataSetTags, collapse=",");
    }

    # chip type
    chipType <- getChipType(this, fullname=FALSE);

    # image set
    set <- "spatial";

    path <- filePath(rootPath, dataSet, dataSetTags, chipType, set);
  }
  path <- Arguments$getWritablePath(path);

  # Argument 'tags':
  tags <- Arguments$getCharacters(tags);

  # Argument 'filename':

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # Update asterisk tags?
  if ("*" %in% tags) {
    idx <- match("*", tags);
    tags[idx] <- field;
    tags <- locallyUnique(tags);
  }

  verbose && enter(verbose, "Writing CEL image to file");

  # Generate the pathname
  if (is.null(fullname)) {
    fullname <- getFullName(this);
  }
  fullname <- paste(c(fullname, tags), collapse=",");
  if (is.null(filename)) {
    filename <- sprintf("%s.%s", fullname, imgFormat);
  }
  pathname <- Arguments$getWritablePathname(filename, path=path);
  verbose && cat(verbose, "Pathname: ", pathname);

  if (!skip || !isFile(pathname)) {
    verbose && enter(verbose, "Getting image");
    img <- getImage(this, ..., verbose=less(verbose));

    verbose && cat(verbose, "Image object:");
    verbose && print(verbose, img);
    verbose && exit(verbose);

    verbose && enter(verbose, "Writing image");
    writeImage(img, file=pathname);
    verbose && exit(verbose);
  }

  verbose && exit(verbose);

  # Return pathname
  pathname;
})



############################################################################
# HISTORY:
# 2012-11-14
# o CLEANUP: writeImage() for AffymetrixCelFile no longer supports
#   sample name aliases.
# 2011-01-30
# o CLEAN UP: getImage(), plotImage() and writeImage() for AffymetrixCelFile
#   no longer depend explicitly on EBImage but instead calls aroma.core's
#   getImage(), display(), and writeImage(), respectively.  The latter
#   methods depend on EBimage for now.
# 2009-09-04
# o Now smoothScatter() is loaded via aroma.core.
# 2009-09-17
# o Now argument 'subset' of plotDensity() of AffymetrixCelFile defaults
#   to NULL (all probes).  Before it was 1/2 (a fraction).
# 2008-03-14
# o Added argument 'readRectFcn' to getImage() allowing one to read data
#   in different ways, e.g. by also read a reference array and return
#   log-ratios.
# 2007-09-17
# o Added write.image()/writeImage() caller to avoid recent EBImage
#   warnings about write.image() being deprecated.
# 2007-06-11
# o Explicit calls to geneplotter::smoothScatter() & EBImage::write.image().
# 2007-06-07
# o BUG FIX: When argument 'transforms' to getImage() of AffymetrixCelFile
#   wasn't a list, then "Error: argument "transform" is missing, with no
#   default" was thrown.  Thanks Karen Vranizan, UC Berkeley for reporting
#   this problem.
# 2007-03-20
# o Now writeImage() uses the file alias as the name if it exists.
# 2007-03-19
# o Now getImage() and writeSpatial() accepts a list of transform functions.
# 2007-03-15
# o Argument '...' to plotDensity() of AffymetrixCelFile and
#   AffymetrixCelSet are no longer passed to identifyCells().
# 2007-02-16
# o Increased the threshold to detect empty rows/columns in getImage().
# 2007-02-06
# o Now writeImage() writes image files to
#   <rootPath>/<dataSet>/<tags>/<chipType>/<set>/.
# 2007-02-03
# o plotDensity() now make sure aroma.light is loaded.
# 2007-01-30
# o Image functions tested with EBImage v1.9.23 on WinXP.
# o Changed the default transform to sqrt().
# o Added 'tags' to writeImage().
# o Renamed plotImageToFile() to writeImage().
# 2007-01-12 /KS
# o TODO? Image methods needed for AffymetrixCelSets?  Or just use lapply?
# o TODO: Annotation of plots generated by EBImage. (are tags enough? /HB)
# o Added getImage(), plotImage() and plotImageToFile(), which use
#   EBImage functionality.
# o Moved image270() and writeSpatial() from AffymetrixCelFile.R.
# 2006-09-26
# o Renamed calcMvsA() to getAm().
# 2006-09-15
# o Added more Rdoc comments.
# o Readded plotDensity().
# o Added stextSize() to annotate with "n=1034".
# 2006-08-27
# o Added plotMvsX() and plotMvsPosition().
# o Added calcMvsA(), plotMvsA(), and smoothScatterMvsA().
# 2006-07-27
# o Added argument 'verbose' to plotDensity().
# 2006-05-29
# o Added Rdoc comments.
# 2006-05-16
# o Created.
############################################################################
