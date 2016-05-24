###########################################################################/**
# @RdocClass SpatialReporter
#
# @title "The SpatialReporter class"
#
# \description{
#  @classhierarchy
#
#  A SpatialReporter generates image files of spatial representations of
#  cell signals for each of the arrays in the input set.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AffymetrixCelSetReporter".}
#   \item{reference}{An optional reference @see "AffymetrixCelFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("SpatialReporter", function(..., reference=NULL) {
  this <- extend(AffymetrixCelSetReporter(...), "SpatialReporter",
    .reference = NULL
  )

  setReference(this, reference);

  this;
})

setMethodS3("as.character", "SpatialReporter", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, paste("Name:", getName(this)));
  s <- c(s, paste("Tags:", paste(getTags(this), collapse=",")));
  s <- c(s, paste("Number of arrays:", length(this)));

  # Reference?
  refFile <- getReference(this);
  if (!is.null(refFile)) {
    s <- c(s, paste("<Relative to reference>"));
    s <- c(s, paste("Name:", getName(refFile)));
    s <- c(s, paste("Tags:", getTags(refFile, collapse=",")));
  }

  colorMaps <- getColorMaps(this);
  if (length(colorMaps) == 0)
    colorMaps <- "<no color maps; set before processing>";
  s <- c(s, paste("Color maps:", paste(colorMaps, collapse="; ")));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  GenericSummary(s);
}, protected=TRUE)


setMethodS3("getReportSet", "SpatialReporter", function(this, ...) {
  "spatial";
}, protected=TRUE)



setMethodS3("getReference", "SpatialReporter", function(this, ...) {
  this$.reference;
}, protected=TRUE)


setMethodS3("setReference", "SpatialReporter", function(this, refFile, ...) {
  # Argument 'refFile':
  if (!is.null(refFile)) {
    refFile <- Arguments$getInstanceOf(refFile, "AffymetrixCelFile");

    ds <- getDataSet(this);
    df <- getOneFile(ds);

    if (!is.element(class(refFile)[1L], class(df))) {
      throw("Cannot set reference. Argument 'refFile' is not of a class compatible with the data set: ", class(refFile)[1]);
    }
  }

  this$.reference <- refFile;
}, protected=TRUE)


setMethodS3("addColorMap", "SpatialReporter", function(this, colorMap, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'colorMap':
  colorMap <- Arguments$getCharacter(colorMap, nchar=c(1,Inf), length=c(1,1));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Parser argument 'colorMap'
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  parts <- strsplit(colorMap, split=",")[[1]];
  tags <- paste(parts, collapse=",");
  n <- length(parts);
  if (n > 3) {
    throw("Argument 'colorMap' must not contain more than three parts: ",
                                                                   tags);
  }
  if (n == 1) {
    transforms <- list();
  } else {
    transforms <- as.list(parts[-n]);
  }
  palette <- parts[n];

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup transforms
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (length(transforms) == 0) {
    transforms <- list(sqrt);
  } else {
    transforms <- lapply(transforms, FUN=function(transform) {
      # Check transform
      if (!exists(transform, mode="function")) {
        throw("Argument 'colorMap' specifies an unknown transform function ('",
                                                 transform, "'): ", colorMap);
      }
      get(transform, mode="function");
    })
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup palette
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(palette)) {
    palette <- gray.colors(256);
  } else if (is.character(palette)) {
    # Parse color palette tag
    pattern <- "^([^_]*)(|[_][0-9]*)$";
    name <- gsub(pattern, "\\1", palette);
    nbrOfColors <- gsub(pattern, "\\2", palette);
    nbrOfColors <- gsub("[_]", "", nbrOfColors);
    nbrOfColors <- as.integer(nbrOfColors);
    if (is.na(nbrOfColors))
      nbrOfColors <- 256;

    # Search for function <palette>() and then <palette>.colors()
    if (!exists(name, mode="function")) {
      name <- sprintf("%s.colors", palette);
      if (!exists(name, mode="function")) {
        throw("Argument 'colorMap' specifies an unknown palette function ('",
                                                   name, "'): ", colorMap);
      }
    }
    fcn <- get(name, mode="function");
    palette <- fcn(nbrOfColors);
  }

  map <- list(list(
    tags = tags,
    transforms = transforms,
    palette = palette
  ));
  names(map) <- map$tags;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Add color map
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  colorMaps <- this$.colorMaps;
  if (is.null(colorMaps))
    colorMaps <- list();
  colorMaps <- c(colorMaps, map);
  this$.colorMaps <- colorMaps;
})


setMethodS3("setColorMaps", "SpatialReporter", function(this, colorMaps=c("sqrt,yellow", "sqrt,rainbow"), ...) {
  this$.colorMaps <- NULL;
  for (colorMap in colorMaps) {
    addColorMap(this, colorMap, ...);
  }
})

setMethodS3("getColorMaps", "SpatialReporter", function(this, parsed=FALSE, ...) {
  colorMaps <- this$.colorMaps;
  if (!parsed) {
    colorMaps <- sapply(colorMaps, .subset2, "tags");
    colorMaps <- unlist(colorMaps);
    colorMaps <- unique(colorMaps);
    colorMaps <- sort(colorMaps);
  }

  colorMaps;
})


setMethodS3("writeImages", "SpatialReporter", function(this, arrays=NULL, aliases=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'arrays':
  cs <- getDataSet(this);
  nbrOfArrays <- length(cs);
  if (is.null(arrays)) {
    arrays <- seq_len(nbrOfArrays);
  } else {
    arrays <- Arguments$getIndices(arrays, max=nbrOfArrays);
    nbrOfArrays <- length(arrays);
  }

  # Argument 'aliases':
  if (!is.null(aliases)) {
    aliases <- Arguments$getCharacters(aliases, length=nbrOfArrays);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Get the path to the image directory
  path <- getPath(this);

  # Relative to a reference?
  refFile <- getReference(this);

  # Get the color maps to be generated
  colorMaps <- getColorMaps(this, parsed=TRUE);
  if (length(colorMaps) == 0) {
    warning("No color maps specified. Nothing to do.");
    return(invisible(path));
  }

  # For each array...
  for (kk in seq_along(arrays)) {
    df <- cs[[arrays[kk]]];

    # Aliases are deprecated
    if (!is.null(aliases)) {
      setAlias(df, aliases[kk]);
    }

    verbose && enter(verbose, sprintf("Array #%d of %d ('%s')",
                                             kk, nbrOfArrays, getName(df)));

    # For each color map...
    for (ll in seq_along(colorMaps)) {
      colorMap <- colorMaps[[ll]];
      tags <- colorMap$tags;
      verbose && enter(verbose, sprintf("Color map #%d ('%s')", ll, tags));
#      verbose && str(verbose, colorMap$transforms);
#      verbose && str(verbose, colorMap$palette);
      writeImage(df, other=refFile, path=path,
                 transforms=colorMap$transforms, palette=colorMap$palette,
                                  tags=tags, ..., verbose=less(verbose, 5));
#      gc <- gc();
      verbose && exit(verbose);
    }
    verbose && exit(verbose);
  }

  invisible(path);
}, private=TRUE)



###########################################################################/**
# @RdocMethod process
#
# @title "Generates image files, scripts and dynamic pages for the explorer"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns nothing.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("process", "SpatialReporter", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Generating ", class(this)[1], " report");

  # Generate bitmap images
  writeImages(this, ..., verbose=less(verbose));

  verbose && exit(verbose);
})


setMethodS3("readRawDataRectangle", "SpatialReporter", function(this, array, ..., field="intensities", transforms=list(), verbose=FALSE) {
  ds <- getDataSet(this);
  df <- ds[[array]];

  y <- readRawDataRectangle(df, fields=field, ..., drop=TRUE, verbose=less(verbose, 5));
##    verbose && str(verbose, y);

  # Relative signals?
  refFile <- getReference(this);
  if (!is.null(refFile)) {
    yR <- readRawDataRectangle(refFile, fields=field, ..., drop=TRUE, verbose=less(verbose, 5));
##    verbose && str(verbose, yR);
    y <- y/yR;
    # Not needed anymore
    yR <- NULL;
##    verbose && str(verbose, y);
  }

  # Transform data
  for (transform in transforms) {
    y <- transform(y);
  }

  y;
}, protected=TRUE)


setMethodS3("calculateMargins", "SpatialReporter", function(this, unshift=TRUE, ..., verbose=FALSE) {
  colMedians <- function(x, ...) {
    x <- t(x);
    rowMedians(x, ...);
  }

  # Read data
  y <- readRawDataRectangle(this, ...);

  # Remove average?
  if (unshift) {
    yAvg <- median(y, na.rm=TRUE);;
    y <- y - yAvg;
  }

  # Calculate margins
  yR <- rowMedians(y, na.rm=TRUE);
  yC <- colMedians(y, na.rm=TRUE);

  list(rowAvgs=yR, colAvgs=yC);
}, protected=TRUE)


setMethodS3("plotMargins", "SpatialReporter", function(this, array, margins=c("rows", "columns"), ..., pch=20, cex=0.7, ylim=NULL, ylab=NULL, rotate=0, verbose=FALSE) {plotMargins
  # Argument 'margins':
#  if (!all(margins %in% formals(plotMargins.SpatialReporter)$margins)) {
#    throw("Unknown value(s) in argument 'margins': ", paste(margins));
#  }

  # Get the array file
  ds <- getDataSet(this);
  df <- ds[[array]];

  yMargins <- calculateMargins(this, array=array, ..., verbose=verbose);
  keep <- is.element(c("rows", "columns"), margins);
  yMargins <- yMargins[keep];

  if (is.null(ylim)) {
    ylim <- c(NA, NA);
    for (kk in seq_along(yMargins)) {
      ylim <- range(c(ylim, range(yMargins[[1]], na.rm=TRUE)));
    }
  }

  if (is.null(ylab)) {
    ylab <- "signal";
  }

  xlabs <- margins;
  mains <- paste("Average signal per", substring(margins, 1, nchar(margins)-1));

  if (length(yMargins) > 1) {
    layout(matrix(seq_along(yMargins), ncol=1));
  }

  for (ff in seq_along(yMargins)) {
    xlab <- xlabs[ff];
    main <- mains[ff];
    y <- yMargins[[ff]];

    x <- seq_along(y);
    if (rotate == -90) {
      x <- rev(x);
    }

    fit <- list();
    fit[[1]] <- .robustSmoothSpline(x, y, spar=0.3);
    fit[[2]] <- .robustSmoothSpline(x, y, spar=0.9);

    if (rotate == 0) {
      plot(x, y, pch=pch, cex=cex, ylim=ylim, xlab=xlab, ylab=ylab);
      abline(h=0, lwd=2, col="gray");
      side <- 3;
    } else if (rotate == 90) {
      mar <- par("mar");
      mar <- mar[c(2,3,4,1)];
      par(mar=mar);
      plot(y, x, pch=pch, cex=cex, xlim=ylim, ylab=xlab, xlab=ylab, axes=FALSE);
      axis(side=1); axis(side=4); box();
      abline(v=0, lwd=2, col="gray");
      for (kk in 1:2) {
        t <- fit[[kk]]$x;
        fit[[kk]]$x <- fit[[kk]]$y;
        fit[[kk]]$y <- t;
      }
      side <- 2;
    } else if (rotate == -90) {
      mar <- par("mar");
      mar <- mar[c(2,3,4,1)];
      par(mar=mar);
      plot(y, x, pch=pch, cex=cex, xlim=ylim, ylab=xlab, xlab=ylab, axes=FALSE);
      axis(side=1); axis(side=4); box();
      abline(v=0, lwd=2, col="gray");
      for (kk in 1:2) {
        t <- fit[[kk]]$x;
        fit[[kk]]$x <- fit[[kk]]$y;
        fit[[kk]]$y <- t;
      }
      side <- 2;
    }
    lines(fit[[1]], col="blue", lwd=3);
    lines(fit[[2]], col="red", lwd=5);

    # Plot annotation
    stext(side=side, pos=0, sprintf("Array: %s", getFullName(df)), cex=0.7);
    stext(side=side, pos=1, sprintf("Chip type: %s", getChipType(df)), cex=0.7);
  }

  invisible(yMargins);
})


##############################################################################
# HISTORY:
# 2008-08-19
# o Now writeImages() takes argument 'arrays'.
# 2008-03-17
# o Added readRawDataRectangle(), calculateMargins(), plotMargins(). Will
#   probably be moved elsewhere.
# o Added support for a reference (file).
# 2007-08-09
# o Now getColorMaps(parsed=FALSE) returns a unique sorted set of color maps.
# 2007-03-19
# o Updated addColorMap() to accept multiple transforms.
# o Created from ArrayExplorer.R.
##############################################################################
