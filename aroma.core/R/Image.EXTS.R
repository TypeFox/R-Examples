##############################################################################
# This source file contains all methods and classes related to the
# Image class (of EBImage).
#
# In aroma.affymetrix, the following methods are used:
# - getImage() for the matrix class.
# - display().
# - writeImage().
##############################################################################

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# EBImage Image class
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
setMethodS3("display", "Image", function(this, ...) {
  EBImage::display(this, ...);
}, protected=TRUE)


setMethodS3("writeImage", "Image", function(x, file, ...) {
  EBImage::writeImage(x, files=file, ...);
}, protected=TRUE)


setMethodS3("getImageData", "Image", function(this, ...) {
  x <- this@.Data;
  x;
}, protected=TRUE)


setMethodS3("setImageData", "Image", function(this, data, ...) {
  this@.Data <- data;
  invisible(this);
}, protected=TRUE)



# given an input Image, transform to new RGB image based
# on list of colours given in argument 'palette'
# and lower and upper bounds given in 'lim'.
setMethodS3("colorize", "Image", function(this, palette=gray.colors(256), lim=c(-Inf,Inf), outlierCol="white", ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  getColorMode <- function(what=c("grayscale", "color"), ...) {
    # Argument 'what':
    what <- match.arg(what);

    if (what == "grayscale") {
      colorMode <- EBImage::Grayscale;
    } else if (what == "color") {
      colorMode <- EBImage::Color;
    }

    colorMode;
  } # getColorMode()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'this':
  colorMode <- EBImage::colorMode(this);
  if (colorMode != EBImage::Grayscale) {
    throw("Cannot colorize() non-Grayscale Image: ", getColorMode("color"));
  }

  # Argument 'lim':
  lim <- Arguments$getDoubles(lim, length=c(2,2));
  stopifnot(lim[1] < lim[2]);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  verbose && enter(verbose, "Colorizing grayscale Image");
  verbose && cat(verbose, "Image before:");
  verbose && cat(verbose, "Dimensions: ", paste(dim(this), collapse="x"));
  verbose && cat(verbose, "Color mode: ", colorMode);

  verbose && enter(verbose, "Creating RGB palette");
  verbose && cat(verbose, "Palette:");
  verbose && str(verbose, palette);


  verbose && enter(verbose, "Standardizing pixel intensities to [0,1]");
  # Data
  x <- getImageData(this);
  dim <- dim(x);
  verbose && cat(verbose, "Before:");
  verbose && print(verbose, summary(as.vector(x)));

  # Outliers
  if (lim[1] > -Inf)
    x[x < lim[1]] <- NA;
  if (lim[2] < +Inf)
    x[x > lim[2]] <- NA;
  x[!is.finite(x)] <- NA;

  # Standardize to [0,1]
  r <- range(x, na.rm=TRUE);
  x <- (x - r[1])/(r[2]-r[1]);
  verbose && cat(verbose, "After:");
  verbose && print(verbose, summary(as.vector(x)));

  # Sanity check
  stopifnot(all(0 <= x & x <= 1, na.rm=TRUE));
  verbose && exit(verbose);


  verbose && enter(verbose, "Binning image signals to [1,n]");
  # Standardize to [0,n] where n is number of bins
  n <- length(palette);
  verbose && cat(verbose, "Number of colors (in palette): n=", n);
  colorIdx <- x*n;
  colorIdx <- as.integer(colorIdx);

  # Bin to [1,n] by moving any values at the left extreme to bin one.
  colorIdx[colorIdx == 0L] <- 1L;

  verbose && print(verbose, summary(as.vector(colorIdx)));

  # Sanity check
  stopifnot(all(1 <= colorIdx & colorIdx <= n, na.rm=TRUE));
  verbose && exit(verbose);


  verbose && enter(verbose, "Mapping grayscale intensities to RGB intensities");
  verbose && cat(verbose, "Palette (R,G,B) in [0,1]x[0,1]x[0,1]:");
  paletteRGB <- col2rgb(palette);
  paletteRGB <- t(paletteRGB);
  paletteRGB <- paletteRGB / 256;
  verbose && str(verbose, paletteRGB);
  verbose && exit(verbose);

  z <- paletteRGB[colorIdx,,drop=FALSE];

  # Sanity check
  stopifnot(dim(z)[1] == prod(dim));
  stopifnot(dim(z)[2] == 3L);
  dim(z) <- c(dim, 3L);

  verbose && cat(verbose, "Remapped pixel intensity:");
  verbose && str(verbose, z);

  # Sanity check
  z <- Arguments$getNumerics(z, range=c(0,1));
  verbose && exit(verbose);


  verbose && enter(verbose, "Creating a EBImage::Image object");
  img <- EBImage::Image(data=z, dim=dim(z), colormode=getColorMode("color"));

  verbose && exit(verbose);

  verbose && exit(verbose);

  img;
}, protected=TRUE)




setMethodS3("interleave", "Image", function(this, what=c("none", "h", "v", "auto"), ..., verbose=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  safeMeans <- function(x) {
    mean(x[is.finite(x)]);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'what':
  what <- match.arg(what);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Interleave
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Nothing todo?
  if (what == "none")
    return(this);

  verbose && enter(verbose, "Interleaving image");

  # Get image data
  z <- getImageData(this);
  verbose && cat(verbose, "z:");
  verbose && str(verbose, z);
  zDim <- dim(z);
  ndim <- length(zDim);
  dim(z) <- zDim;

  # Sanity check
  stopifnot(ndim == 2 || ndim == 3);
  if (ndim == 2) {
    dim(z) <- c(zDim, 1L);
  }

  # if only PM locations have signal, add a fake row
  if (what == "auto") {
    verbose && enter(verbose, "Infering horizontal, vertical, or no interleaving");
    n <- 2*(nrow(z) %/% 2);
    idxOdd <- seq(from=1, to=n, by=2);
    zOdd <- z[idxOdd,,,drop=FALSE];
    zEven <- z[idxOdd+1,,,drop=FALSE];
    hOdd <- safeMeans(abs(zOdd));
    hEven <- safeMeans(abs(zEven));
    verbose && printf(verbose, "hOdd=%.2g\n", hOdd);
    verbose && printf(verbose, "hEven=%.2g\n", hEven);
    hRatio <- log(hOdd/hEven);
    if (is.na(hRatio)) hRatio <- Inf;
    verbose && printf(verbose, "hRatio=%.2g\n", hRatio);

    n <- 2*(ncol(z) %/% 2);
#    n <- max(n, 40);  # Infer from the first 40 rows.
    idxOdd <- seq(from=1, to=n, by=2);
    zOdd <- z[,idxOdd,,drop=FALSE];
    zEven <- z[,idxOdd+1,,drop=FALSE];
    vOdd <- safeMeans(abs(zOdd));
    vEven <- safeMeans(abs(zEven));
    verbose && printf(verbose, "vOdd=%.2g\n", vOdd);
    verbose && printf(verbose, "vEven=%.2g\n", vEven);
    vRatio <- log(vOdd/vEven);
    if (is.na(vRatio)) vRatio <- Inf;
    verbose && printf(verbose, "vRatio=%.2g\n", vRatio);

    what <- "none";
    if (abs(vRatio) > abs(hRatio)) {
      if (abs(vRatio) > 0.25) {
        what <- "v"
      }
    } else {
      if (abs(hRatio) > 0.25) {
        what <- "h";
      }
    }
    verbose && cat(verbose, "what: ", what);
    verbose && exit(verbose);
  }

  isUpdated <- FALSE;
  if (what == "h") {
    idxOdd <- seq(from=1, to=2*(nrow(z) %/% 2), by=2);
    z[idxOdd,,] <- z[idxOdd+1,,,drop=FALSE];
  } else if (what == "v") {
    idxOdd <- seq(from=1, to=2*(ncol(z) %/% 2), by=2);
    z[,idxOdd,] <- z[,idxOdd+1,,drop=FALSE];
  } else {
    isUpdated <- FALSE;
  }

  if (ndim == 2) {
    z <- z[,,1,drop=TRUE];
  }

  # Update?
  if (isUpdated) {
    this <- setImageData(this, z);
  }

  verbose && exit(verbose);

  this;
}, protected=TRUE)


setMethodS3("rescale", "Image", function(this, scale=1, blur=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'scale':
  scale <- Arguments$getNumeric(scale, range=c(0,Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Nothing to do?
  if (scale == 1)
    return(this);

  verbose && enter(verbose, "Rescaling image");
  verbose && sprintf(verbose, "Scale: %.2g\n", scale);
  img <- EBImage::resize(this, w=scale*dim(this)[1], blur=blur);
  verbose && exit(verbose);

  img;
}, protected=TRUE)




# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# Images from matrices
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
setMethodS3("getImage", "matrix", function(z, ..., palette=NULL) {
  img <- as.GrayscaleImage(z, ...);

  if (!is.null(palette)) {
    img <- colorize(img, palette=palette, ...);
  }

  img;
}, protected=TRUE)



setMethodS3("createImage", "matrix", function(z, dim=NULL, colorMode=c("gray", "color"), ..., class=NULL, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'z':
  z <- Arguments$getNumerics(z, range=c(0,1));

  # Argument 'class':
  knownClasses <- c("png::array", "EBImage::Image");
  if (is.null(class)) {
    class <- getOption(aromaSettings, "output/ImageClasses", knownClasses);
  }
  class <- match.arg(class, choices=knownClasses, several.ok=TRUE);

  # Argument 'dim':
  if (!is.null(dim)) {
    dim <- Arguments$getIntegers(dim, range=c(0,Inf), length=c(2,2));
  }

  # Argument 'colorMode':
  colorMode <- match.arg(colorMode);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Creating image object");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create an EBImage Image object?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  img <- NULL;

  for (kk in seq_along(class)) {
    verbose && enter(verbose, sprintf("Class #%d ('%s') of %d",
                                        kk, class[kk], length(class)));

    if (class[kk] == "EBImage::Image") {
      tryCatch({
        if (colorMode == "gray") {
          colormode <- EBImage::Grayscale;
        } else if (colorMode == "color") {
          colormode <- EBImage::Color;
        }
        z <- t(z);
        if (is.null(dim)) {
          dim <- dim(z);
        }
        img <- EBImage::Image(data=z, dim=dim, colormode=colormode);
      }, error = function(ex) {
        verbose && print(verbose, ex);
      })
    } else if (class[kk] == "png::array") {
      if (is.null(dim)) {
        dim <- dim(z);
      }

      tryCatch({
        # You can create a 'RasterImage' object without the
        # 'png' package being installed, but we will need it
        # later if the image should be saved to file.
        # If only requested to display on screen or save by other
        # means, the 'png' package is not needed.  This calls
        # for an option to specify "what the purpose is".
        # /HB 2011-02-24
        img <- RasterImage(z);
        if (colorMode == "color") {
          img <- colorize(z);
        }
      }, error = function(ex) {
        verbose && print(verbose, ex);
      })
    }

    # Success?
    if (!is.null(img)) {
      verbose && cat(verbose, "Image was successfully created (using '%s').",
                                                                  class[kk]);
      verbose && exit(verbose);
      break;
    }

    verbose && exit(verbose);
  } # for (kk ...)

  if (is.null(img)) {
    throw("Failed to create image object trying several methods: ", paste(class, colllapse=", "));
  }

  verbose && exit(verbose);

  img;
}, protected=TRUE)



###########################################################################/**
# @set "class=matrix"
# @RdocMethod as.GrayscaleImage
#
# @title "Creates a Grayscale (Color) Image from a matrix file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{z}{A KxN @matrix.}
#   \item{transforms}{A @list of transform @functions.}
#   \item{interleaved}{A @character string specifying how the image data
#     should be interleaved, if at all.}
#   \item{scale}{A @numeric scale factor in (0,+Inf) for resizing the
#     imaging. If \code{1}, no resizing is done.}
#   \item{...}{Passed to \code{colorize()} for the object created.}
#   \item{verbose}{A @logical or a @see "R.utils::Verbose" object.}
# }
#
# \value{
#   Returns a bitmap image object.
# }
#
# @author "HB, KS"
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword internal
#*/###########################################################################
setMethodS3("as.GrayscaleImage", "matrix", function(z, transforms=NULL, interleaved=c("none", "h", "v", "auto"), scale=1, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'transforms':
  if (!is.null(transforms)) {
    if (!is.list(transforms)) {
      transforms <- list(transforms);
    }
  }

  for (transform in transforms) {
    if (!is.function(transform)) {
      throw("Argument 'transforms' contains a non-function: ",
                                                        mode(transform));
    }
  }

  # Argument 'interleaved':
  interleaved <- match.arg(interleaved);

  # Argument 'scale':
  scale <- Arguments$getNumeric(scale, range=c(0,Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Creating Image object from matrix");

  # Transform signals?
  for (transform in transforms) {
    dim <- dim(z);
    z <- transform(z);
    dim(z) <- dim;
  }

  # Transform into [0,1]
  # EBImage: "Grayscale values are assumed to be in the range [0,1],
  # although this is not a requirement in sense of data storage. Although,
  # many image processing functions will assume data in this range or will
  # generate invalid results for the data out of this range."
  verbose && enter(verbose, "Rescaling to [0,1]");
  r <- range(z, na.rm=TRUE, finite=TRUE);
  z <- (z - r[1])/(r[2]-r[1]);
  verbose && summary(verbose, as.vector(z));
  verbose && exit(verbose);

  verbose && enter(verbose, "Censor +/-Inf to NA");
  z[is.infinite(z)] <- NA;
  verbose && summary(verbose, as.vector(z));
  verbose && exit(verbose);

  # Create an Image object
  # Sanity check
  stopifnot(all(z >= 0, na.rm=TRUE));
  stopifnot(all(z <= 1, na.rm=TRUE));
  img <- createImage(z, colorMode="gray", ...);

  verbose && cat(verbose, "Create image object:");
  verbose && print(verbose, img);

  # if only PM locations have signal, add a fake row?
  img <- interleave(img, what=interleaved);

  # Scale (if scale == 1, it does nothing)
  img <- rescale(img, scale=scale);

  verbose && exit(verbose);

  img;
}, protected=TRUE)



############################################################################
# HISTORY:
# 2013-09-11
# o CLEANUP: Backward-compatibility code for EBImage (< 2.1.23) was
#   no longer needed.
# 2011-03-06
# o BUG FIX: Update colorize() for Image to work with "newer"
#   EBImage (>= 3.9.7), it may break older version fo EBImage.
# o ROBUSTNESS: Now asserting 'z' is in [0,1] in createImage(z, ...).
# o CLEANUP: Dropped deprecated as.TrueColorImage() for Image and matrix.
# 2011-03-05
# o Now createImage() for matrix creates images of type "png::array"
#   by default (used to be "EBImage::Image").
# 2011-02-27
# o BUG FIX: colorize() for Image would throw '<simpleError in ...: could
#   not find function "colorMode">', because the colorMode() function
#   needs to be explicitly imported after the recent package cleanups.
# o BUG FIX: colorize() for Image try to call createImage() using a vector
#   and not a matrix.
# 2011-02-24
# o GENERALIZATION: Now the default for createImage() for matrix is to
#   test to create images according to aroma settings option
#   'output/ImageClasses'.
# o BUG FIX: createImage() for matrix would not result the first possible
#   image created (when testing different image classes) but instead
#   continue trying to create image for all possible classes.
#   For instance, this meant that although you had the 'EBImage' package
#   installed, but not the 'png' package, it would still in the end try
#   to (also) use 'png' package.  If writing PNG images to file, say via
#   ArrayExplorer, this would result in "Error in loadNamespace(name) :
#   there is no package called 'png'".  Thanks Richard Beyer at
#   University of Washington for reporting on this.
# 2011-01-31
# o Added createImage() with argument 'class' for specifying what type of
#   bitmap image should be created.  Currently only Image objects of
#   EBImage can be created.
# o CLEAN UP: Removed require("EBImage") from resize().
# o CLEAN UP: Removed non-used local safeMean() from as.GrayscaleImage().
# o CLEAN UP: Removed explicit dependencies on EBImage from interleave().
# o Added getImageData() and setImageData() for the Image class.
# o Made rescale() and interleave() protected.
# o CLEAN UP: Removed deprecated internal rgbTransform() for Image.
# o Moved all methods related to the Image class to one source file.
#   All history has been merged accordingly.
# 2011-01-30
# o Added writeImage() in order to clean up aroma.affymetrix, cf.
#   writeImage() for AffymetrixCelFile.
# 2010-06-22
# o BUG FIX: as.GrayscaleImage(..., transforms=NULL) for 'matrix' would
#   throw "Exception: Argument 'transforms' contains a non-function: NULL".
# 2009-05-16
# o Now argument 'scale' of as.GrayscaleImage() is validated using
#   Arguments$getNumeric() [not getDouble()].
# o Now rescale() for Image uses Arguments$getNumerics(), not
#   getDoubles(), where possible.  This will save memory in some cases.
# 2009-05-10
# o Forgot argument 'verbose' for interleave().
# o BUG FIX: interleave() for Image gave 'Error in z[idxOdd,, ] : incorrect
#   number of dimensions'.  The internal image structure is now a 2-dim
#   matrix, again(?!?).
# 2008-10-16
# o BUG FIX: Tried to turn a function passed by 'transforms' to a list using
#   as.list() and not list().
# 2008-04-13
# o Added asserting to colorize() that the Image has color mode Grayscale.
# o Renamed rgbTransform() to colorize().
# 2008-04-14
# o Added display().
# 2008-05-10
# o BUG FIX: interleave() for Image gave 'Error in z[idxOdd, ] : incorrect
#   number of dimensions'.  The internal image structure is a 3-dim array.
# 2008-03-14
# o Added rescale() from getImage() of AffymetrixCelFile.
# o Added interleave() from from getImage() of AffymetrixCelFile.
# o Added getImage().
# o Added as.GrayscaleImage() and as.TrueColorImage().
# o Created from getImage() of AffymetrixCelFile.
# 2007-06-11
# o Explicit calls to EBImage::Image() etc in colorize().
# 2007-01-30 /HB
# o Improved the speed of colorize() even further by dirt simple binning.
# o Move to Image.rgbTransform.R.
# o Improved speed of colorize() on the internal binning in rgbTransform().
# 2007-01-12 /KS
# o TOFIX: THIS [colorize()] IS CURRENTLY TOO SLOW.
# o Created. colorize() is used by getImage() in AffymetrixCelFile.
############################################################################
