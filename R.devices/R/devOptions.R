###########################################################################/**
# @RdocFunction devOptions
# @alias getDevOption
#
# @title "Gets the default device options"
#
# \description{
#  @get "title" as given by predefined devices options adjusted for
#  the default arguments of the device function.
# }
#
# @synopsis
#
# \arguments{
#   \item{type}{A @character string or a device @function specifying
#      the device to be queried.  May also be a @vector, for querying
#      device options for multiple devices.}
#   \item{custom}{If @TRUE, also the default settings specific to this
#      function is returned. For more details, see below.}
#   \item{special}{A @logical.  For more details, see below.}
#   \item{inherits}{If @TRUE, the global option is used if the
#      type-specific is not set (or @NULL).}
#   \item{drop}{If @TRUE and only one device type is queried, then
#      a @list is returned, otherwise a @matrix.}
#   \item{options}{Optional named @list of settings.}
#   \item{...}{Optional named arguments for setting new defaults.
#      For more details, see below.}
#   \item{reset}{If @TRUE, the device options are reset to R defaults.}
# }
#
# \value{
#   If \code{drop=TRUE} and a single device is queries, a named @list is
#   returned, otherwise a @matrix is returned.
#   If a requested device is not implemented/supported on the current system,
#   then "empty" results are returned.
#   If options were set, that is, if named options were specified via
#   \code{...}, then the list is returned invisibly, otherwise not.
# }
#
# \details{
#  If argument \code{special} is @TRUE, then the 'width' and 'height'
#  options are adjusted according to the rules explained for
#  argument 'paper' in @see "grDevices::pdf", @see "grDevices::postscript",
#  and @see "grDevices::xfig".
# }
#
# \section{Setting new defaults}{
#  When setting device options, the \code{getOption("devOptions")[[type]]}
#  option is modified.  This means that for such options to be effective,
#  any device function needs to query also such options, which for instance
#  is done by @see "devNew".
#
#  Also, for certain devices (eps, pdf, postscript, quartz, windows and x11),
#  builtin R device options are set.
# }
#
# @examples "../incl/devOptions.Rex"
#
# @author
#
# @keyword device
# @keyword utilities
#*/###########################################################################
devOptions <- function(type=NULL, custom=TRUE, special=TRUE, inherits=FALSE, drop=TRUE, options=list(), ..., reset=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local setups
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  devList <- list(
    ## Global device options
    "*"=NA_character_
  )
  ## All known and supported graphics devices on this system
  devList <- c(devList, devAll())
  knownTypes <- names(devList)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (.Platform$OS.type == "windows") {
    # To please R CMD check
    windows.options <- NULL; rm(list="windows.options");
    x11.options <- grDevices::windows.options;
    X11.options <- x11.options;
  }


  # Global options; only used for its ability to reset them
  "*.options" <- function(..., reset=FALSE) {
    if (reset) {
      devOptions(type="*", sep=",", path="figures", field=NULL, force=TRUE)
    }
  }


  # A template for a dummy device options function.
  getNnnOptions <- function(type, fallbacks=FALSE, ...) {
    optList <- list(
      "*"="*.options",
      eps="ps.options",
      jpeg2="ps.options",
      pdf="pdf.options",
      png2="ps.options",
      postscript="ps.options",
      quartz="quartz.options",
      windows="windows.options",
      x11="X11.options",
      X11="x11.options"
    );

    if (fallbacks) {
      optList <- c(optList, list(
        "CairoPDF"="pdf.options",
        "CairoPS"="ps.options",
        "CairoX11"="x11.options"
      ))
    }

    dummy <- function(...) { list(); }

    if (!is.element(type, names(optList))) {
      return(dummy);
    }

    key <- optList[[type]];

    # Sanity check
    stopifnot(length(key) == 1L);

    # Does the nnn.function() already exists?
    envir <- getNamespace("grDevices");
    if (exists(key, envir=envir, mode="function")) {
      fcn <- get(key, envir=envir, mode="function");
      return(fcn);
    } else if (exists(key, mode="function")) {
      fcn <- get(key, mode="function");
      return(fcn);
    }

    # If not, create either a real one or a dummy one...
    typeC <- capitalize(type);
    keyE <- sprintf(".%senv", typeC);
    if (exists(keyE, envir=envir, mode="environment")) {
      # A real one
      envir <- get(keyE, envir=envir, mode="environment");
      fcn <- function(..., reset=FALSE) {
        keyO <- sprintf(".%s.Options", typeC);
        opts <- get(keyO, envir=envir, mode="list");
        if (reset) {
          keyD <- sprintf("%s.default", keyO);
          opts <- get(keyD, envir=envir, mode="list");
          assign(keyD, value=opts, envir=envir);
        }
        opts;
      };
    } else {
      # A dummy
      fcn <- dummy;
    }

    fcn;
  } # getNnnOptions()


  # See argument 'paper' in help("pdf"), help("postscript"), and
  # help("xfig").
  getSpecialDimensions <- function(options=list(), sizes=names(paperSizes), ...) {
    paperSizes <- list(
      a4        = c( 8.27, 11.69),
      a4r       = c(11.69,  8.27),
      executive = c( 7.25, 10.5 ),
      legal     = c( 8.5 , 14   ),
      letter    = c( 8.5 , 11   ),
      USr       = c(11   , 8.5  )
    );

    paper <- tolower(options$paper);
    if (paper == "default") {
      paper <- getOption("papersize", "a4");
    }

    if (paper != "special") {
      paperSizes <- paperSizes[sizes];
      dim <- paperSizes[[paper]];

      # Replace "special" 0:s with NA:s, to indicate they are missing
      dim[dim == 0L] <- NA_real_;
    } else {
      dim <- c(options$width, options$height);
    }

    dim;
  } # getSpecialDimensions()


  getDevOptions <- function(type, ...) {
    opts <- getOption("devOptions", list());
    opts <- opts[[type]];
    if (is.null(opts)) {
      opts <- list();
    }
    opts;
  } # getDevOptions()


  setDevOptions <- function(.type, ..., reset=FALSE) {
    devOpts <- getOption("devOptions", list());
##    str(list(.type=.type, devOpts=devOpts));
    oopts <- opts <- devOpts[[.type]];
    if (is.null(opts)) opts <- list();
    # Sanity check
    stopifnot(is.list(opts));
##    str(list(opts=opts));
    if (reset) {
      opts <- NULL;
    } else {
      args <- list(...);
      if (length(args) > 0L) {
        for (key in names(args)) {
          value <- args[[key]];
##          str(list(key=key, value=value))
          opts[[key]] <- value;
        }
      }
    }
    devOpts[[.type]] <- opts;
    options(devOptions=devOpts);
    invisible(oopts);
  } # setDevOptions()


  getDeviceFunctions <- function(type, default=function() {}, ...) {
    # Find all device functions.  If non-existent, return the default
    devs <- devList[[type]];

    parts <- strsplit(devs, split="::", fixed=TRUE);
    devs <- lapply(parts, FUN=function(s) {
      if (length(s) > 1L) {
        # Device package may not be installed, e.g. Cairo and JavaGD.
        envir <- tryCatch({
          # Attach it as well, so that the device function is found
          # when needed inside devNew() at do.call().
          suppressWarnings({
            if (require(s[1L], character.only=TRUE, quietly=TRUE)) {
              # Return namespace
              getNamespace(s[1L]);
            } else {
              emptyenv();
            }
          })
        }, error = function(ex) emptyenv());
        s <- s[-1L];
      } else {
        envir <- as.environment(-1L);
      }

      if (exists(s[1L], envir=envir, mode="function")) {
        get(s[1L], envir=envir, mode="function");
      } else {
        default;
      }
    });

    devs;
  } # getDeviceFunctions()


  findDeviceFunction <- function(fcn, ...) {
    types <- names(devList)

    ## (i) Scan primary device functions
    for (type in types) {
      devs <- getDeviceFunctions(type, default=NULL)
      if (identical(fcn, devs[[1]])) return(type)
    }

    ## (ii) Scan non-primary device functions
    for (type in types) {
      devs <- getDeviceFunctions(type, default=NULL)
      for (dev in devs[-1]) {
        if (identical(fcn, dev)) return(type)
      }
    }

    "<function>"
  } # findDeviceFunction()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'type':
  if (is.null(type)) {
    type <- names(devList)
  }

  # Argument 'options':
  if (!is.list(options)) {
    throw("Argument 'options' must be a list: ", class(options)[1L]);
  }
  nopts <- length(options);
  if (nopts > 0L && is.null(names(options))) {
    throw("Argument 'options' must be a named list.");
  }

  # Argument 'inherits':
  inherits <- as.logical(inherits);

  # Argument 'reset':
  reset <- as.logical(reset);

  # Additional arguments
  args <- list(...);
  nargs <- length(args);
  if (nargs > 0L && is.null(names(args))) {
    throw("Optional ('...') arguments must be named.");
  }

  # Append/overwrite 'options' with the named argument, e.g. field=NULL
  for (key in names(args)) {
    options[key] <- list(args[[key]]);
  }
  nopts <- length(options);
  # Not needed anymore
  args <- nargs <- NULL;

  # Argument 'type':
  if (missing(type) || length(type) == 0L) {
    if (nopts > 0L) {
      throw("Cannot set device options. Argument 'type' is missing or NULL. Should be one of: ", paste(sQuote(knownTypes), collapse=", "));
    }

    res <- devOptions(type=knownTypes, custom=custom, special=special, inherits=inherits, drop=drop, reset=reset);
    if (reset) {
      return(invisible(res));
    } else {
      return(res);
    }
  } else if (is.character(type)) {
    # Expand by regexp matching, iff any
    type <- .devTypeName(type, pattern=TRUE, knownTypes=knownTypes);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Vectorized call?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (length(type) > 1L) {
    if (nopts > 0L) {
      throw("Cannot set device options for more than one devices at the time: ", hpaste(sQuote(type), collapse=", "));
    }

    # Support vector of 'type':s
    types <- type;
    res <- lapply(types, FUN=devOptions, inherits=inherits, drop=TRUE);
    fields <- lapply(res, FUN=function(opts) names(opts));
    fields <- unique(unlist(fields, use.names=FALSE));
    opts <- lapply(res, FUN=function(x) x[fields]);
    opts <- unlist(opts, use.names=FALSE, recursive=FALSE);
    dim(opts) <- c(length(fields), length(types));
    dimnames(opts) <- list(fields, types);
    opts <- t(opts);
    ##class(opts) <- c("DeviceOptions", class(opts));
    return(opts);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # A single device at this point
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.function(type)) {
    # Try to find name of device function
    type <- findDeviceFunction(fcn=type);
  }
  if (is.character(type)) {
    if (!is.element(type, knownTypes)) {
      throw(sprintf("Device type %s is not known/supported on this operating system/platform. Supported devices types are: %s", sQuote(type), paste(sQuote(setdiff(knownTypes, "*")), collapse=", ")))
    }
  }

  if (!is.element(type, names(devList))) {
    throw("Cannot query/modify device options. Unknown device: ", type);
  }

  if (inherits) {
    throw("Argument inherits=TRUE is not yet supported.");
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locate the nnn.options() function for this type of device
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nnn.options <- getNnnOptions(type);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locate the list of device functions used by this type of device
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # If non-existent, return a dummy
  devs <- getDeviceFunctions(type, default=function() {});


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reset?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (reset) {
    res <- devOptions(type=type, special=special, drop=drop, reset=FALSE);
    setDevOptions(type, reset=TRUE);

    # Only for certain devices...
    nnn.options(reset=TRUE);

    return(invisible(res));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Assign user arguments, iff possible
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (nopts > 0L) {
    do.call(setDevOptions, args=c(list(.type=type), options));

    # Only for certain devices...
    # (a) Try to set all options at once
    ok <- tryCatch({
      do.call(nnn.options, args=options);
      TRUE
    }, error = function(ex) FALSE)

    # (b) Otherwise, try to set them one by one
    if (!ok) {
      for (kk in seq_along(options)) {
        res <- tryCatch({
          do.call(nnn.options, args=options[kk]);
        }, error = function(ex) {})
      }
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get builtin device options, if available
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Only for certain devices...
  nnn.options <- getNnnOptions(type, fallbacks=TRUE);
  opts <- nnn.options();


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get (nested) device formals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  defArgs <- lapply(rev(devs), FUN=formals);
  defArgs <- Reduce(append, defArgs);
  # Drop '...'
  defArgs <- defArgs[names(defArgs) != "..."];
  # Drop missing
  # Drop overridden values
  defArgs <- defArgs[!duplicated(names(defArgs), fromLast=TRUE)];


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Merge arguments that either are not in the predefined set of
  # device options, or ones that replaced the default value.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  optsT <- c(opts, defArgs);


  # Include custom options specific to devOptions()?
  if (custom) {
    devOpts <- getDevOptions(type);
    optsT <- c(optsT, devOpts);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Build the concensus of all options
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (a) Keep all non-duplicated options
  dups <- duplicated(names(optsT));
  idxsA <- which(!dups);
  optsA <- optsT[idxsA];

  # (b) Among duplicates, keep those with values
  idxsB <- which(dups);
  idxsB <- idxsB[!sapply(optsT[idxsB], FUN=is.symbol)];
  optsB <- optsT[idxsB];

  opts <- append(optsA, optsB);

  # Drop overridden values
  opts <- opts[!duplicated(names(opts), fromLast=TRUE)];


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Adjust for special cases?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (special) {
    if (is.element(type, c("eps", "postscript", "CairoPS"))) {
      sizes <- c("a4", "executive", "legal", "letter");
      dim <- getSpecialDimensions(opts, sizes);
    } else if (type == "xfig") {
      sizes <- c("a4", "legal", "letter");
      dim <- getSpecialDimensions(opts, sizes);
    } else if (is.element(type, c("pdf", "CairoPDF"))) {
      sizes <- c("a4", "a4r", "executive", "legal", "letter", "USr");
      dim <- getSpecialDimensions(opts, sizes);
    } else {
      dim <- NULL;
    }

    if (!is.null(dim)) {
      opts$width <- dim[1L];
      opts$height <- dim[2L];
    }
  }

  # Return a list?
  if (!drop) {
    names <- names(opts);
    dim(opts) <- c(1L, length(opts));
    dimnames(opts) <- list(type, names);
    ## Possibly, class(opts) <- c("DeviceOptions", class(opts))
  }

  # Return invisibly?
  if (nopts > 0L) {
    invisible(opts);
  } else {
    opts;
  }
} # devOptions()


getDevOption <- function(type, name, default=NULL, inherits=TRUE, ..., old=TRUE) {
  # BACKWARD COMPATIBILITY:
  if (isTRUE(old)) old <- sprintf("devEval/args/%s", name)
  if (is.character(old)) {
    .importOldGlobalOption(name, from=old, remove=TRUE)
  }

  # Get options
  if (is.character(type) && type != "*") {
    value <- devOptions(type=type, ...)[[name]]
  } else {
    value <- NULL
  }

  # Fallback to global options?
  if (is.null(value) && inherits) {
    value <- devOptions(type="*", ...)[[name]]
  }

  # Use default value?
  if (is.null(value)) value <- default

  value
} # getDevOption()


# BACKWARD COMPATIBILITY for old-style R options().
.importOldGlobalOption <- function(name, from, remove=TRUE, ...) {
  # Nothing to do?
  if (!is.element(from, names(options()))) return()

  # (a) Get old-style option value
  opt <- options(from);

  # (b) Assign to global device options instead
  names(opt) <- name;
  do.call(devOptions, args=list(type="*", options=opt));

  # (c) Delete old-style option and never look back
  if (remove) setOption(from, NULL);
} # .importOldGlobalOption()


############################################################################
# HISTORY:
# 2014-09-16
# o Added support for regular expression matching for argument 'type'.
# o BUG FIX: devOptions() returned the incorrect options for device
#   types "eps", "jpg2" and "png2" if package was not attached.
# 2014-09-15
# o Added "favicon" as a known type.
# 2014-09-12
# o Added support for global options via devOptions("*").
# o BUG FIX: On Windows, devOptions() assumed that the 'grDevices'
#   package was attached.
# o BUG FIX: devOptions(type, name=NULL) did not assign the option NULL.
# 2014-09-12
# o Added getDevOption().
# 2013-12-08
# o BUG FIX: devOptions(types) would drop all options for combinations
#   devices types that have identical sets of options, e.g.
#   types=c("png", "png") or types=c("bmp", "png").
# 2012-07-26
# o Added arguments 'options' and 'drop' to devOptions().
# 2012-07-24
# o Now devOptions() supports a vector of devices types.  If so, then
#   a matrix where each row specifies a device type.  The union of
#   all possible fields are returned as columns.  Non-existing fields
#   are returned as NULL.  If no device type is specified (or NULL), then
#   all are returned.
# o Now devOptions() also supports "win.metafile".
# 2012-04-30
# o Now devOptions() returns options invisibly if some options were set,
#   otherwise not, e.g. devOptions() versus devOptions("png", width=1024).
# 2012-02-26
# o GENERALIZATION: Now devOptions() accepts passing a device function
#   in addition a sting, e.g. devOptions(png) and devOptions("png").
# 2011-11-07
# o Added 'quarts' to the list of (possible) devices.
# o BUG FIX: devOptions() assumed that all devices exist on
#   all platforms, causing it to give an error on some.
# 2011-11-05
# o Created.
############################################################################
