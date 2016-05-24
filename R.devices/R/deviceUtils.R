###########################################################################/**
# @RdocFunction devIsOpen
#
# @title "Checks if zero or more devices are open or not"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{which}{An index (@numeric) @vector or a label (@character) @vector.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a named @logical @vector with @TRUE if a device is open,
#   otherwise @FALSE.
# }
#
# @examples "../incl/deviceUtils.Rex"
#
# @author
#
# @keyword device
# @keyword utilities
#*/###########################################################################
devIsOpen <- function(which=dev.cur(), ...) {
  # Nothing to do?
  if (length(which) == 0L) {
    res <- logical(0L);
    names(res) <- character(0L);
    return(res);
  }

  devList <- .devList();
  devs <- devList[which];

  labels <- names(devs);
  isOpen <- sapply(devs, FUN=function(dev) !is.null(dev) && nzchar(dev));
  isOpen <- isOpen & !is.na(labels);
  isOpen;
} # devIsOpen()





###########################################################################/**
# @RdocFunction devList
#
# @title "Lists the indices of the open devices named by their labels"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{interactiveOnly}{If @TRUE, only open interactive/screen devices
#     are returned.}
#   \item{dropNull}{If @TRUE, the "null" device (device index 1) is
#     not returned.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a named @integer @vector.
# }
#
# @author
#
# \seealso{
#   \code{\link[grDevices:dev]{dev.list}()}.
# }
#
# @keyword device
# @keyword utilities
#*/###########################################################################
devList <- function(interactiveOnly=FALSE, dropNull=TRUE, ...) {
  devList <- .devList();

  # Return only opened devices
  isOpen <- sapply(devList, FUN=function(dev) (dev != ""));
  names(isOpen) <- names(devList);
  idxs <- which(isOpen);

  # Exclude the "null" device?
  if (dropNull) {
    idxs <- idxs[-1L];
  }

  # Include only interactive devices?
  if (interactiveOnly) {
    types <- unlist(devList[idxs]);
    keep <- devIsInteractive(types);
    idxs <- idxs[keep];
  }

  idxs;
} # devList()



###########################################################################/**
# @RdocFunction devGetLabel
#
# @title "Gets the labels of zero or more devices"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{which}{An index (@numeric) @vector or a label (@character) @vector.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character @vector.
#   If a device does not exist, an error is thrown.
# }
#
# @author
#
# \seealso{
#   @see "devSetLabel" and @see "devList".
# }
#
# @keyword device
# @keyword utilities
#*/###########################################################################
devGetLabel <- function(which=dev.cur(), ...) {
  devList <- devList(dropNull=FALSE);
  if (is.numeric(which)) {
    devs <- devList[match(which, devList)];
  } else {
    devs <- devList[which];
  }
  labels <- names(devs);
  unknown <- which[is.na(labels)];
  if (length(unknown) > 0L) {
    known <- names(devList(dropNull=FALSE));
    if (length(known) == 0L) known <- "<none>";
    throw(sprintf("Cannot get device label. No such device: %s (known devices: %s)", paste(sQuote(unknown), collapse=", "), paste(sQuote(known), collapse=", ")));
  }
  labels;
} # devGetLabel()



###########################################################################/**
# @RdocFunction devSetLabel
#
# @title "Sets the label of a device"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{which}{An index (@numeric) or a label (@character).}
#   \item{label}{A @character string.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns nothing.
# }
#
# @author
#
# \seealso{
#   @see "devGetLabel" and @see "devList".
# }
#
# @keyword device
# @keyword utilities
#*/###########################################################################
devSetLabel <- function(which=dev.cur(), label, ...) {
  # Argument 'which':
  if (length(which) != 1L) {
    throw("Argument 'which' must be a scalar: ", paste(which, collapse=", "));
  }

  devList <- .devList();
  if (is.numeric(which)) {
    idx <- which;
  } else {
    idx <- .devListIndexOf(which);
  }

  # Unknown devices?
  if (devList[[idx]] == "") {
    known <- names(devList(dropNull=FALSE));
    if (length(known) == 0L) known <- "<none>";
    throw(sprintf("Cannot set device label. No such device: %s (known devices: %s)", paste(sQuote(which), collapse=", "), paste(sQuote(known), collapse=", ")));
  }

  # Update the label
  if (is.null(label))
    label <- "";
  names(devList)[idx] <- label

  assign(".Devices", devList, envir=baseenv());
}





###########################################################################/**
# @RdocFunction devSet
#
# @title "Activates a device"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{which}{An index (@numeric) or a label (@character).
#     If neither, then a label corresponding to the checksum of
#     \code{which} as generated by @see "digest::digest" is used.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns what \code{\link[grDevices:dev]{dev.set}()} returns.
# }
#
# @author
#
# \seealso{
#   @see "devOff" and @see "devDone".
#   Internally, \code{\link[grDevices:dev]{dev.set}()} is used.
# }
#
# @keyword device
# @keyword utilities
#*/###########################################################################
devSet <- function(which=dev.next(), ...) {
  args <- list(...);

  # Argument 'which':
  if (!is.numeric(which) || length(which) != 1L) {
    if (length(which) != 1L || !is.character(which)) {
      # To please R CMD check
      requireNamespace("digest") || throw("Package not loaded: digest");
      which <- digest::digest(which);
    }

    if (is.character(which)) {
      args$label <- which
      idx <- .devListIndexOf(which, error=FALSE)
      # If not existing, open the next available one
      if (is.na(idx)) {
        which <- .devNextAvailable()
      } else {
        devList <- devList(dropNull=FALSE)
        which <- devList[[idx]]
      }
    }
  }

  # Argument 'which':
  if (length(which) != 1L) {
    throw("Argument 'which' must be a scalar: ", paste(which, collapse=", "));
  }

  if (which < 2L || which > 63L) {
    throw("Cannot set device: Argument 'which' is out of range [2,63]: ", which);
  }

  if (devIsOpen(which)) {
    # Active already existing device
    return(dev.set(which));
  }

  # Identify set devices that needs to be opened inorder for
  # the next device to get the requested index
  if (which == 2L) {
    toBeOpened <- c();
  } else {
    toBeOpened <- setdiff(2:(which-1L), dev.list());
  }

  len <- length(toBeOpened);
  if (len > 0L) {
    tempfiles <- sapply(toBeOpened, FUN=function(...) tempfile());

    # Make sure to close all temporary devices when exiting function
    on.exit({
      for (kk in seq_along(toBeOpened)) {
        dev.set(toBeOpened[kk]);
        dev.off();
        if (file.exists(tempfiles[kk])) file.remove(tempfiles[kk]);
      }
    }, add=TRUE);

    # Create a dummy temporary postscript device (which is non-visible)
    for (kk in seq_along(toBeOpened)) {
      postscript(file=tempfiles[kk]);
    }
  }

  # Open the device
  res <- do.call(devNew, args=args);

  invisible(res);
} # devSet()




###########################################################################/**
# @RdocFunction devOff
#
# @title "Closes zero or more devices"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{which}{An index (@numeric) @vector or a label (@character) @vector.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns \code{\link[grDevices:dev]{dev.cur}()}.
# }
#
# @author
#
# \seealso{
#   @see "devDone".
#   Internally, \code{\link[grDevices:dev]{dev.off}()} is used.
# }
#
# @keyword device
# @keyword utilities
#*/###########################################################################
devOff <- function(which=dev.cur(), ...) {
  # Nothing to do?
  if (length(which) == 0L) return(dev.cur());

  # Only close each device once
  which <- unique(which);
  if (length(which) > 1L) {
    lapply(which, FUN=devOff);
    return(dev.cur());
  }

  # Nothing to do?
  if (!devIsOpen(which)) {
    return(dev.cur());
  }

  # Identify device
  which <- devSet(which);

  # Reset the label
  devSetLabel(which, label=NULL);

  # Close device
  dev.off(which);

  return(dev.cur());
} # devOff()




###########################################################################/**
# @RdocFunction devDone
#
# @title "Closes zero or more open devices except screen (interactive) devices"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{which}{An index (@numeric) @vector or a label (@character) @vector.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns (invisibly) \code{\link[grDevices:dev]{dev.cur}()}.
# }
#
# @author
#
# \seealso{
#   @see "devOff".
#   @see "grDevices::dev.interactive".
# }
#
# @keyword device
# @keyword utilities
#*/###########################################################################
devDone <- function(which=dev.cur(), ...) {
  # Nothing to do?
  if (length(which) == 0L) return(dev.cur());

  # Only close each device once
  which <- unique(which);
  if (length(which) > 1L) {
    lapply(which, FUN=devDone);
    return(dev.cur());
  }

  # Nothing to do?
  if (!devIsOpen(which)) {
    return(invisible(dev.cur()));
  }

  # Do nothing?
  if (is.numeric(which) && length(which) == 1L && which <= 1L) {
    return(invisible(dev.cur()));
  }

  which <- devSet(which);
  if (which != 1L) {
    type <- tolower(names(which));
    type <- gsub(":.*", "", type);

    knownInteractive <- deviceIsInteractive();
    knownInteractive <- tolower(knownInteractive);
    isOnScreen <- (is.element(type, knownInteractive));
    if (!isOnScreen)
      devOff(which);
  }

  return(invisible(dev.cur()));
} # devDone()


###########################################################################/**
# @RdocFunction devIsInteractive
#
# @title "Checks whether a device type is interactive or not"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{types}{A @character @vector.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @logical @vector with @TRUE if the device type is interactive,
#   otherwise @FALSE.
# }
#
# @author
#
# \seealso{
#   Internally, @see "grDevices::deviceIsInteractive" is used.
# }
#
# @keyword device
# @keyword utilities
#*/###########################################################################
devIsInteractive <- function(types, ...) {
  # Known interactive devices
  knownInteractive <- grDevices::deviceIsInteractive();
  knownInteractive <- c(knownInteractive, "CairoWin", "CairoX11");
##  knownInteractive <- c(knownInteractive, "Cairo");

  # Return all known?
  if (missing(types)) return(knownInteractive);

  # Nothing to do?
  if (length(types) == 0L) {
    res <- logical(0L);
    names(res) <- character(0L);
    return(res);
  }

  if (length(types) > 1L) {
    res <- sapply(types, FUN=devIsInteractive);
    if (is.character(types)) names(res) <- types;
    return(res);
  }


  # Sanity check
  stopifnot(length(types) == 1L);

  # Investigate one type below
  type0 <- type <- types;

  if (is.function(type)) {
    for (name in knownInteractive) {
      if (exists(name, mode="function")) {
        dev <- get(name, mode="function");
        if (identical(dev, type)) {
          res <- TRUE;
          names(res) <- name;
          return(res);
        }
      }
    }
    res <- FALSE;
    names(res) <- NA_character_;
  } else {
    type <- as.character(type);
    # Device type aliases?
    type <- .devTypeName(type);
    type <- tolower(type);

    # An known one?
    res <- is.element(type, tolower(knownInteractive));
    names(res) <- type0;
  }

  res;
} # devIsInteractive()



devAll <- local({
  # Memoize results, because this function is called many many
  # times either directly or indirectly by various functions.
  .devAll <- NULL

  isFALSE <- function(x) identical(FALSE, unname(x))

  base_capabilities <- local({
    res <- base::capabilities()
    function(force=FALSE) {
      if (force) res <<- base::capabilities()
      res
    }
  })

  supports <- function(type, pkg="grDevices") {
    capabilities <- function() character(0L)
    if (isPackageInstalled(pkg)) {
      if (pkg == "Cairo") {
        ns <- getNamespace("Cairo")
        Cairo.capabilities <- get("Cairo.capabilities", envir=ns, mode="function")
        capabilities <- function() {
          res <- Cairo.capabilities()
          names <- sprintf("Cairo%s", toupper(names(res)))
          names <- gsub("WIN", "Win", names, fixed=TRUE)
          names(res) <- names
          res
        }
      } else if (pkg == "grDevices") {
        capabilities <- base_capabilities
      }
    }
    x <- capabilities()
    !isFALSE(x[type])
  } # supports()


  function(force=FALSE, ...) {
    res <- .devAll
    if (force || is.null(res)) {
      if (force) base_capabilities(force=TRUE)

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Setup all possible devices
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      res <- list()

      ## grDevices
      res <- c(res, list(
        bmp          = "grDevices::bmp",
        jpeg         = "grDevices::jpeg",
        pdf          = "grDevices::pdf",
        pictex       = "grDevices::pictex",
        png          = "grDevices::png",
        postscript   = "grDevices::postscript",
        quartz       = "grDevices::quartz",
        svg          = c("grDevices::svg", "needs::cairo"),
        cairo_pdf    = c("grDevices::cairo_pdf", "needs::cairo"),
        cairo_ps     = c("grDevices::cairo_ps", "needs::cairo"),
        tiff         = "grDevices::tiff",
        win.metafile = "grDevices::win.metafile",
        xfig         = "grDevices::xfig"
      ))
      res <- c(res, list(
        windows    = "grDevices::windows",
        x11        = "grDevices::x11",
        X11        = "grDevices::X11"
      ))

      # R.devices
      res <- c(res, list(
        eps     = c("R.devices::eps",
                    "grDevices::postscript"),
        favicon = c("R.devices::favicon",
                    "grDevices::png"),
        jpeg2   = c("R.devices::jpeg2",
                    "grDevices::bitmap", "grDevices::postscript"),
        png2    = c("R.devices::png2",
                    "grDevices::bitmap", "grDevices::postscript")
      ))

      ## Cairo package
      res <- c(res, list(
        CairoPDF  = c("Cairo::CairoPDF",
                      "grDevices::pdf"),
        CairoPS   = c("Cairo::CairoPS",
                      "grDevices::postscript"),
        CairoPNG  = c("Cairo::CairoPNG",
                      "grDevices::png"),
        CairoJPEG = c("Cairo::CairoJPEG",
                      "grDevices::jpeg"),
        CairoTIFF = c("Cairo::CairoTIFF",
                      "grDevices::tiff"),
        CairoSVG  = c("Cairo::CairoSVG",
                      "grDevices::svg")
      ))
      res <- c(res, list(
        CairoWin  = c("Cairo::CairoWin",
                      "grDevices::windows"),
        CairoX11  = c("Cairo::CairoX11",
                      "grDevices::x11")
      ))

      ## JavaGD
      ## JavaGD=c("JavaGD::JavaGD")


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Drop devices this system is not capable of
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      for (type in names(res)) {
        # Assume it is supported, unless...
        supported <- TRUE

        for (fcn in res[[type]]) {
          pattern <- "^(.+)::(.+)$"
          pkg <- gsub(pattern, "\\1", fcn)
          name <- gsub(pattern, "\\2", fcn)

          if (pkg == "needs") {
            if (!supports(name)) {
              supported <- FALSE
              break
            }
            next
          }

          if (!isPackageInstalled(pkg)) {
            supported <- FALSE
            break
          }

          ns <- getNamespace(pkg)
          if (!exists(name, envir=ns, mode="function")) {
            supported <- FALSE
            break
          }

          if (!supports(type, pkg=pkg)) {
            supported <- FALSE
            break
          }
        } # for (fcn ...)

        # Drop any 'needs::*'
        res[[type]] <- grep("^needs::", res[[type]], value=TRUE, invert=TRUE)

        # Not supported?
        if (!supported) res[[type]] <- NULL
      } # for (type ...)


      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Order by name
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      o <- order(names(res))
      res <- res[o]

      # Memoize
      .devAll <<- res
    } # if (force ...)

    res
  }
}) # devAll()



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BEGIN: Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
.devList <- function() {
  if (exists(".Devices")) {
    devList <- get(".Devices");
  } else {
    devList <- list("null device");
  }

  labels <- names(devList);
  if (is.null(labels)) {
    labels <- paste("Device", seq_along(devList), sep=" ");
    names(devList) <- labels;
    assign(".Devices", devList, envir=baseenv());
  } else {
    # Update the names
    labels <- names(devList);
    idxs <- which(nchar(labels) == 0L);
    if (length(idxs) > 0L) {
      labels[idxs] <- paste("Device", idxs, sep=" ");
    }
    names(devList) <- labels;
  }

  devList;
} # .devList()

## Gets the devList() index of a device by label
.devListIndexOf <- function(labels, error=TRUE) {
  # Nothing to do?
  if (length(labels) == 0L) {
    res <- integer(0L);
    names(res) <- character(0L);
    return(res);
  }

  devList <- devList(dropNull=FALSE);
  idxs <- match(labels, names(devList));
  names(idxs) <- labels;

  # Sanity check
  if (error) {
    if (any(is.na(idxs)))
      throw("No such device: ", paste(labels[is.na(idxs)], collapse=", "));
  }

  idxs
} # .devListIndexOf()


.devNextAvailable <- function() {
  # All open devices
  devList <- dev.list();

  if (length(devList) == 0L)
    return(2L);

  devPossible <- seq(from=2L, to=max(devList)+1L);
  devFree <- setdiff(devPossible, devList);

  devFree[1L];
} # .devNextAvailable()


.devTypeName <- function(types, pattern=FALSE, knownTypes=names(devAll()), ...) {
  # Nothing todo?
  if (!is.character(types)) {
    return(types);
  }

  names <- names(types);
  if (is.null(names)) {
    names(types) <- types;
  }

  # Match to known set of device types by regular expression?
  if (pattern) {
    types <- as.list(types);
    for (kk in seq_along(types)) {
      typeKK <- types[[kk]];
      pattern <- sprintf("^%s$", typeKK);
      idxs <- grep(pattern, knownTypes);
      if (length(idxs) > 0L) {
        typesKK <- knownTypes[idxs];
        names(typesKK) <- rep(typeKK, times=length(typesKK));
        types[[kk]] <- typesKK;
      } else {
        names(types[[kk]]) <- typeKK;
      }
    } # for (kk ...)
    types <- unlist(types, use.names=TRUE);
  }

  # Common aliases
  names <- names(types);
  types[types == "jpg"] <- "jpeg";
  types[types == "ps"] <- "postscript";
  names(types) <- names;

  types;
} # .devTypeName()

.devTypeExt <- function(types, ...) {
  if (is.function(types)) {
    types <- .devTypeNameFromFunction(types)
  }
  types <- as.character(types)
  exts <- types

  ## Cairo package
  pattern <- "^Cairo(JPEG|PDF|PNG|PS|SVG|TIFF)$";
  idxs <- grep(pattern, exts);
  exts[idxs] <- tolower(gsub(pattern, "\\1", exts[idxs]));

  ## cairo_* devices
  pattern <- "^cairo_(pdf|ps)$";
  exts <- gsub(pattern, "\\1", exts);

  ## Recognize types of any case. Always return in lower case.
  exts <- tolower(exts);

  # Common type-to-extension conversions
  exts[exts == "win.metafile"] <- "wmf";
  exts[exts == "png2"] <- "png";
  exts[exts == "jpeg"] <- "jpg";
  exts[exts == "jpeg2"] <- "jpg";
  exts[exts == "postscript"] <- "ps";

  exts
} # .devTypeExt()


.devTypeNameFromFunction <- function(fcn, knownTypes=R.devices:::devAll(), ...) {
  stopifnot(length(fcn) == 1, is.function(fcn))
  knownFcns <- lapply(knownTypes, FUN=function(x) eval(parse(text=x[1])))
  name <- NULL
  for (name in names(knownFcns)) {
    if (identical(fcn, knownFcns[[name]])) return(name)
  }
  NA_character_
} # .devTypeNameFromFunction()


.devEqualTypes <- (function() {
   # Recorded known results
   known <- list();

   function(type, other, args=list()) {
     # (a) Same types?
     if (identical(unname(type), unname(other))) return(TRUE);

     # (b) A known equality?
     if (is.character(type) && is.character(other)) {
       res <- known[[type]][other];
       if (is.logical(res)) return(res);
     }

     # (c) Comparing to a device function?
     if (is.function(type) && is.character(other)) {
       if (!exists(other, mode="function")) return(FALSE);
       otherT <- get(other, mode="function");
       if (identical(unname(type), unname(otherT))) return(TRUE);
     } else if (is.function(other) && is.character(type)) {
       if (!exists(type, mode="function")) return(FALSE);
       typeT <- get(type, mode="function");
       if (identical(unname(typeT), unname(other))) return(TRUE);
     }

     # (d) Check if temporarily opening the requested type actually
     # creates a device matching the existing one.
     typeT <- tryCatch({
       local({
         do.call(type, args=args)
         on.exit(dev.off(), add=TRUE)
       })
       names(dev.cur())
     }, error=FALSE)
     if (identical(typeT, FALSE)) return(FALSE)

     if (is.function(other)) {
       other <- tryCatch({
         local({
           do.call(other, args=args)
           on.exit(dev.off(), add=TRUE)
           names(dev.cur())
         })
       }, error=FALSE)
       if (identical(other, FALSE)) return(FALSE)
     }

     res <- (typeT == other);

     # Record result to avoid opening next?
     if (is.character(type) && is.character(other)) {
       known[[type]][other] <<- res;
     }

     res;
  }
})() # .devEqualTypes()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# END: Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


############################################################################
# HISTORY:
# 2015-12-15
# o CLARIFICATION: Renamed .devIndexOf() to .devListIndexOf().
# 2014-10-17
# o SPEEDUP: Made devAll() memoize base::capabilities() results, which
#   gives a significant speedup when for instance an X11 server times out.
# 2014-09-30
# o BUG FIX: devDone() would close some devices despite them being
#   on screen/interactive devices, e.g. an x11 device.
# 2014-09-17
# o Added devAll().
# 2014-09-16
# o Now devIsInteractive() without arguments returns all known
#   interactive devices.
# 2014-09-11
# o Added .devTypeExt().
# 2014-04-27
# o Added .devEqualTypes().
# 2013-10-29
# o ROBUSTESS/BUG FIX: devSet(which) where 'which' is a very large number
#   could leave lots of stray temporary devices open when error "too many
#   open devices" occurred.  Now all temporary devices are guaranteed to
#   be closed also when there is an error.
# 2013-10-28
# o BUG FIX: dev(Get|Set)Label(which) would not handle the case when
#   the device specified by an numeric 'which' and there is a
#   gap in the device list.
# o Added argument 'interactiveOnly' to devList().
# o ROBUSTNESS: Now devSet() is guaranteed to close all temporary
#   devices it opens.
# 2013-10-15
# o BUG FIX: devSet(key), where 'key' is a non-integer object (which is
#   coerced to a device label via digest()), stopped working due to a too
#   conservative test.
# 2013-09-24
# o CONSISTENCY: Now devList() returns an empty integer vector
#   (instead of NULL) if no open devices exists.
# o Now devOff() and devDone() checks if device is opened before trying
#   to close it.  This avoids opening and closing of non-opened devices.
# o ROBUSTNESS: The device functions that are not vectorize do now
#   throw an informative error if passed a vector.
# o GENERALIZATION: Vectorized devIsOpen(), devGetLabel(), and
#   devIsInteractive().
# o Vectorized internal .devIndexOf().
# o Added argument 'dropNull' to devList().
# 2013-08-27
# o Added devIsInteractive().
# o Added .devTypeName().
# 2012-11-18
# o Replaced all stop() with throw().
# 2012-04-30
# o Moved devNew() to devNew.R.
# o Moved devEval() to devEval.R.
# 2011-11-05
# o Now the default 'width' is inferred from devOptions() is 'height'
#   is not given and aspectRatio != 1.
# 2011-03-16
# o Now R.archive:ing is only done if the R.archive package is loaded.
# o DOCUMENTATION: The title of devDone() was incorrect.
# 2008-10-26
# o Now argument 'which' to devSet() can be any object.  If not a single
#   numeric or a single character string, then a checksum character string
#   is generated using digest::digest(which).
# 2008-10-16
# o Now devDone(which=1) does nothing.  Before it gave an error.
# 2008-08-01
# o Added devList() and removed devLabels().
# o Added internal .devNextAvailable().
# o Added argument 'error=TRUE' to internal .devIndexOf().
# 2008-07-31
# o Now devSet(idx) opens a new device with index 'idx' if not already
#   opened.
# 2008-07-29
# o Using term 'label' instead of 'name' everywhere, e.g. devLabels().
#   This was changed because the help pages on 'dev.list' etc. already
#   use the term 'name' for a different purpose, e.g. 'windows'.
# 2008-07-18
# o Created.
############################################################################
