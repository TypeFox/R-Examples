###########################################################################/**
# @RdocFunction devEval
# @alias devDump
#
# @title "Opens a new graphics device, evaluate (graphing) code, and closes device"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{type}{Specifies the type of graphics device to be used.
#    The device is created and opened using @see "devNew".
#    Multiple types may be specified.}
#   \item{expr}{The @expression of graphing commands to be evaluated.}
#   \item{initially, finally}{Optional @expression:s to be evaluated
#    before and after \code{expr}. If \code{type} specifies multiple
#    devices, these optional @expression:s are only evaluated ones.}
#   \item{envir}{The @environment where \code{expr} should be evaluated.}
#   \item{name, tags, sep}{The fullname name of the image is specified
#     as the name with optional \code{sep}-separated tags appended.}
#   \item{ext}{The filename extension of the image file generated, if any.
#    By default, it is inferred from argument \code{type}.}
#   \item{...}{Additional arguments passed to @see "devNew".}
#   \item{filename}{The filename of the image saved, if any.
#    By default, it is composed of arguments \code{name}, \code{tags},
#    \code{sep}, and \code{ext}.  See also below.}
#   \item{path}{The directory where then image should be saved, if any.}
#   \item{field}{An optional @character string specifying a specific
#     field of the named result @list to be returned.}
#   \item{onIncomplete}{A @character string specifying what to do with
#     an image file that was incompletely generated due to an interrupt
#     or an error.}
#   \item{force}{If @TRUE, and the image file already exists, then it is
#     overwritten, otherwise not.}
#   \item{which}{A @vector of devices to be copied.  Only applied if
#     argument \code{expr} is missing.}
#   \item{.exprAsIs, .allowUnknownArgs}{(Internal use only).}
# }
#
# \value{
#   Returns a @see "DevEvalFileProduct" if the device generated an
#   image file, otherwise an @see "DevEvalProduct".
#   If argument \code{field} is given, then the field of the
#   @see "DevEvalProduct" is returned instead.
#   \emph{Note that the default return value may be changed in future releases.}
# }
#
# \section{Generated image file}{
#   If created, the generated image file is saved in the directory
#   specfied by argument \code{path} with a filename consisting of
#   the \code{name} followed by optional comma-separated \code{tags}
#   and a filename extension given by argument \code{ext}.
#
#   By default, the image file is only created if the \code{expr}
#   is evaluated completely.  If it is, for instance, interrupted
#   by the user or due to an error, then any incomplete/blank image
#   file that was created will be removed.  This behavior can be
#   turned of using argument \code{onIncomplete}.
# }
#
# @examples "../incl/devEval.Rex"
#
# @author
#
# \seealso{
#   To change default device parameters such as the width or the height,
#   @see "devOptions".
#   @see "devNew".
# }
#
# @keyword device
# @keyword utilities
#*/###########################################################################
devEval <- function(type=getOption("device"), expr, initially=NULL, finally=NULL, envir=parent.frame(), name=NULL, tags=NULL, sep=getDevOption(type, "sep", default=","), ..., ext=NULL, filename=NULL, path=getDevOption(type, "path", default="figures"), field=getDevOption(type, name="field"), onIncomplete=c("remove", "rename", "keep"), force=getDevOption(type, "force", default=TRUE), which=dev.cur(), .exprAsIs=FALSE, .allowUnknownArgs=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Vectorized version
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Nothing to do?
  if (length(type) == 0L) {
    return(structure(character(0L), class=class(DevEvalProduct())));
  }

  # Parse 'type' in case multiple types is specified in one string
  if (is.character(type)) {
    type <- unlist(strsplit(type, split=","), use.names=FALSE);
    type <- trim(type);
  }

  # Argument 'expr':
  hasExpr <- !missing(expr);

  # Argument '.allowUnknownArgs':
  ## If TRUE, it tells devNew() that it's ok to silently ignore
  ## arguments not recognized by certain device functions.
  ## Used when multiple device types are used in one devEval() call.
  .allowUnknownArgs <- as.logical(.allowUnknownArgs);


  ## Expand device type names by regexp matching, iff any
  type <- .devTypeName(type, pattern=TRUE);


  ## SPECIAL CASE: Handle calls like toNnn({ expr }) where the actual plot
  ## expression is actually passed via the first argument which is 'name'
  ## and not via 'expr', e.g. toX11({ plot(1:3) }).
  if (!hasExpr) {
    # Was the expression passed implicitly via 'name' instead?
    # In order to infer this, we have to inspect 'name' in the
    # parent frame:
    nameClass <- eval(expression(class(substitute(name))), envir=parent.frame());
    isExpr <- is.element(nameClass, c("call", "{", "("))
    if (isTRUE(isExpr)) {
      # This avoid the plot expression to be evaluated here.
      delayedAssign("expr", name);
      hasExpr <- TRUE;
      nameOrg <- NULL;
    } else {
      nameOrg <- name;
    }
  } else {
    nameOrg <- name;
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (a) No plot code expression?  Then copy multiple input devices...
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!hasExpr && length(which) > 1L) {
    res <- list();
    for (kk in seq_along(which)) {
      idx <- which[kk];
      devSet(idx);
      res[[kk]] <- devEval(type=type, initially=NULL, finally=NULL, envir=envir, name=nameOrg, tags=tags, sep=sep, ..., ext=ext, filename=filename, path=path, field=field, onIncomplete=onIncomplete, force=force, which=idx, .allowUnknownArgs=TRUE);
    } # for (idx ...)
    names(res) <- names(which);
    return(res);
  } # if (length(which) > 1L)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (b) Plot code expression, but with multiple output types/devices?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (length(type) > 1L) {
    types <- type;

    # Evaluate 'initially' only once
    eval(initially, envir=envir);

    if (hasExpr) {
      # Expression must be substitute():d to avoid the being evaluated here
      if (!.exprAsIs) expr <- substitute(expr);

      # Evaluate 'expr' once per graphics device
      res <- lapply(types, FUN=function(type) {
        devEval(type=type, expr=expr, initially=NULL, finally=NULL, envir=envir, name=nameOrg, tags=tags, sep=sep, ..., ext=ext, filename=filename, path=path, field=field, onIncomplete=onIncomplete, force=force, .exprAsIs=TRUE, .allowUnknownArgs=TRUE);
      });
    } else {
      # Evaluate 'expr' once per output device
      res <- lapply(types, FUN=function(type) {
        devEval(type=type,            initially=NULL, finally=NULL, envir=envir, name=nameOrg, tags=tags, sep=sep, ..., ext=ext, filename=filename, path=path, field=field, onIncomplete=onIncomplete, force=force, which=which, .allowUnknownArgs=TRUE);
      });
    } # if (hasExpr)
    names(res) <- types;

    # Evaluate 'finally' only once
    eval(finally, envir=envir);

    return(res);
  } # if (length(type) > 1L)

  # Sanity check
  stopifnot(length(type) == 1L);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (c) Use first successful device type among several options?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  types <- type
  if (is.character(types)) {
    types <- unlist(strsplit(types, split="|", fixed=TRUE));
    types <- trim(types);
    types <- unique(types);
  }

  if (length(types) > 1L) {
    res <- NULL;
    errors <- list(); # Record any errors

    # Evaluate 'initially' only once
    eval(initially, envir=envir);

    if (hasExpr) {
      # Expression must be substitute():d to avoid the being evaluated here
      if (!.exprAsIs) expr <- substitute(expr);

      # Evaluate 'expr' once per graphics device
      for (type in types) {
        tryCatch({
          res <- devEval(type=type, expr=expr, initially=NULL, finally=NULL, envir=envir, name=nameOrg, tags=tags, sep=sep, ..., ext=ext, filename=filename, path=path, field=field, onIncomplete=onIncomplete, force=force, .exprAsIs=TRUE, .allowUnknownArgs=TRUE);
          attr(res, "type") <- type;
          break
        }, error = function(ex) {
          errors <<- c(errors, list(ex))
        })
      } # for (type ...)
    } else {
      # Evaluate 'expr' once per output device
      for (type in types) {
        tryCatch({
          res <- devEval(type=type,            initially=NULL, finally=NULL, envir=envir, name=nameOrg, tags=tags, sep=sep, ..., ext=ext, filename=filename, path=path, field=field, onIncomplete=onIncomplete, force=force, which=which, .allowUnknownArgs=TRUE);
          attr(res, "type") <- type;
          break
        }, error = function(ex) {
          errors <<- c(errors, list(ex))
        })
      } # for (type ...)
    } # if (hasExpr)

    # Evaluate 'finally' only once
    eval(finally, envir=envir);

    # Successfully created output?
    if (!is.null(res)) {
      return(res);
    }

    msg <- sprintf("None of the specified device types (%s) generated output", paste(sQuote(types), collapse=", "));
    if (length(errors) > 0L) {
      reasons <- sapply(errors, FUN=function(ex) ex$message);
      reasons <- unique(reasons);
      reasons <- sprintf("(%d) %s", seq_along(reasons), reasons);
      reasons <- paste(reasons, collapse="; ");
      msg <- sprintf("%s. The reasons reported were: %s", msg, reasons)
    }
    throw(msg);
  } # if (length(type) > 1L)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (d) Plot code expression and a single output type/device
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sanity check
  stopifnot(length(type) == 1L);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'type':
  # Sanity check
  if (length(type) != 1L) {
    throw("Argument 'type' must be a single object: ", length(type));
  }
  if (is.function(type)) {
  } else {
    type <- as.character(type);
    type <- .devTypeName(type);
  }

  # An interactive/non-file device?
  isInteractive <- devIsInteractive(type);

  # Argument 'name', 'tags' and 'sep':
  # Default 'name' value
  if (is.null(nameOrg) && !isInteractive) {
    if (hasExpr) {
      nameOrg <- "Rplot"; # Backward compatible
    } else {
      nameOrg <- devGetLabel(which);
    }
  }

  # Construct the full name
  fullname <- paste(c(nameOrg, tags), collapse=sep);
  fullname <- unlist(strsplit(fullname, split=sep, fixed=TRUE));
  fullname <- sub("^[\t\n\f\r ]*", "", fullname); # trim tags
  fullname <- sub("[\t\n\f\r ]*$", "", fullname); #
  fullname <- fullname[nchar(fullname) > 0L];     # drop empty tags
  fullname <- paste(fullname, collapse=sep);

  if (!isInteractive) {
    # Argument 'ext':
    if (is.null(ext)) {
      # Infer filename extension from the type
      if (is.character(type)) {
        ext <- .devTypeExt(type);
      } else {
        ext <- .devTypeExt(substitute(type));
      }
    }

    # Argument 'filename' & 'path':
    if (is.null(filename)) {
      filename <- sprintf("%s.%s", fullname, ext);
    }
    pathname <- Arguments$getWritablePathname(filename, path=path);
  }


  # Argument 'field':
  if (!is.null(field)) {
    field <- Arguments$getCharacter(field);
  }

  # Argument 'onIncomplete':
  onIncomplete <- match.arg(onIncomplete);

  # Argument 'force':
  force <- Arguments$getLogical(force);

  devIdx <- NULL;
  if (!isInteractive) {
    # Make sure the currently open device, iff any, is still the active
    # one when returning from this function.
    devListEntry <- devList();
    devCur <- dev.cur();
    on.exit({
      if (devCur != 1L) devSet(devCur);

      # Assert that no temporarily opened devices are left behind,
      devListExit <- setdiff(devList(), devIdx);
      devListDiff <- setdiff(devListExit, devListEntry);
      if (length(devListDiff) > 0L) {
        types <- unlist(.devList());
        types <- types[devListDiff];
        throw("Detected new graphics devices that was opened but not closed while executing devEval(): ", paste(sprintf("%s (%s)", sQuote(devListDiff), types), collapse=", "));
      }
    }, add=TRUE)
  }

  # Result object
  if (isInteractive) {
    res <- DevEvalProduct(fullname, type=type);
  } else {
    res <- DevEvalFileProduct(pathname, type=type);
  }

  done <- FALSE;
  if (force || isInteractive || !isFile(pathname)) {
    # Open device
    if (isInteractive) {
      devIdx <- devNew(type, which=fullname, ..., .allowUnknownArgs=.allowUnknownArgs);
    } else {
      devIdx <- devNew(type, pathname, ..., .allowUnknownArgs=.allowUnknownArgs);

      on.exit({
        # Make sure to close the device (the same that was opened),
        # unless it's an interactive device.
        if (!is.null(devIdx)) {
          devDone(devIdx);
        }

        # Archive file?
        if (isPackageLoaded("R.archive")) {
          # To please R CMD check
          getArchiveOption <- archiveFile <- NULL;
          if (getArchiveOption("devEval", FALSE)) archiveFile(pathname);
        }

        if (!done && isFile(pathname)) {
          if (onIncomplete == "remove") {
            # Remove incomplete image file
            file.remove(pathname);
          } else if (onIncomplete == "rename") {
            # Rename incomplete image file to denote that
            # it is incomplete.
            pattern <- "^(.*)[.]([^.]*)$";
            if (regexpr(pattern, pathname) != -1L) {
              fullname <- gsub(pattern, "\\1", pathname);
              ext <- gsub(pattern, "\\2", pathname);
              fmtstr <- sprintf("%s,INCOMPLETE_%%03d.%s", fullname, ext);
            } else {
              fmtstr <- sprintf("%s,INCOMPLETE_%%03d", pathname);
            }

            # Try to rename
            maxTries <- 999L
            for (kk in seq_len(maxTries)) {
              pathnameN <- sprintf(fmtstr, kk);
              if (isFile(pathnameN)) next;
              resT <- file.rename(pathname, pathnameN);
              # Done?
              if (resT && isFile(pathnameN)) {
                pathname <- pathnameN;
                break;
              }
            } # for (kk ...)

            # Failed to rename?
            if (!isFile(pathnameN) && isFile(pathname)) {
              throw(sprintf("Failed to rename incomplete image file (after %d tries): %s", pathname, kk));
            }
          } # if (onIncomplete == ...)
        } # if (!done && isFile(...))
      }, add=TRUE); # on.exit()
    } # if (isInteractive)

    # Evaluate 'initially', 'expr' and 'finally' (in that order)
    eval(initially, envir=envir);
    if (hasExpr) {
      eval(expr, envir=envir);
    } else {
      # Copy the device specified by argument 'which'.
      # Assert that device is interactive
      devList <- devList(interactiveOnly=TRUE);
      if (is.numeric(which)) {
        ok <- is.element(which, devList);
      } else {
        ok <- is.element(which, names(devList));
      }
      if (!ok) {
        types <- unlist(.devList());
        type <- types[which];
        throw(sprintf("Cannot copy a %s device - only interactive/screen devices are supported: %s", sQuote(type), sQuote(which)));
      }

      devSet(which);
      dev.copy(which=devIdx);
    }
    eval(finally, envir=envir);

    done <- TRUE;
  }

  # Close it here to make sure the image file is created.
  # This is needed if field="dataURI" (which triggers a
  # reading the image file).
  if (done) {
    devDone(devIdx);
    devIdx <- NULL;
  }

  # Subset?
  if (!isInteractive && !is.null(field)) {
    res <- res[[field]];
  }

  res;
} # devEval()


devDump <- function(type=c("png", "pdf"), ..., path=NULL, envir=parent.frame(), field=NULL, which=devList(interactiveOnly=TRUE)) {
  if (is.null(path)) {
    # Timestamp, e.g. 2011-03-10_041359.032
    timestamp <- format(Sys.time(), "%Y-%m-%d_%H%M%OS3", tz="");
    path <- getDevOption("*", "path", default="figures")
    path <- file.path(path, timestamp);
  }

  devEval(type=type, ..., path=path, envir=envir, field=field, which=which);
} # devDump()


############################################################################
# HISTORY:
# 2015-02-01
# o ROBUSTNESS: When calling devEval() on multiple devices then
#   underlying graphics device functions no longer gives an error
#   of unknown arguments, which means it is now possible to pass
#   device specific arguments also when in multi-device calls, e.g.
#   devEval(c("x11", "png"), res=100, { plot(1:10) }).
# 2015-01-21
# o BUG FIX: Renaming incomplete files gave an error if there were already
#   files renamed for similar reasons.
# o BUG FIX: devEval(type=<device function>, ..., ext) did not work.
# 2014-09-16
# o Added support regexpr device type names.
# 2014-09-12
# o Now arguments 'sep', 'path', 'field' and 'force' uses devOptions("*").
#   If old-style R options are set they are used first.
# 2014-09-11
# o Now devEval(..., ext=NULL) does a better job of inferring the default
#   filename extension from the device type.
# 2014-08-29
# o Added support for devEval() to try multiple device types one-by-one
#   until success, e.g. devEval("png|jpg|bmp", ...).
# 2014-04-28
# o By code inspection, it's now possible to do toX11({ plot(1:2) }), where
#   argument 'name' is not given and the plot expression is the first
#   argument (which will be assigned to 'name').  Internally, this is
#   handled by devEval(), instead of each single toNnn() function.
# 2014-04-27
# o BUG FIX: Now devEval("windows", { plot(1:10) }) no longer gives
#   "Error: Detected new graphics devices that was opened but not
#   closed while executing devEval(): '2' (windows)".
# 2014-01-02
# o Now the timestamp of the default path for devDump() is in the
#   current time zone.
# 2013-10-29
# o Now devDump() by default outputs to figures/<timestamp>/
# 2013-10-28
# o ROBUSTNESS: Now devEval() detects and gives an informative error
#   if one tries to copy a non-interactive/non-screen device.
# o Added devDump() which is short for devEval(c("png", "pdf"), ...,
#   which=devList()).
# o If 'expr' is missing, devEval() copies the current active device
#   and devEval(which=devList()) copies all open devices.
# 2013-09-27
# o BUG FIX: devEval() could generate "Error in devEval(type = "...",
#   name = name, ..., field = field) : object 'done' not found".
# 2013-09-25
# o Added arguments 'initially' and 'finally'.
# o Vectorized devEval().
# o Updated the formal defaults of several devEval() arguments to be NULL.
#   Instead, NULL for such arguments are translated to default internally.
#   This makes it easier/possible to vectorize devEval().
# 2013-09-24
# o ROBUSTNESS: Now devEval() gives an error if argument 'type' is not
#   of length one.  However, eventually devEval() will be vectorized.
# 2013-08-27
# o Now devEval() utilizes devIsInteractive().
# 2013-08-17
# o BUG FIX/ROBUSTNESS: Argument 'ext' of devEval() can now be inferred
#   from argument 'type' also when 'type' is passed via a string variable.
# 2013-07-29
# o Now devEval() returns a list of class 'DevEvalFile'.
# 2013-07-15
# o Added argument 'sep' to devEval() together with an option to set
#   its default value.
# 2013-04-04
# o ROBUSTNESS: Now devEval() does a better job of making sure to close
#   the same device as it opened.  Previously it would close the current
#   active device, which would not be the correct one if for instance
#   other devices had been open in the meanwhile/in parallel.
# 2013-02-23
# o Now argument 'field' for devEval() defaults to
#   getOption("devEval/args/field", NULL).
# 2012-08-21
# o DOCUMENTATION: Link to devNew() help was broken.
# 2012-04-30
# o Extracted devEval() to devEval.R.
# 2012-04-05
# o Now devEval() returns the correct pathname if onIncomplete="rename"
#   and image file was incomplete.  Now it also gives an error, if it
#   for some reason fails to rename the file in such cases.
# 2012-04-04
# o Now it is possible to have devEval() rename incompletely generated
#   image files, by using argument onIncomplete="rename".
# 2012-02-28
# o Added argument 'onIncomplete' to devEval().  The default is now to
#   remove any half-generated image files so that no incomplete/blank
#   image files are left behind.
# 2011-10-31
# o Added argument 'field' to devEval().
# 2011-04-12
# o Now devEval("jpg", ...) is recognized as devEval("jpeg", ...).
# 2011-03-29
# o Now argument 'force' of devEval() defaults to
#   getOption("devEval/args/force", TRUE).
# 2011-03-18
# o Now devEval() does a better job of "cleaning up" 'name' and 'tags'.
# o Now argument 'path' of devEval() defaults to
#   getOption("devEval/args/path", "figures/").
# 2011-03-16
# o Now R.archive:ing is only done if the R.archive package is loaded.
# 2011-03-09
# o Added support for automatic file archiving in devEval().
# 2011-02-20
# o Changed argument 'force' of devEval() to default to TRUE.
# o Added argument 'par' to devNew() allowing for applying graphical
#   settings at same time as the device is opened, which is especially
#   useful when using devEval().
# 2011-02-14
# o Now devEval() returns a named list.
# 2011-02-13
# o Added devEval().
############################################################################
