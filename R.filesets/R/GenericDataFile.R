###########################################################################/**
# @RdocClass GenericDataFile
#
# @title "The abstract GenericDataFile class"
#
# \description{
#  @classhierarchy
#
#  A GenericDataFile is an object refering to a data file on a file system.
#  Note that this class is abstract and can not be instanciated, but
#  instead you have to use one of the subclasses or the generic
#  @seemethod "fromFile" method.
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{The filename of the file.}
#   \item{path}{An optional path to the file.}
#   \item{mustExist}{If @TRUE, an exception is thrown if the file does
#     not exists, otherwise not.}
#   \item{...}{Not used.}
#   \item{.onUnknownArgs}{A @character string specifying what should occur
#      if there are unknown arguments in \code{...}.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Filename convention}{
#   The filename of an \code{GenericDataFile} is structured as follows:
#   \itemize{
#    \item{filename}{: \code{"sample001,a,b,c.CEL"}
#       (this follows the \R convention, but not the Unix convention)}
#    \item{fullname}{: \code{"sample001,a,b,c"}}
#    \item{name}{: \code{"sample001"}}
#    \item{tags}{: \code{c("a", "b", "c")}}
#    \item{extension}{: \code{"CEL"}}
#   }
# }
#
# @author
#
# \seealso{
#   An object of this class is typically part of an @see "GenericDataFileSet".
# }
#*/###########################################################################
setConstructorS3("GenericDataFile", function(filename=NULL, path=NULL, mustExist=!is.na(filename), ..., .onUnknownArgs=c("error", "warning", "ignore")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(filename)) {
    ## Developers mode: Test default NA_character_ early on. /HB 2015-12-23
    filename <- getOption("R.filesets/GenericDataFile/args/filename", NULL)
    ## FIXME: In a future release, make NULL an invalid value. /HB 2015-12-23
    ## throw("Argument 'filename' must be a single character string: NULL")
  } else if (length(filename) != 1L) {
    throw("Argument 'filename' must be a single character string: ", length(filename))
  }

  if (!is.null(filename)) {
    pathname <- Arguments$getReadablePathname(filename, path=path, absolutePath=TRUE, mustExist=mustExist)
    # Assert that it is not pointing to a directory
    if (isDirectory(pathname)) {
      throw("The specified pathname is a directory: ", pathname)
    }
  } else {
    pathname <- NULL
  }

  # Arguments '...':
  args <- list(...);

  # Ignore any argument called 'recursive'
  keep <- which(regexpr("^recursive$", names(args)) == -1L);
  args <- args[keep];

  if (length(args) > 0L) {
    if (is.element(.onUnknownArgs, c("error", "warning"))) {
      argsStr <- paste(names(args), collapse=", ");
      msg <- sprintf("Unknown arguments: %s", argsStr);
      if (.onUnknownArgs == "error") {
        throw(msg);
      } else if (.onUnknownArgs == "warning") {
        warning(msg);
      }
    }
  }

  this <- extend(Object(), c("GenericDataFile", uses("FullNameInterface")),
    "cached:.checksum" = NULL,
    .prevModifiedOn = NULL,
    .pathname = pathname,
    .attributes = list()
  );

  # Update time stamps
  hasBeenModified(this);

  this;
}, abstract=TRUE)



setMethodS3("clone", "GenericDataFile", function(this, clear=TRUE, ...) {
  object <- NextMethod("clone");
  if (clear)
    clearCache(object);
  object;
}, protected=TRUE)



###########################################################################/**
# @RdocMethod equals
# @alias equals
#
# @title "Checks if a file equals another"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{other}{The other @see "GenericDataFile" to be compared to.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns @TRUE if the file equals the other, otherwise @FALSE.
#   If @FALSE, attributes are added specifying the pathnames of the two
#   files compared, and the reason for them being different.
# }
#
# \details{
#   The two files compared are equal if they have the same pathname.
#
#   The two files compared are \emph{not} equal if:
#   \itemize{
#    \item Argument \code{other} is not a @see "GenericDataFile", or
#    \item their file sizes differ, or
#    \item their file checksums differ.
#   }
#
#   If none of the above occurs, the two files are considered equal.
#
#   Note that subclasses use refined rules.
# }
#
# @author
#
# \seealso{
#   @seemethod "getFileSize".
#   @seemethod "getChecksum".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("equals", "GenericDataFile", function(this, other, ...) {
  # Default values
  notEqual <- FALSE;
  attr(notEqual, "thisFile") <- getPathname(this);
  if (!inherits(other, "GenericDataFile")) {
    msg <- sprintf("The 'other' is not a GenericDataFile: %s",
                                                 class(other)[1]);
    attr(notEqual, "reason") <- msg;
    return(notEqual);
  }

  attr(notEqual, "otherFile") <- getPathname(other);
  msg <- NULL;

  if (identical(getPathname(this), getPathname(other)))
    return(TRUE);

  value <- getFileSize(this);
  valueOther <- getFileSize(other);
  if (!identical(value, valueOther)) {
    msg <- sprintf("The file sizes differ: %.0f != %.0f",
                                          value, valueOther);
    attr(notEqual, "reason") <- msg;
    return(notEqual);
  }

  value <- getChecksum(this);
  valueOther <- getChecksum(other);
  if (!identical(value, valueOther)) {
    msg <- sprintf("The checksums differ: %s != %s",
                                          value, valueOther);
    attr(notEqual, "reason") <- msg;
    return(notEqual);
  }

  TRUE;
})



###########################################################################/**
# @RdocMethod as.character
#
# @title "Returns a short string describing the file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("as.character", "GenericDataFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);

  # Name and tags of file set
  s <- c(s, sprintf("Name: %s", getName(this)));
  tags <- getTags(this, collapse=",");
  if (!is.null(tags)) {
    s <- c(s, sprintf("Tags: %s", tags));
  }

  # Full names of file set
  s <- c(s, sprintf("Full name: %s", getFullName(this)));

  # Pathname
  pathname <- getPathname(this, absolute=FALSE, NULLasNA=TRUE)
  if (!is.na(pathname)) {
    pathnameA <- getPathname(this, absolute=TRUE);
    if (nchar(pathnameA, type="chars") < nchar(pathname, type="chars")) {
      pathname <- pathnameA;
    }
  }
  s <- c(s, paste("Pathname: ", pathname, sep=""));

  # File size
  fileSize <- getFileSize(this, "units");
  if (!is.na(fileSize)) {
    fileSizeB <- sprintf("%.0f bytes", getFileSize(this, "numeric"));
    if (fileSizeB != fileSize) {
      fileSize <- sprintf("%s (%s)", fileSize, fileSizeB);
    }
  }
  s <- c(s, sprintf("File size: %s", fileSize));

  s <- c(s, sprintf("RAM: %.2f MB", objectSize(this)/1024^2));

  GenericSummary(s);
}, protected=TRUE)




###########################################################################/**
# @RdocMethod getPathname
#
# @title "Gets the pathname of the file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{absolute}{If @TRUE, the absolute pathname is returned,
#     otherwise the relative.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns the pathname as @character string
#   (or @NULL if an "empty" file).
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getPathname", "GenericDataFile", function(this, absolute=FALSE, ...) {
  ## Secret argument. /HB 2015-12-23
  args <- list(...)
  NULLasNA <- isTRUE(args$NULLasNA)

  pathname <- this$.pathname
  if (is.null(pathname)) {
    if (NULLasNA) pathname <- NA_character_
  } else {
    if (absolute) {
      pathname <- getAbsolutePath(pathname)
    } else {
      pathname <- getRelativePath(pathname)
    }
  }
  pathname
})



###########################################################################/**
# @RdocMethod getPath
#
# @title "Gets the path (directory) of the file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Optional arguments passed to @seemethod "getPathname".}
# }
#
# \value{
#   Returns a @character string (@NA if an "empty" file).
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getPath", "GenericDataFile", function(this, ...) {
  pathname <- getPathname(this, ..., NULLasNA=TRUE)
  dirname(pathname)
})




###########################################################################/**
# @RdocMethod getFilename
# @alias getFilename
#
# @title "Gets the filename of the file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Optional arguments passed to @seemethod "getPathname".}
# }
#
# \value{
#   Returns a @character string (@NA if an "empty" file).
# }
#
# \details{
#  The filename of a file is the pathname excluding any path.
#  For instance, the filename of \code{path/to/foo,a.2,b.ext} is
#  \code{foo,a.2,b.ext}.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFilename", "GenericDataFile", function(this, ...) {
  pathname <- getPathname(this, NULLasNA=TRUE)
  basename(pathname)
})



###########################################################################/**
# @RdocMethod getDefaultFullName
#
# @title "Gets the full name of the file"
#
# \description{
#   @get "title", that is the filename without the filename extension.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character string (@NA if "empty" file).
# }
#
# \details{
#  The full name of a file is the filename excluding any
#  extension (and period).
#  For instance, the full name of \code{path/to/foo,a.2,b.ext} is
#  \code{foo,a.2,b}.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getDefaultFullName", "GenericDataFile", function(this, ...) {
  # TO BE REMOVED. /HB 2013-10-15
  args <- list(...);
  if (is.element("aliased", names(args))) {
    .Defunct(msg="Argument 'aliased' of getDefaultFullName() for GenericDataFile is deprecated.  Use fullname translators instead.");
  }

  filename <- getFilename(this, ...);
  if (is.null(filename))
    return(as.character(NA));
  pattern <- getExtensionPattern(this);
  fullname <- gsub(pattern, "", filename);

  fullname;
}, protected=TRUE)


setMethodS3("getFilenameExtension", "GenericDataFile", abstract=TRUE, protected=TRUE)


setMethodS3("getOutputExtension", "GenericDataFile", function(...) {
  getFilenameExtension(...);
}, protected=TRUE)


setMethodS3("getExtensionPattern", "GenericDataFile", function(this, ..., default="\\.([^.]+)$", force=FALSE) {
  pattern <- this$.extensionPattern;
  if (force || is.null(pattern)) {
    # Argument 'default':
    if (!is.null(default)) {
      default <- Arguments$getRegularExpression(default);
    }

    # Default pattern is anything following the last period
    pattern <- default;

##   pattern <- toAsciiRegExprPattern(pattern); # Don't handle [.] and ~
    this$.extensionPattern <- pattern;
  }
  pattern;
}, static=TRUE, protected=TRUE)


setMethodS3("setExtensionPattern", "GenericDataFile", function(this, pattern=NULL, ...) {
  # Argument 'pattern':
  if (!is.null(pattern)) {
    pattern <- Arguments$getRegularExpression(pattern);
  }

  this$.extensionPattern <- pattern;
  invisible(this);
}, protected=TRUE)




###########################################################################/**
# @RdocMethod getExtension
#
# @title "Gets the filename extension"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to @seemethod "getFilename" and
#     @seemethod "getDefaultFullName".}
# }
#
# \value{
#   Returns a @character string (which is of length zero if "empty" file).
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getExtension", "GenericDataFile", function(this, ...) {
  filename <- getFilename(this, ...);
  fullname <- getDefaultFullName(this, ...);
  # Drop <fullname> and a possible '.'.
  substring(filename, first=nchar(fullname, type="chars") + 2L);
})




###########################################################################/**
# @RdocMethod getFileType
#
# @title "Gets the file type of a file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character in lower case letters
#   (which is of length zero if "empty" file).
# }
#
# \details{
#   By default, this methods returns the filename extension, but subclasses
#   may override this.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFileType", "GenericDataFile", function(this, ...) {
  pattern <- "(.*)[.]([abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0-9]+)$";
  filename <- getFilename(this, ...);
  ext <- gsub(pattern, "\\2", filename);
  tolower(ext);
})


###########################################################################/**
# @RdocMethod isFile
# @alias isFile
#
# @title "Checks if this is an existing file"
#
# \description{
#   @get "title" and not a directory.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @TRUE if an existing file (and not a directory),
#   otherwise @FALSE.
# }
#
# @author
#
# \seealso{
#   @see "R.utils::isFile".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("isFile", "GenericDataFile", function(this, ...) {
  pathname <- getPathname(this, NULLasNA=TRUE)
  isFile(pathname)
})


setMethodS3("is.na", "GenericDataFile", function(x) {
  pathname <- getPathname(x, NULLasNA=TRUE)
  is.na(pathname)
}, appendVarArgs=FALSE)



###########################################################################/**
# @RdocMethod validate
# @alias validate
#
# @title "Validates the content of a file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   If the file is invalid, then an error is thrown.
#   If the files is valid, then @TRUE is returned.
#   Otherwise, @NA is returned, which happens if the file was not validated.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("validate", "GenericDataFile", function(this, ...) {
  NA;
}, protected=TRUE)



###########################################################################/**
# @RdocMethod getFileSize
#
# @title "Gets the size of a file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{what}{A @character string specifying the data type returned.
#   If \code{"numeric"}, then a @numeric value is returned.
#   If \code{"units"}, then a human-readable @character string is returned.
#  }
#  \item{sep}{A @character string.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @numeric or a @character string.
#   A missing value (@NA) is returned if the file does not exist.
# }
#
# @author
#
# \seealso{
#   @see "base::file.info".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFileSize", "GenericDataFile", function(this, what=c("numeric", "units"), sep="", ...) {
  # Argument 'what':
  what <- match.arg(what);

  pathname <- getPathname(this, NULLasNA=TRUE)
  fileSize <- file.info2(pathname)$size

  if (what == "numeric")
    return(fileSize);

  if (is.na(fileSize))
    return(fileSize);

  .asIEC(fileSize)
})


###########################################################################/**
# @RdocMethod getCreatedOn
#
# @title "Gets when the file was created"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a \code{POSIXct} time stamp.
# }
#
# @author
#
# \seealso{
#   @see "base::file.info".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getCreatedOn", "GenericDataFile", function(this, ...) {
  pathname <- getPathname(this, NULLasNA=TRUE)
  file.info(pathname)$ctime
}, protected=TRUE)


###########################################################################/**
# @RdocMethod getLastModifiedOn
#
# @title "Gets when the file was last modified"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a \code{POSIXct} time stamp.
# }
#
# @author
#
# \seealso{
#   @see "base::file.info".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getLastModifiedOn", "GenericDataFile", function(this, ...) {
  pathname <- getPathname(this, NULLasNA=TRUE)
  file.info(pathname)$mtime
}, protected=TRUE)



###########################################################################/**
# @RdocMethod getLastAccessedOn
#
# @title "Gets when the file was last accessed"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a \code{POSIXct} time stamp.
# }
#
# @author
#
# \seealso{
#   @see "base::file.info".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getLastAccessedOn", "GenericDataFile", function(this, ...) {
  pathname <- getPathname(this, NULLasNA=TRUE)
  file.info(pathname)$atime
}, protected=TRUE)



###########################################################################/**
# @RdocMethod hasBeenModified
#
# @title "Checks whether the file has been modified"
#
# \description{
#   @get "title" since last time checked.
# }
#
# @synopsis
#
# \arguments{
#  \item{unknown}{The @logical value returned if the timestamp for the
#   previous modification, if any, is unknown.}
#  \item{update}{If @TRUE, the internal check timestamp is updated after
#   calling this function, otherwise not.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns @TRUE, @FALSE, or the value of argument \code{unknown}.
# }
#
# @author
#
# \seealso{
#   @see "base::file.info".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("hasBeenModified", "GenericDataFile", function(this, update=TRUE, unknown=TRUE, ...) {
  lastModifiedOn <- getLastModifiedOn(this);
  prevModifiedOn <- this$.prevModifiedOn;

  # Unknown modification timestamp on file?
  if (is.na(lastModifiedOn) || lastModifiedOn == 0L) {
    res <- unknown;
    attr(res, "lastModifiedOn") <- lastModifiedOn;
    attr(res, "prevModifiedOn") <- prevModifiedOn;
    return(res);
  }

  if (is.null(prevModifiedOn)) {
    res <- unknown;
    attr(res, "lastModifiedOn") <- lastModifiedOn;
    attr(res, "prevModifiedOn") <- prevModifiedOn;
  } else {
    res <- (lastModifiedOn > prevModifiedOn);
    attr(res, "lastModifiedOn") <- lastModifiedOn;
    attr(res, "prevModifiedOn") <- prevModifiedOn;
  }

  if (update) this$.prevModifiedOn <- lastModifiedOn;

  res;
}, protected=TRUE)



###########################################################################/**
# @RdocMethod fromFile
#
# @title "Defines a GenericDataFile from a pathname"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{filename, path}{The filename and the path to the file.  The file
#    must exist, otherwise an exception is thrown.}
#  \item{unknown}{The @logical value returned if the timestamp for the
#   previous modification, if any, is unknown.}
#  \item{...}{Not used.}
#  \item{recursive}{If TRUE, ...}
#  \item{verbose}{...}
#  \item{.checkArgs}{(Internal) If FALSE, validation of file existance and
#   arguments are skipped.}
# }
#
# \value{
#   Returns a @see "GenericDataFile" (or a subclass thereof).
# }
#
# @author
#
# \seealso{
#   @see "base::file.info".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fromFile", "GenericDataFile", function(static, filename, path=NULL, ..., recursive=TRUE, verbose=FALSE, .checkArgs=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  if (.checkArgs) {
    # Argument 'recursive':
    recursive <- Arguments$getLogical(recursive);

    # Argument 'filename' and 'path':
    pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=TRUE);
  } else {
    pathname <- filename;
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Use subclasses to setup the file?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (recursive) {
    # Get all known subclasses (bottom up)
    className <- class(static)[1L];
    clazz <- Class$forName(className, envir=parent.frame());
    knownSubclasses <- getKnownSubclasses(clazz);
    knownSubclasses <- rev(knownSubclasses);
    for (className in knownSubclasses) {
      clazz <- Class$forName(className);

      # Try reading the file using the static fromFile() method of each class
      static <- getStaticInstance(clazz);
      tryCatch({
        # Are we calling the same fromFile() instance multiple times here?
        # Slow? /HB 2010-01-30
        res <- fromFile(static, filename=pathname, .checkArgs=FALSE);
        return(res);
      }, error = function(ex) {})
    }
  }

  # If not "read" above, just create an instance as is.
  res <- newInstance(static, filename=pathname, ...);

  res;
}, static=TRUE, protected=TRUE)




###########################################################################/**
# @RdocMethod copyTo
#
# @title "Safely copies a file to a new pathname"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{filename, path}{The filename and the path for the destination file.
#   The default is to use the same filename as the source file.}
#  \item{...}{Additional arguments passed to @see "R.utils::copyFile".}
# }
#
# \value{
#   Returns a @see "GenericDataFile" (of the same class as the source file)
#   refering to the new file copy.
# }
#
# \details{
#   In order to minimize the risk for corrupt copies, the
#   @see "R.utils::copyFile" method of \pkg{R.utils} is used, which
#   provides several protection against user, system and file errors.
# }
#
# @author
#
# \seealso{
#   To link to a @see "GenericDataFile", see @seemethod "linkTo".
#   To rename a @see "GenericDataFile", see @seemethod "renameTo".
#   Internally @see "R.utils::copyFile" is used.
#   @seeclass
# }
#*/###########################################################################
setMethodS3("copyTo", "GenericDataFile", function(this, filename=getFilename(this), path=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'filename' and 'path':
  pathname <- Arguments$getWritablePathname(filename, path=path)

  # Fail-safe copying
  pathnameS <- getPathname(this, NULLasNA=TRUE)
  copyFile(pathnameS, pathname, ...)

  # Create object of the same class.
  res <- newInstance(this, pathname)

  res
}, protected=TRUE) # copyTo()



###########################################################################/**
# @RdocMethod linkTo
# @alias linkFrom
#
# @title "Creates a symbolic file link"
#
# \description{
#   @get "title" to a @see "GenericDataFile" at/from a given
#   destination pathname.
# }
#
# @synopsis
#
# \arguments{
#  \item{filename, path}{The filename and the path for the link.
#   The default is to use the same filename as the source file.}
#   \item{skip}{If @TRUE and a file with the same name as argument
#     \code{link} already exists, then the nothing is done.}
#   \item{overwrite}{If @TRUE, an existing link file is overwritten,
#     otherwise not.}
#  \item{...}{Additional arguments passed to @see "R.utils::createLink".}
# }
#
# \value{
#   Returns a @see "GenericDataFile" (of the same class as the source file)
#   refering to the file via the link.
# }
#
# \section{Required privileges on Windows}{
#   In order for this method to succeed on Windows,
#   the client/R session must run with sufficient privileges.
#   See @see "R.utils::createLink" for more details.
# }
#
# @author
#
# \seealso{
#   To copy a @see "GenericDataFile", see @seemethod "copyTo".
#   Internally @see "R.utils::createLink" is used.
#   @seeclass
# }
#*/###########################################################################
setMethodS3("linkTo", "GenericDataFile", function(this, filename=getFilename(this), path=NULL, skip=!overwrite, overwrite=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'filename' and 'path':
  pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=!skip)

  # Create link
  pathnameS <- getPathname(this, NULLasNA=TRUE)
  createLink(target=pathnameS, link=pathname, skip=skip, overwrite=overwrite, ...)

  # Create object of the same class.
  res <- newInstance(this, pathname)

  res
}, protected=TRUE) # linkTo()



###########################################################################/**
# @RdocMethod renameTo
#
# @title "Renames/moves a file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{filename, path}{The filename and the path for the destination file.
#   The default is to use the same filename as the source file.
#   The destination pathname must not be the same as the source file,
#   otherwise an exception is thrown.}
#  \item{...}{Additional arguments passed to @see "R.utils::renameFile".}
#  \item{verbose}{...}
# }
#
# \value{
#   Returns the soure @see "GenericDataFile".
# }
#
# @author
#
# \seealso{
#   Internally @see "R.utils::renameFile".
#   @seemethod "copyTo".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("renameTo", "GenericDataFile", function(this, filename=getFilename(this), path=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'filename' and 'path':
  pathname <- Arguments$getWritablePathname(filename, path=path);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Assert that the file exists
  if (!isFile(this)) {
    throw("Cannot rename file. Source file does not exist: NA");
  }

  # Nothing to do?
  srcPathname <- getPathname(this, NULLasNA=TRUE)
  if (identical(pathname, srcPathname))
    return(this);

  # Assert that file is not overwritten by mistake.
  pathname <- Arguments$getWritablePathname(pathname, mustNotExist=TRUE);

  verbose && enter(verbose, "Renaming ", class(this)[1], " pathname");
  verbose && cat(verbose, "Source: ", srcPathname);
  verbose && cat(verbose, "Destination: ", pathname);

  verbose && enter(verbose, "Renaming file");
  ## FIXME: Should we use R.utils::renameFile() instead? /HB 2014-01-06
  res <- renameFile(srcPathname, pathname, ...);
  if (!res) {
    throw("Failed to rename file: ", srcPathname, " -> ", pathname);
  }
  verbose && exit(verbose);

  # Update GenericDataFile object
  this$.pathname <- pathname;

  verbose && exit(verbose);

  invisible(this);
}, protected=TRUE)



###########################################################################/**
# @RdocMethod getChecksum
#
# @title "Gets the checksum of a file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{write}{If @TRUE or @NA and a checksum file does not exists, then
#    a checksum file is created, iff possible.  If @NA and the file could
#    not be created, then it falls back to @FALSE, but if @TRUE an error
#    is thrown.  If @FALSE and no checksum file exists, the checksum is
#    calculated on the fly.}
#  \item{force}{If @FALSE, the file exists and has not be modified since,
#    then the cached checksum is returned.}
#  \item{verbose}{...}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character string, which can be @NA if file is missing.
# }
#
# @author
#
# \seealso{
#   Internally @see "digest::digest" is used.
#   @seemethod "readChecksum".
#   @seemethod "writeChecksum".
#   @seemethod "compareChecksum".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getChecksum", "GenericDataFile", function(this, write=NA, force=FALSE, verbose=FALSE, ...) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  checksum <- this$.checksum;
  if (force || is.null(checksum) || hasBeenModified(this, update=FALSE)) {
    if (isFile(this)) {
      dfZ <- NULL
      if (is.na(write)) {
        dfZ <- tryCatch({
          getChecksumFile(this, force=force)
        }, error=function(ex) { NULL })
      } else if (write) {
        dfZ <- getChecksumFile(this, force=force)
      } else if (hasChecksumFile(this)) {
        dfZ <- getChecksumFile(this, force=force)
      }

      if (isFile(dfZ)) {
        checksum <- readChecksum(dfZ)
      } else {
        verbose && enter(verbose, "Calculating checksum");
        pathname <- getPathname(this);
        checksum <- digest(pathname, file=TRUE);
        verbose && exit(verbose);
      }
    } else {
      checksum <- NA_character_
    }

    ## Update checksum
    this$.checksum <- checksum
  }

  checksum;
})


###########################################################################/**
# @RdocMethod writeChecksum
#
# @title "Write the file checksum to a checksum file"
#
# \description{
#   @get "title" having the same filename with suffix \code{.md5} added.
# }
#
# @synopsis
#
# \arguments{
#  \item{skip}{If @TRUE, an already written checksum file is skipped.}
#  \item{...}{Not used.}
#  \item{verbose}{...}
# }
#
# \value{
#   Returns (invisibly) the pathname to the checksum file.
#   An exception is thrown if the file does not exist.
# }
#
# @author
#
# \seealso{
#   @seemethod "validateChecksum".
#   @seemethod "compareChecksum".
#   @seemethod "readChecksum".
#   @seemethod "getChecksum".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("writeChecksum", "GenericDataFile", function(this, ..., skip=FALSE) {
  dfZ <- getChecksumFile(this, force=!skip, ...)
  invisible(getPathname(dfZ))
})



###########################################################################/**
# @RdocMethod readChecksum
#
# @title "Reads the value of the corresponding checksum file"
#
# \description{
#   @get "title", if existing.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character string, which can be @NA if file is missing.
#   An exception is thrown if the file does not exist, and hence not
#   the checkum file.
# }
#
# \details{
#   The content of the checksum file is trimmed from comment lines,
#   whitespaces and then validated that the remaining part contains a
#   hexadecimal value.
# }
#
# @author
#
# \seealso{
#   @seemethod "validateChecksum".
#   @seemethod "compareChecksum".
#   @seemethod "writeChecksum".
#   @seemethod "getChecksum".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("readChecksum", "GenericDataFile", function(this, ...) {
  if (!isFile(this)) {
    throw("Cannot read stored checksum. File does not exist: NA");
  }

  if (!hasChecksumFile(this)) {
    throw("Cannot read stored checksum. No checksum file available: ", getPathname(this))
  }

  dfZ <- getChecksumFile(this)
  readChecksum(dfZ)
}, protected=TRUE)


###########################################################################/**
# @RdocMethod compareChecksum
#
# @title "Compares the file checksum with the value of the checksum file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
#  \item{verbose}{...}
# }
#
# \value{
#   Returns @TRUE if the file checksum is identical to the stored value
#   in the corresponding checksum file, otherwise @FALSE.  @FALSE is
#   also returned if the checksum file does not exist.
#   An exception is thrown if the file does not exist.
# }
#
# @author
#
# \seealso{
#   @seemethod "validateChecksum".
#   @seemethod "readChecksum".
#   @seemethod "writeChecksum".
#   @seemethod "getChecksum".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("compareChecksum", "GenericDataFile", function(this, ...) {
  if (!isFile(this)) {
    throw("Cannot compare checksum. File does not exist: NA");
  } else if (!hasChecksumFile(this)) {
    return(FALSE)
  }

  tryCatch({
    validateChecksum(this, ...)
    TRUE
  }, error=function(ex) FALSE)
})


###########################################################################/**
# @RdocMethod validateChecksum
#
# @title "Asserts that the file checksum matches the one of the checksum file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   An exception is thrown if not, that is, if @seemethod "compareChecksum"
#   returns @FALSE.
# }
#
# @author
#
# \seealso{
#   @seemethod "validateChecksum".
#   @seemethod "readChecksum".
#   @seemethod "writeChecksum".
#   @seemethod "getChecksum".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("validateChecksum", "GenericDataFile", function(this, ...) {
  if (!isFile(this)) {
    throw("Cannot validate checksum. File does not exist: NA");
  } else if (!hasChecksumFile(this)) {
    throw("Cannot validate checksum. No checksum file available: ", getPathname(this))
  }

  dfZ <- getChecksumFile(this)
  res <- validate(dfZ)
  invisible(res)
})





###########################################################################/**
# @RdocMethod gzip
# @aliasmethod gunzip
# @aliasmethod isGzipped
# @alias gzip
# @alias gunzip
# @alias isGzipped
#
# @title "Compresses/uncompresses a file"
#
# \description{
#   @get "title" using gzip compression.
#   When compressing (uncompressing), the new filename has suffix \code{.gz}
#   appended (removed), which is also used to test if a file is gzip'ed
#   or not.
# }
#
# \usage{
#  @usage gzip,GenericDataFile
#  @usage gunzip,GenericDataFile
#  @usage isGzipped,GenericDataFile
# }
#
# \arguments{
#  \item{...}{Arguments passed to "R.utils::gzip" and "R.utils::gunzip",
#    respectively.}
#  \item{verbose}{...}
# }
#
# \value{
#   Returns (invisibly) the updated pathname.
#   When compressing (uncompressing), an exception is thrown if the file
#   is already compressed (not compressed).
#   An exception is thrown if the file does not exist.
# }
#
# @author
#
# \seealso{
#   Internally @see "R.utils::gzip" is used.
#   @seeclass
# }
#*/###########################################################################
setMethodS3("gzip", "GenericDataFile", function(this, ...) {
  if (!isFile(this)) {
    throw("Cannot gzip file. File does not exist: NA");
  }

  pathname <- getPathname(this);
  if (isGzipped(pathname)) {
    throw("File is already gzip'ed: ", pathname);
  }

  outPathname <- gzip(pathname, ...);

  # Attributes that may be of interest
  ## temporary <- attr(outPathname, "temporary");
  ## nbrOfBytes <- attr(outPathname, "nbrOfBytes");

  this$.pathname <- outPathname;

  # Return the output file
  invisible(outPathname);
}, protected=TRUE)



setMethodS3("gunzip", "GenericDataFile", function(this, ...) {
 # Argument 'this':
  if (!isFile(this)) {
    throw("Cannot gunzip file. File does not exist: NA");
  }
  pathname <- getPathname(this);
  if (!isGzipped(pathname)) {
    throw("File is not gzip'ed: ", pathname);
  }

  # Decompress
  outPathname <- gunzip(pathname, ...);

  # Attributes that may be of interest
  ## temporary <- attr(outPathname, "temporary");
  ## nbrOfBytes <- attr(outPathname, "nbrOfBytes");

  this$.pathname <- outPathname;

  invisible(outPathname);
}, protected=TRUE)


setMethodS3("isGzipped", "GenericDataFile", function(this, ...) {
  pathname <- getPathname(this, NULLasNA=TRUE)
  isGzipped(pathname, ...);
}, protected=TRUE)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ODDS AND ENDS
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("renameToUpperCaseExt", "GenericDataFile", function(static, pathname, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # isFileCase() is a case-sensitive isFile() for Windows
  isFileCase <- function(pathname, ...) {
    # Non-case sensitive check
    if (!isFile(pathname))
      return(FALSE);

    # There can still be a case-difference
    path <- dirname(pathname);
    filename <- basename(pathname);
    filenames <- list.files(path=path, all.files=TRUE);
    res <- grep(filename, filenames, fixed=TRUE);
    res <- (length(res) >= 1L);
    res;
  }


  # Identify the filename extension
  ext <- gsub("^.*[.]", "", pathname);

  # No filename extension?  Do nothing.
  if (identical(ext, pathname))
    return(pathname);

  # Generate the pathname with lower- and upper-case extensions
  extL <- tolower(ext);
  pathnameL <- gsub(sprintf("[.]%s$", ext), sprintf(".%s", extL), pathname);

  extU <- toupper(ext);
  pathnameU <- gsub(sprintf("[.]%s$", ext), sprintf(".%s", extU), pathname);

  # Does a lower-case filename exist? If not, nothing to do.
  if (!isFileCase(pathnameL))
    return(pathnameU);

  # Can we rename it?
  if (identical(pathname, pathnameL) && isFileCase(pathnameU)) {
    throw("Cannot rename pathname to have upper-case filename extension, because such a file already exists: ", pathnameU);
  }

  # Try to rename the file
  res <- file.rename(pathnameL, pathnameU);
  if (res) {
    msg <- paste("Renamed file to have an upper-case filename extension:", pathname);
    warning(msg);
  } else {
    throw("Failed to rename file such that it gets an upper-case filename extension (try to rename the file manually): ", pathname);
  }

  pathnameU;
}, static=TRUE, protected=TRUE)



############################################################################
# HISTORY:
# 2015-12-23
# o Now is.na() for GenericDataFile returns TRUE also if pathname is NULL.
# 2015-12-13
# o Now argument 'filename' must be a single string - not NULL.
# 2014-08-26
# o Now gzip()/gunzip() for GenericDataFile returns the output pathname.
# 2014-02-28
# o DOCUMENTATION: Added an Rd section on privileges required on Windows
#   for linkTo() to work.
# 2014-01-07
# o Added is.na() for GenericDataFile.
# 2014-01-06
# o Added linkTo() for GenericDataFile.
# o CLEANUP: Made copyTo() a light-weight wrapper around copyFile(),
#   which takes care of all the validation.
# o copyTo() for GenericDataFile passes '...' to R.utils::copyFile().
# 2014-01-05
# o CLEANUP: copyTo() and renameTo() flr GenericDataFile had verbose
#   output enabled by default.
# 2013-06-20
# o Now isGzipped() for GenericDataFile ignores the case of the filename
#   extension when testing whether the file is gzipped or not.
# 2013-05-22
# o CLEANUP: Now getChecksum() for GenericDataFile calls
#   digest::digest(..., file=TRUE) instead of digest2(..., file=TRUE),
#   which would give identical results.
# 2012-12-06
# o Added argument 'absolute=FALSE' to getPathname() for GenericDataFile.
# 2012-12-03
# o Added isGzipped() to GenericDataFile.
# 2012-11-28
# o Now GenericDataFile() retrieves the file time stamps such that
#   hasBeenModified() returns a correct value also when first called,
#   and not only TRUE just in case.
# o Now declaring the .checksum field as a 'cached' field, which will
#   be cleared if clearCache() or gc() is called on the object.
# 2012-10-30
# o Added validate() to GenericDataFile.
# 2011-11-19
# o Now more methods for GenericDataFile handle a so called "empty" file,
#   which is a file with pathname equal to NULL.
# 2011-04-04
# o BUG FIX: equals() for GenericDataFile would consider to files not
#   to be equal only if their checksums was equal, and vice versa.
#   Also, when creating the message string explaining why they differ
#   an error would have been thrown.
# 2010-02-13
# o Added argument '.onUnknownArgs' to GenericDataFile().
# 2010-01-31
# o DOCUMENTATION: Added Rdoc comments to most methods.
# o CLEAN UP: Made readChecksum() of GenericDataFile protected.
# 2010-01-04
# o Now setExtensionPattern(..., pattern=NULL) of GenericDataFile works.
# o Added argument 'default="\\.([^.]+)$"' to getExtensionPattern() of
#   GenericDataFile.
# 2009-12-30
# o BUG FIX: Now GenericDataFile(pathname) throws an error if 'pathname'
#   is refering to a directory.
# o Now GenericDataFile(NA, mustExist=FALSE) is a valid object.  Made all
#   methods aware of such missing files.
# 2009-10-02
# o CLEAN UP: Removed setFullName() for GenericDataFile, because there
#   is not a "default" on.
# o Now setFullNameTranslator(...) for GenericDataFile dispatches on the
#   'translator' argument (2nd) to call setFullNameTranslatorBy<class>().
#   setFullNameTranslatorByFunction() and setFullNameTranslatorByNULL()
#   are defined by default.
# 2009-05-19
# o Now setFullNameTranslator() for GenericDataFile asserts that the
#   fullname translator function accepts also argument 'file'.
# 2009-04-23
# o BUG FIX: as.character() would throw 'Error in sprintf("%d", getFileSize
#   (db, "numeric")) : use format %f, %e, %g or %a for numeric objects' if
#   file size is returned as a double, which happens for very large files.
# 2009-02-26
# o Now hasTags(..., tags) splits the 'tags' argument.
# 2009-02-23
# o (Re-)Added space between number and unit for RAM and file size.
# o Now as.character() of GenericDataFile also reports the exact file size
#   in case the file size is reported in kB, MB, etc.  It also tries to
#   report the relative pathname rather than the absolute.
# o Now getChecksum() of GenericDataFile caches results unless the file
#   has been modified since last time.
# o Added hasBeenModified() to GenericDataFile.
# 2008-09-18
# o Now readChecksum() does some validation.  It is also possible to have
#   commented rows in the checksum file.
# o Added argument 'skip' to writeChecksum().
# 2008-07-23
# o Added getCreatedOn(), getLastModifiedOn(), getLastAccessedOn().
# 2008-07-17
# o Added setFullName(), and setName().
# 2008-06-11
# o BUG FIX: Argument 'recursive' is ignored in the GenericDataFile
#   constructor.
# 2008-06-06
# o Added argument 'recursive' to fromFile().
# 2008-05-16
# o Added support for full name translation of GenericDataFile:s.
# 2008-05-15
# o Added gzip() and gunzip().
# o Update equals() to also compare classes, file sizes, and checksums.
# 2008-05-11
# o Now static fromFile() always creates an instance.
# 2008-05-09
# o Moved private get/seLabel() from AffymetrixFile to GenericDataFile.
# o Moved the attributes features from AffymetrixFile to GenericDataFile.
# 2008-03-22
# o Added 'aliased' to getFullName().
# 2007-09-25
# o Added isFile() to test if the file exists or not.
# 2007-09-15
# o Added renameTo().
# o Now copyTo() utilizes fileCopy(), which is fail safe.
# 2007-09-14
# o Extracted GenericDataFile from AffymetrixFile.
# 2007-09-13
# o Added missing setAttributesByTags().
# 2007-08-09
# o Added static renameToUpperCaseExt().
# 2007-03-20
# o Added getAlias() and setAlias().  Note, getName() etc are still
#   unaffected by these.
# 2007-03-05
# o Added setAttributesByTags(), which now also tries to coerce values.
# o Added support for (in-memory) attributes.
# 2007-02-07
# o Added getChecksum(), writeChecksum(), readChecksum(), and
#   compareChecksum() and validateChecksum(). I did this because I noticed
#   by chance that some of my CEL files transferred via an external HDD got
#   corrupt probe signals.
# 2007-01-14
# o Added a test for "unknown" (=unused) arguments to constructor.
# 2007-01-07
# o Added hasTags() and hasTag().
# 2006-11-02
# o Added getFullName(), getTags() and redefined getName().
# 2006-09-15
# o Added stextSize().
# 2006-08-27
# o Added stextLabel() and stextLabels(). stext is for "side text", cf.
#   mtext for "margin text". stext() is slightly more convenient than mtext
#   when it comes to different font sizes.
# o Added copyTo().
# 2006-08-14
# o Added abstract fromFile().
# 2006-08-11
# o Created from AffymetrixDataFile in order to represent CDF files too.
############################################################################
