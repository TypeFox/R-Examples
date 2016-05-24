###########################################################################/**
# @RdocClass AbstractFileArray
#
# @title "Class representing a persistent array stored in a file"
#
# \description{
#  @classhierarchy
#
#  Note that this is an abstract class, i.e. it is not possible to create
#  an object of this class but only from one of its subclasses.
#  For a vector data type, see @see "FileVector".
#  For a matrix data type, see @see "FileMatrix".

# }
#
# @synopsis
#
# \arguments{
#  \item{filename}{The name of the file storing the data.}
#  \item{path}{An optional path where data should be stored.}
#  \item{storageMode}{The storage @see "mode" of the data elements.}
#  \item{bytesPerCell}{The number of bytes each element (cell) takes up
#      on the file system.  If \code{NULL}, it is inferred from the
#      \code{storageMode} argument.}
#  \item{dim}{A @numeric @vector specifying the dimensions of the array.}
#  \item{dimnames}{An optional @list of dimension names.}
#  \item{dimOrder}{The order of the dimensions.}
#  \item{comments}{An optional @character string of arbitrary length.}
#  \item{nbrOfFreeBytes}{The number of "spare" bytes after the comments
#      before the data section begins.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
#
# }
#
# \details{
#   The purpose of this class is to be able to work with large arrays
#   in \R without being limited by the amount of memory available.
#   Data is kept on the file system and elements are read and written
#   whenever queried.
# }
#
# \section{Maximum number of elements}{
#   It is only the header that is kept in memory, not the data, and
#   therefore the maximum length of a array that can be allocate, is limited
#   by the amount of available space on the file system.
#   Since element names (optional) are stored in the header,
#   these may also be a limiting factor.
# }
#
# \section{Element names}{
#   The element names are stored in the header and are currently read and
#   written to file one by one.  This may slow down the performance
#   substantially if the dimensions are large.  For optimal opening
#   performance, avoid names.
#
#   For now, do \emph{not} change names after file has been allocated.
# }
#
# \section{File format}{
#   The file format consist of a header section and a data section.
#   The header contains information about the file format, the length
#   and element names of the array, as well as data type
#   (storage @see "mode"), the size of each element.
#   The data section, which follows immediately after the header section,
#   consists of all data elements with non-assigned elements being
#   pre-allocated with zeros.
#
#   For more details, see the source code.
# }
#
# \section{Limitations}{
#   The size of the array in bytes is limited by the maximum file size
#   of the file system.
#   For instance, the maximum file size on a Windows FAT32 system is
#   4GB (2GB?).  On Windows NTFS the limit is in practice ~16TB.
# }
#
# @author
#
# \references{
#  [1] New Technology File System (NTFS), Wikipedia, 2006
#      \url{http://en.wikipedia.org/wiki/NTFS}.
# }
#
# @visibility public
#*/###########################################################################
setConstructorS3("AbstractFileArray", function(filename=NULL, path=NULL, storageMode=c("integer", "double"), bytesPerCell=1, dim=NULL, dimnames=NULL, dimOrder=NULL, comments=NULL, nbrOfFreeBytes=4096) {
##   # ROBUSTNESS/WORKAROUND: For now, package attaches the 'R.oo' package.
##   # This is needed due to what appears to be a bug in how R.oo
##   # finalizes Object:s assuming R.oo is/can be attached.  Until that
##   # is resolved, we make sure R.oo is attached. /HB 2013-09-21
##   pkg <- "R.oo";
##   suppressPackageStartupMessages(require(pkg, character.only=TRUE, quietly=TRUE)) || throw("Package not loaded: ", pkg);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'filename' and 'path':
  if (!is.null(filename)) {
    pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=FALSE);
  } else {
    pathname <- NULL;
  }

  # Argument 'storageMode':
  storageMode <- match.arg(storageMode);

  # Argument 'nbrOfFreeBytes':
  nbrOfFreeBytes <- Arguments$getNumeric(nbrOfFreeBytes, range=c(0,Inf));

  # Argument 'bytesPerCell':
  bytesPerCell <- Arguments$getInteger(bytesPerCell, range=c(1,8));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create header
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  header <- NULL;
  if (!is.null(pathname)) {
    if (!isFile(pathname)) {
      if (is.null(dimOrder))
        dimOrder <- seq(along=dim);

      # Create header
      header <- list(
        magic = "",
        version = "0.1.1",
        .commentsOffset = NULL,
        .dataOffset = NULL,
        storageMode = as.character(storageMode),
        bytesPerCell = as.integer(bytesPerCell),
        dim = dim,
        dimOrder = dimOrder,
        dimnames = dimnames,
        comments = as.character(comments),
        nbrOfFreeBytes = as.double(nbrOfFreeBytes)
      );
    }
  }

  this <- extend(Object(), "AbstractFileArray",
    .pathname = as.character(pathname),
    con = NULL,
    header = header
  )

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(pathname)) {
    open(this);
  }

  this;
}, abstract=TRUE)



###########################################################################/**
# @RdocMethod as.character
#
# @title "Returns a short string describing the file array"
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
setMethodS3("as.character", "AbstractFileArray", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- paste(class(this)[1], ":", sep="");
  s <- paste(s, " Pathname: ", getPathname(this), ".", sep="");
  s <- paste(s, " Opened?: ", isOpen(this), ".", sep="");
  fileSize <- getFileSize(this);
  s <- paste(s, " File size: ", fileSize, " bytes", sep="");
  if (fileSize > 1024^3) {
    s <- paste(s, sprintf(" (%.3g Gb)", fileSize/1024^3), sep="");
  } else if (fileSize > 1024^2) {
    s <- paste(s, sprintf(" (%.3g Mb)", fileSize/1024^2), sep="");
  } else if (fileSize > 1024) {
    s <- paste(s, sprintf(" (%.3g kb)", fileSize/1024), sep="");
  }
  s <- paste(s, ".", sep="");
  s <- paste(s, " Dimensions: ", paste(dim(this), collapse="x"), ".", sep="");
  s <- paste(s, " Number of elements: ", length(this), ".", sep="");
  s <- paste(s, " Bytes per cell: ", getBytesPerCell(this), ".", sep="");
  s;
})




###########################################################################/**
# @RdocMethod isOpen
#
# @title "Checks whether the data file of the file array is open or not"
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
#  Returns nothing.
# }
#
# @author
#
# \seealso{
#   @seemethod "open".
#   @seemethod "close".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("isOpen", "AbstractFileArray", function(this, ...) {
  if (is.null(this$con))
    return(FALSE);

  # Why the tryCatch()? /HB 2013-09-21
  res <- FALSE;
  tryCatch({
    res <- isOpen(this$con);
##    res <- base::isOpen(this$con);
  }, error = function(ex) {
    print(ex)
    this$con <- NULL;
  })

  res;
}, createGeneric=FALSE)



###########################################################################/**
# @RdocMethod open
#
# @title "Opens a connection to the data file of the file array"
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
#  Returns (invisibly) itself.
# }
#
# @author
#
# \seealso{
#   @seemethod "isOpen".
#   @seemethod "close".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("open", "AbstractFileArray", function(con, ...) {
  # To please R CMD check
  this <- con;

  if (isOpen(this))
    throw("File already open.");

  pathname <- getPathname(this);
  fileExists <- isFile(pathname);
  if (fileExists) {
    con <- NULL;
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Try 1: Open the file for read and updating
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tryCatch({
      pathname <- Arguments$getWritablePathname(pathname);
      # Open existing file: read header.
      # Note, if open="w+b", the file is trucated to size zero.
      con <- file(pathname, open="r+b");
    }, error = function(ex) {});

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Try 2: Open the file for reading only ("better than nothing")
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (is.null(con)) {
      pathname <- Arguments$getReadablePathname(pathname);
      con <- file(pathname, open="rb");
    }
    this$con <- con;
    if (is.null(this$header)) {
      this$header <- readHeader(this);
    }
  } else {
    pathname <- Arguments$getWritablePathname(pathname);
    # Create a new file: open conection, create header and data section.
    this$con <- file(pathname, open="w+b");
    this <- writeHeader(this);
    writeEmptyData(this);
    flush(this);
  }

  invisible(this);
})



###########################################################################/**
# @RdocMethod close
#
# @title "Closes a connection to the data file of the file array"
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
#  Returns @TRUE if the file was closed.  If the file is not opened,
#  an exception is thrown.
# }
#
# @author
#
# \seealso{
#   @seemethod "isOpen".
#   @seemethod "open".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("close", "AbstractFileArray", function(con, ...) {
  # To please R CMD check
  this <- con;

  if (!isOpen(this)) {
    this$con <- NULL;
    throw("File not opened.");
  }

  con <- this$con;
  if (!is.null(con) && isOpen(con)) {
##  if (!is.null(con) && base::isOpen(con)) {
    flush(con);
    close(con);
    con <- NULL;
  }

  this$con <- NULL;
  invisible(TRUE);
})


###########################################################################/**
# @RdocMethod finalize
#
# @title "Internal: Clean up when file array is deallocated from memory"
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
#  Returns nothing.
# }
#
# \details{
#   The finalizer of a file array makes sure to close the file connection,
#   if it is open.
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
setMethodS3("finalize", "AbstractFileArray", function(this, ...) {
  if (isOpen(this))
    close(this);
}, protected=TRUE)



###########################################################################/**
# @RdocMethod delete
#
# @title "Deletes the file array from the file system"
#
# \description{
#  @get "title".
#  If the file array is open, it is first closed.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns (invisibly) @TRUE if the file was successfully deleted (or did
#   not exist in the first place), otherwise @FALSE.
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
setMethodS3("delete", "AbstractFileArray", function(this, ...) {
  if (isOpen(this))
    close(this);

  if (!isFile(getPathname(this)))
    return(invisible(TRUE));

  # Delete the actual file
  pathname <- getPathname(this);
  pathname <- Arguments$getWritablePathname(pathname);
  invisible(file.remove(pathname));
})




###########################################################################/**
# @RdocMethod flush
#
# @title "Internal: Flushes the write buffer"
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
#   Returns nothing.
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
setMethodS3("flush", "AbstractFileArray", function(con, ...) {
  # To please R CMD check
  this <- con;

  if (isOpen(this)) {
    flush(this$con);
  }
}, protected=TRUE)




###########################################################################/**
# @RdocMethod clone
#
# @title "Clones a file array"
#
# \description{
#  @get "title" including the file on the file system.
# }
#
# @synopsis
#
# \arguments{
#   \item{copyData}{If @TRUE, also the data file is copied, otherwise not.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns the new @see "AbstractFileArray" object.
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
setMethodS3("clone", "AbstractFileArray", function(con, copyData=TRUE, ...) {
  # To please R CMD check
  this <- con;

  clone <- NextMethod("clone");
  clone$con <- NULL;

  # Generate new filename
  path <- getPath(this);
  name <- getName(this);
  cloneId <- getCloneNumber(this);
  if (is.na(cloneId))
    cloneId <- 0;
  ext <- getExtension(this);

  MAX.COUNT <- 256;
  count <- 0;
  ready <- FALSE;
  while(!ready && count < MAX.COUNT) {
    cloneId <- cloneId + 1;
    basename <- sprintf("%s.%04d.%s", name, as.integer(cloneId), ext);
    pathname <- file.path(path, basename);
    if (!isFile(pathname)) {
      # Copy file?
      if (copyData) {
        pathname <- Arguments$getWritablePathname(pathname);
        if (file.copy(getPathname(this), pathname))
          ready <- TRUE;
      } else {
        # ...otherwise, the structure will be create when the opened!
        ready <- TRUE;
      }
    }
    count <- count + 1;
  }

  if (!ready)
    throw("Failed to copy file: ", getPathname(this));

  clone$.pathname <- pathname;

  # Open connection?
  if (isOpen(this)) {
    open(clone);
  }

  clone;
})



###########################################################################/**
# @RdocMethod getPathname
#
# @title "Gets the full pathname to the data file"
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
#   Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seemethod "getPath".
#   @seemethod "getBasename".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("getPathname", "AbstractFileArray", function(this, ...) {
  file.path(getPath(this), getBasename(this));
})


###########################################################################/**
# @RdocMethod getPath
#
# @title "Gets the path (directory) where the data file lives"
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
#   Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seemethod "getPathname".
#   @seemethod "getBasename".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("getPath", "AbstractFileArray", function(this, ...) {
  dirname(this$.pathname);
})


###########################################################################/**
# @RdocMethod getBasename
#
# @title "Gets the basename (filename) of the data file"
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
#   Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seemethod "getPathname".
#   @seemethod "getBasename".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("getBasename", "AbstractFileArray", function(this, ...) {
  basename(this$.pathname);
})


###########################################################################/**
# @RdocMethod getName
#
# @title "Gets the name of the file array"
#
# \description{
#  @get "title" from the filename of the data file.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seemethod "getBasename".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("getName", "AbstractFileArray", function(this, ...) {
  name <- getBasename(this);

  # Strip extension.
  name <- gsub("[.][^.]*$", "", name);

  pattern <- "^(.*)[.]([^.]*)$";
  if (regexpr(pattern, name) == -1)
    return(name);

  res <- gsub(pattern, "\\1", name);

  res;
})



###########################################################################/**
# @RdocMethod getExtension
#
# @title "Gets the filename extension of the file array"
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
#   Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seemethod "getName".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("getExtension", "AbstractFileArray", function(this, ...) {
  name <- getBasename(this);

  pattern <- "^(.*)[.]([^.]*)$";
  if (regexpr(pattern, name) == -1)
    return(NA);

  res <- gsub(pattern, "\\2", name);

  res;
})


###########################################################################/**
# @RdocMethod getFileSize
#
# @title "Gets the size of the file array"
#
# \description{
#  @get "title" (in bytes).
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns the number of bytes as a @number.
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
setMethodS3("getFileSize", "AbstractFileArray", function(this, ...) {
  file.info(getPathname(this))$size;
})



###########################################################################/**
# @RdocMethod getCloneNumber
#
# @title "Gets the clone number of the file array"
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
#   Returns an @integer.
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
setMethodS3("getCloneNumber", "AbstractFileArray", function(this, ...) {
  name <- getBasename(this);

  # Strip extension.
  name <- gsub("[.][^.]*$", "", name);

  pattern <- "^(.*)[.]([^.]*)$";
  if (regexpr(pattern, name) == -1)
    return(as.integer(NA));

  res <- gsub(pattern, "\\2", name);

  as.integer(res);
}, protected=TRUE)



###########################################################################/**
# @RdocMethod getSizeOfComments
#
# @title "Gets the number of bytes the comments occupies"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{comments}{The comments to be measured.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns an @integer.
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
setMethodS3("getSizeOfComments", "AbstractFileArray", function(this, comments=this$header$comments, ...) {
  count <- 0;

  str <- comments;
  if (length(str) == 0) {
    len <- 0;
  } else {
    len <- nchar(str);
  }

  count <- 4;

  if (len > 0) {
    count <- count + len + 1;
  }

  count;
}, protected=TRUE);



###########################################################################/**
# @RdocMethod getComments
#
# @title "Gets the comments for this file array"
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
#   Returns a @character @vector.
# }
#
# @author
#
# \seealso{
#   @seemethod "setComments".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("getComments", "AbstractFileArray", function(this, ...) {
  this$header$comments;
})



###########################################################################/**
# @RdocMethod setComments
#
# @title "Sets the comments for this file array"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{comments}{A @character @vector of new comments.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns (invisibly) the old comments.
# }
#
# @author
#
# \seealso{
#   @seemethod "getComments".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("setComments", "AbstractFileArray", function(this, comments=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'comments':
  comments <- Arguments$getCharacter(comments);


  header <- this$header;
  oldComments <- header$comments;

  if (length(comments) == 0)
    return(invisible(oldComments));

  # Assert that the new comments can be added
  offset <- header$.commentsOffset;
  nbrOfFreeBytes <- header$.dataOffset - (offset + getSizeOfComments(this, comments=comments));
  if (nbrOfFreeBytes < 0) {
    throw("Not enough of space to set new comments because that would overwrite data by ", nbrOfFreeBytes, " byte(s): ", comments);
  }

  # Update the header in memory
  this$header$comments <- comments;
  this$header$nbrOfFreeBytes <- nbrOfFreeBytes;

  # Finally, update the file
  writeHeaderComments(this);

  invisible(oldComments);
})


setMethodS3("writeHeaderComments", "AbstractFileArray", function(this, ...) {
  writeString <- function(con, str) {
    count <- 0;

    if (length(str) == 0) {
      len <- 0;
    } else {
      len <- nchar(str);
    }
    writeBin(con=con, as.integer(len), size=4);
    count <- 4;

    if (len > 0) {
      writeChar(con=con, str);
      count <- count + len + 1;
    }

    count;
  } # writeString()

  header <- this$header;
  con <- this$con;
  offset <- header$.commentsOffset;

  nbrOfFreeBytes <- header$.dataOffset - (offset + getSizeOfComments(this));
  if (nbrOfFreeBytes < 0) {
    throw("Internal error: Cannot write header comments because that would overwrite data by ", nbrOfFreeBytes, " byte(s).");
  }

  # Comments
  seek(con, where=offset, rw="write");
  invisible(writeString(con, header$comments));
}, protected=TRUE)



###########################################################################/**
# @RdocMethod writeHeader
#
# @title "Writes the header of a file array to file"
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
#  Returns nothing.
# }
#
# @author
#
# \seealso{
#   @seemethod "writeHeader".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("writeHeader", "AbstractFileArray", function(this, ...) {
  writeString <- function(con, str) {
    count <- 0;

    if (length(str) == 0) {
      len <- 0;
    } else {
      len <- nchar(str);
    }
    writeBin(con=con, as.integer(len), size=4);
    count <- 4;

    if (len > 0) {
      writeChar(con=con, str);
      count <- count + len + 1;
    }

    count;
  } # writeString()

  header <- this$header;
  con <- this$con;
  offset <- 0;

  # MAGIC
  seek(con, where=offset, rw="write");
  writeChar(con=con, class(this)[1]);
  offset <- offset + 64;

  # File format version string
  seek(con, where=offset, rw="write");
  writeChar(con=con, header$version);
  offset <- offset + 16;

  # Offset to data section (double = 8 bytes)
  # To be written later
  offsetDataOffset <- offset;
  offset <- offset + 8;

  # Offset to comments (double = 8 bytes)
  # To be written later
  offsetCommentsOffset <- offset;
  offset <- offset + 8;

  # Number of dimensions (double = 8 bytes)
  nbrOfDims <- length(header$dim);
  seek(con, where=offset, rw="write");
  writeBin(con=con, as.double(nbrOfDims), size=8);
  offset <- offset + 8;

  # Length of each dimension
  for (kk in seq(length=nbrOfDims)) {
    # Length of dimension 'kk' (double = 8 bytes)
    seek(con, where=offset, rw="write");
    writeBin(con=con, as.double(header$dim[kk]), size=8);
    offset <- offset + 8;
  }

  # Order of dimensions
  seek(con, where=offset, rw="write");
  writeBin(con=con, as.double(header$dimOrder), size=8);
  offset <- offset + 8 * nbrOfDims;

  # Bytes per cell
  seek(con, where=offset, rw="write");
  writeBin(con=con, as.integer(header$bytesPerCell), size=4);
  offset <- offset + 4;

  # R storage mode
  seek(con, where=offset, rw="write");
  offset <- offset + writeString(con, header$storageMode);

  # Dimension names
  for (kk in seq(length=nbrOfDims)) {
    # Number of names for dimension 'kk' (double = 8 bytes)
    names <- header$dimnames[[kk]];
    nbrOfNames <- length(names);
    seek(con, where=offset, rw="write");
    writeBin(con=con, as.double(nbrOfNames), size=8);
    offset <- offset + 8;

    for (ll in seq(length=nbrOfNames)) {
      offset <- offset + writeString(con, names[ll]);
    }
    names <- NULL; # Not needed anymore
  }

  # Update 'commentOffset' field
  this$header$.commentsOffset <- as.double(offset);

  # Comments
  offset <- offset + writeString(con, header$comments);

  # Free space (random bytes)
  seek(con, where=offset, rw="write");
  offset <- offset + header$nbrOfFreeBytes;

  # Update 'dataOffset' field
  this$header$.dataOffset <- as.double(offset);

  # Offset to data section
  seek(con, where=offsetCommentsOffset, rw="write");
  writeBin(con=con, as.double(this$header$.commentsOffset), size=8);

  # Offset to data section
  seek(con, where=offsetDataOffset, rw="write");
  writeBin(con=con, as.double(this$header$.dataOffset), size=8);

  invisible(this);
}, protected=TRUE) # writeHeader()



###########################################################################/**
# @RdocMethod readHeader
#
# @title "Read the header of a file array data file"
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
#  Returns @list structure.
# }
#
# @author
#
# \seealso{
#   @seemethod "writeHeader".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("readHeader", "AbstractFileArray", function(this, ...) {
  readString <- function(con) {
    len <- readBin(con=con, what=integer(0), size=4);
    count <- 4;

    if (len > 0) {
      str <- readChar(con=con, nchars=len);
      count <- count + len + 1;
    }

    attr(str, "count") <- as.integer(count);
    str;
  } # readString()

  con <- this$con;
  offset <- 0;

  header <- list();

  # MAGIC
  seek(con, where=offset, rw="read");
  value <- readChar(con=con, nchars=64);
  value <- substring(value, 1);
  header$magic <- value;
  offset <- offset + 64;

  # File format version string
  seek(con, where=offset, rw="read");
  value <- readChar(con=con, nchars=16);
  value <- substring(value, 1);
  header$version <- value;
  offset <- offset + 16;

  # Offset to data section (double = 8 bytes)
  seek(con, where=offset, rw="read");
  header$.dataOffset <- readBin(con=con, what=double(0), size=8);
  offset <- offset + 8;

  # Offset to comment section (double = 8 bytes)
  seek(con, where=offset, rw="read");
  header$.commentsOffset <- readBin(con=con, what=double(0), size=8);
  offset <- offset + 8;

  # Number of dimension (double = 8 bytes)
  seek(con, where=offset, rw="read");
  nbrOfDims <- readBin(con=con, what=double(0), size=8);
  offset <- offset + 8;

  # Length of each dimension
  dim <- vector("double", nbrOfDims);
  for (kk in seq(length=nbrOfDims)) {
    seek(con, where=offset, rw="read");
    dim[kk] <- readBin(con=con, what=double(0), size=8);
    offset <- offset + 8;
  }
  header$dim <- dim;

  # Order of dimensions
  seek(con, where=offset, rw="read");
  header$dimOrder <- readBin(con=con, what=double(0), n=nbrOfDims, size=8);
  offset <- offset + 8 * nbrOfDims;

  # Bytes per cell
  seek(con, where=offset, rw="read");
  header$bytesPerCell <- readBin(con=con, what=integer(0), size=4);
  offset <- offset + 4;

  # R storage mode
  seek(con, where=offset, rw="read");
  str <- readString(con=con);
  offset <- offset + attr(str, "count");
  header$storageMode <- as.character(str);

  # Dimension names
  dimnames <- vector("list", nbrOfDims);
  for (kk in seq(length=nbrOfDims)) {
    # Number of names for dimension 'kk' (double = 8 bytes)
    seek(con, where=offset, rw="read");
    nbrOfNames <- readBin(con=con, what=double(0), size=8);
    offset <- offset + 8;

    names <- vector("character", nbrOfNames);
    for (ll in seq(length=nbrOfNames)) {
      seek(con, where=offset, rw="read");
      str <- readString(con=con);
      offset <- offset + attr(str, "count");
      names[ll] <- str;
    }
    dimnames[[kk]] <- names;
    names <- NULL; # Not needed anymore
  }
  header$dimnames <- dimnames;
  dimnames <- NULL; # Not needed anymore

  # Comments
  seek(con, where=offset, rw="read");
  str <- readString(con=con);
  offset <- offset + attr(str, "count");
  header$comments <- as.character(str);

  # Number of free bytes before the data section
  header$nbrOfFreeBytes <- header$.dataOffset - offset;

  header;
}, protected=TRUE)  # readHeader()




###########################################################################/**
# @RdocMethod writeEmptyData
#
# @title "Writes an empty data section to the data file of a file array"
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
#  Returns nothing.
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
setMethodS3("writeEmptyData", "AbstractFileArray", function(this, ...) {
  con <- this$con;

  # 1. Move to the last file position.
  pos <- getDataOffset(this);
  pos <- pos + getSizeOfData(this);
  seek(con, where=pos-1, rw="write");

  # 2. Write an empty element
  writeBin(con=con, as.integer(0), size=1);
}, protected=TRUE)



###########################################################################/**
# @RdocMethod getStorageMode
#
# @title "Gets the storage mode of the file array"
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
setMethodS3("getStorageMode", "AbstractFileArray", function(this, ...) {
  this$header$storageMode;
})


###########################################################################/**
# @RdocMethod getBytesPerCell
#
# @title "Gets the number of bytes per element in a file array"
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
#  Returns an @integer.
# }
#
# \details{
#  In \R, an @integer
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
setMethodS3("getBytesPerCell", "AbstractFileArray", function(this, ...) {
  as.integer(this$header$bytesPerCell);
})



###########################################################################/**
# @RdocMethod getDataOffset
#
# @title "Gets file position of the data section in a file array"
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
#  Returns a @double (since the largest @integer is only 2147483648-1).
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
setMethodS3("getDataOffset", "AbstractFileArray", function(this, ...) {
  dataOffset <- this$header$.dataOffset;
  if (is.null(dataOffset)) {
    dataOffset <- readHeader(this)$.dataOffset;
    this$header$.dataOffset <- dataOffset;
  }
  dataOffset;
})



###########################################################################/**
# @RdocMethod getSizeOfData
#
# @title "Gets the size of the data section in bytes"
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
#  Returns a @double (since the largest @integer is only 2147483648-1).
# }
#
# @author
#
# \seealso{
#   @seemethod "getBytesPerCell".
#   @seemethod "length".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("getSizeOfData", "AbstractFileArray", function(this, ...) {
  getBytesPerCell(this)*length(this);
})


###########################################################################/**
# @RdocMethod getDimensionOrder
#
# @title "Gets the order of dimension"
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
#  Returns an @numeric @vector.
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
setMethodS3("getDimensionOrder", "AbstractFileArray", function(this, ...) {
  this$header$dimOrder;
})




###########################################################################/**
# @RdocMethod dim
#
# @title "Gets the dimension of the file array"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \value{
#  Returns a @double @vector of length two.
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
setMethodS3("dim", "AbstractFileArray", function(x) {
  # To please R CMD check
  this <- x;

  this$header$dim;
}, appendVarArgs=FALSE)


###########################################################################/**
# @RdocMethod length
#
# @title "Gets the number of elements in a file array"
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
#  Returns a @double (since the largest @integer is only 2147483648-1).
# }
#
# @author
#
# \seealso{
#   @seemethod "dim".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("length", "AbstractFileArray", function(x) {
  # To please R CMD check
  this <- x;

  prod(dim(this));
}, appendVarArgs=FALSE)


###########################################################################/**
# @RdocMethod dimnames
#
# @title "Gets the dimension names of a file array"
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
#  Returns a @character @vector.
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
setMethodS3("dimnames", "AbstractFileArray", function(x) {
  # To please R CMD check
  this <- x;

  this$header$dimnames;
}, appendVarArgs=FALSE)



###########################################################################/**
# @RdocMethod as.vector
#
# @title "Returns the elements of a file array as an R vector"
#
# \description{
#  @get "title", that is, imported into memory (if possible).
# }
#
# @synopsis
#
# \arguments{
#   \item{mode}{Not used.}
# }
#
# \value{
#  Returns a @vector.
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
setMethodS3("as.vector", "AbstractFileArray", function(x, mode="any") {
  # To please R CMD check
  this <- x;

  readAllValues(this);
}, appendVarArgs=FALSE)




###########################################################################/**
# @RdocMethod readAllValues
#
# @title "Reads all values in the file array"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{mode}{The storage mode to read.}
#   \item{size}{The number of bytes each values allocates on file.}
#   \item{offset}{The file offset to the first value on file.}
#   \item{...}{Not used.}
#   \item{.checkArgs}{If @TRUE, arguments are validated, otherwise not.}
# }
#
# \value{
#  Returns a @numeric @vector.
# }
#
# @author
#
# \seealso{
#   @see "readValues".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("readAllValues", "AbstractFileArray", function(this, mode=getStorageMode(this), size=getBytesPerCell(this), offset=getDataOffset(this), ..., .checkArgs=FALSE) {
  con <- this$con;

  .seekCon(con=con, where=offset, rw="read");

  what <- mode;
  n <- prod(this$header$dim);
## .Internal() calls are no longer allowed. /HB 2012-04-15
##  .Internal(readBin(con, what, n, size, TRUE, FALSE));
  readBin(con=con, what=what, n=n, size=size, signed=TRUE);
})




###########################################################################/**
# @RdocMethod readContiguousValues
#
# @title "Reads sets of contiguous values in the file array"
#
# \description{
#  @get "title".
#  A set of \emph{contiguous} values are values that are connecting without
#  a break.  It is much faster to read contiguous sequences at once than read
#  each value separately.
# }
#
# @synopsis
#
# \arguments{
#   \item{indices}{A @numeric @vector of start (first) indices in each of
#     the contiguous sets.}
#   \item{length}{A @numeric @vector specifying the length of each of the
#     contiguous sets.}
#   \item{mode}{The storage mode to read.}
#   \item{size}{The number of bytes each values allocates on file.}
#   \item{offset}{The file offset to the first value on file.}
#   \item{...}{Not used.}
#   \item{.checkArgs}{If @TRUE, arguments are validated, otherwise not.}
# }
#
# \value{
#  Returns a @numeric @vector.
# }
#
# @author
#
# \seealso{
#   @seemethod "readAllValues" and @see "readValues".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("readContiguousValues", "AbstractFileArray", function(this, indices, lengths=1, mode=getStorageMode(this), size=getBytesPerCell(this), offset=getDataOffset(this), ..., .checkArgs=TRUE) {
  con <- this$con;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  n <- length(this);
  if (.checkArgs) {
    # Argument 'indices':
    # (Allow double indices as well in case the 32-bit integer are to short)
    indices <- Arguments$getNumerics(indices, range=c(1,n));
    nbrOfIndices <- length(indices);

    # Argument 'lengths':
    n <- length(this);
    # (Allow double length as well in case the 32-bit integer are to short)
    lengths <- Arguments$getNumerics(lengths, range=c(1,n));
    nbrOfLengths <- length(lengths);
    if (nbrOfLengths < nbrOfIndices) {
      lengths <- rep(lengths, length.out=nbrOfIndices);
    } else if (nbrOfLengths < nbrOfIndices) {
      throw("Argument 'length' is longer than argument 'indices': ",
                                            nbrOfLengths, " > ", nbrOfIndices);
    }
  } else {
    nbrOfIndices <- length(indices);
  }

  # Nothing to do?
  if (nbrOfIndices == 0)
    return(values);

  # Calculate file offsets for all starting positions
  fileOffsets <- as.double(offset + (indices-1)*size);
  indices <- NULL; # Not needed anymore

  # Allocate vector for values;
  values <- vector(mode, sum(lengths));

  # Current position in vector of values.
  pos <- 0;

  # Call internal functions directly; less overhead (almost twice as fast)
  what <- mode;
  ## origin <- pmatch("start", c("start", "current", "end"));
  ## rw <- pmatch("read", c("read", "write"), 0);
  for (kk in seq(length=nbrOfIndices)) {
    # Move to start position
    # seek(con=con, where=fileOffsets[kk], rw="read");
    ## .Internal(seek(con, fileOffsets[kk], origin, rw));
    .seekCon(con=con, where=fileOffsets[kk], rw="read");

    # The number of values to read at this position
    n <- lengths[kk];

    # Where to store the values
    idx <- pos + 1:n;

    # Read the sequence of values at this position
    values[idx] <- readBin(con=con, what=what, n=n, size=size, signed=TRUE);

    # Next position
    pos <- pos + n;
  }
  fileOffsets <- NULL; # Not needed anymore

  values;
})


###########################################################################/**
# @RdocMethod readValues
#
# @title "Reads individual values in the file array"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{indices}{A @numeric @vector of indices.}
#   \item{mode}{The storage mode to read.}
#   \item{size}{The number of bytes each values allocates on file.}
#   \item{offset}{The file offset to the first value on file.}
#   \item{...}{Not used.}
#   \item{.checkArgs}{If @TRUE, arguments are validated, otherwise not.}
# }
#
# \value{
#  Returns a @numeric @vector.
# }
#
# @author
#
# \seealso{
#   @seemethod "readContiguousValues" and @see "readValues".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("readValues", "AbstractFileArray", function(this, indices=NULL, mode=getStorageMode(this), size=getBytesPerCell(this), offset=getDataOffset(this), order=FALSE, ..., .checkArgs=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (.checkArgs) {
    # Argument 'size':
    size <- Arguments$getInteger(size, range=c(1,8));

    # The number of elements to read
    if (is.null(indices)) {
      return(readAllValues(this, mode=mode, size=size, offset=offset, ...));
    }
  }

  con <- this$con;
  what <- mode;

  ni <- length(indices);

  # Ordered reading for faster access
  res <- sort(indices, index.return=TRUE);
  o <- order(res$ix);
  indices <- res$x;
  res <- NULL; # Not needed anymore

  # Allocate array for return value
  mode <- getStorageMode(this);
  values <- vector(mode=mode, length=ni);

  fileOffsets <- as.double(offset + (indices-1)*size);
  indices <- NULL; # Not needed anymore

  # Call internal functions directly; less overhead (almost twice as fast)
  ## origin <- pmatch("start", c("start", "current", "end"));
  ## rw <- pmatch("read", c("read", "write"), 0);
  for (kk in seq(length=ni)) {
    # seek(con=con, where=fileOffsets[kk], rw="read");
    ## .Internal(seek(con, fileOffsets[kk], origin, rw));
    .seekCon(con=con, where=fileOffsets[kk], rw="read");
    values[kk] <- readBin(con=con, what=what, n=1L, size=size, signed=TRUE);
  }
  fileOffsets <- NULL; # Not needed anymore

  values <- values[o];
  o <- NULL; # Not needed anymore

  values;
})



###########################################################################/**
# @RdocMethod writeAllValues
#
# @title "Writes all values to a file array"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{values}{A @numeric @vector of values to be assigned.}
#   \item{mode}{The storage mode to be used.}
#   \item{size}{A @integer specifying the number of bytes per cell.}
#   \item{offset}{A @integer specifying the file offset of the data section.}
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
#   @seemethod "writeValues".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("writeAllValues", "AbstractFileArray", function(this, values, mode=getStorageMode(this), size=getBytesPerCell(this), offset=getDataOffset(this), ...) {
  con <- this$con;
  # Make sure values are of the correct data type.
  storage.mode(values) <- mode;

  # Allocate vector for values;
  seek(con=con, where=offset, rw="write");
  writeBin(con=con, values, size=size);
}, protected=TRUE)



###########################################################################/**
# @RdocMethod writeValues
#
# @title "Writes values to a file array"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{indices}{An @integer @vector of indices to be updated.
#     If @NULL, all cells are considered.}
#   \item{values}{A @numeric @vector of values to be assigned.}
#   \item{mode}{The storage mode to be used.}
#   \item{size}{A @integer specifying the number of bytes per cell.}
#   \item{offset}{A @integer specifying the file offset of the data section.}
#   \item{order}{If @TRUE, the data is reordered before being written.}
#   \item{...}{Additional arguments passed to @seemethod "writeAllValues".}
# }
#
# \value{
#   Returns nothing.
# }
#
# @author
#
# \seealso{
#   @seemethod "writeAllValues".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("writeValues", "AbstractFileArray", function(this, indices=NULL, values, mode=getStorageMode(this), size=getBytesPerCell(this), offset=getDataOffset(this), order=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##   # Call internal writeBin() to avoid overhead.
##   # From R v2.10.0
##   writeBinX <- function(object, con, size, ...) {
##     .Internal(writeBin(values[kk], con, size, FALSE, FALSE));
##   }

  # The number of elements to write
  if (is.null(indices)) {
    return(writeAllValues(this, values=values, mode=mode, size=size, offset=offset, ...));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'size':
  size <- Arguments$getInteger(size, range=c(1,8));

  con <- this$con;

  # The number of elements to read
  if (is.null(indices)) {
    indices <- seq(along=values);
  } else {
  }

  ni <- length(indices);

  # Ordered writing for faster access
  res <- sort(indices, index.return=TRUE);
  indices <- res$x;
  values <- values[res$ix];
  res <- NULL; # Not needed anymore

  # Make sure values are of the correct data type.
  storage.mode(values) <- mode;

  fileOffsets <- as.double(offset + (indices-1)*size);

  # Call internal functions directly; less overhead
  ## origin <- pmatch("start", c("start", "current", "end"));
  ## rw <- pmatch("write", c("read", "write"), 0);
  for (kk in seq(length=ni)) {
    # seek(con=con, where=fileOffsets[kk], rw="write");
    ## .Internal(seek(con, fileOffsets[kk], origin, rw));
    .seekCon(con=con, where=fileOffsets[kk], rw="write");
    ## writeBinX(values[kk], con=con, size=size);
    writeBin(values[kk], con=con, size=size, useBytes=FALSE);
  }
}, protected=TRUE)



############################################################################
# HISTORY:
# 2013-09-21
# o ROBUSTNESS/WORKAROUND: For now, package attaches the 'R.oo' package.
#   This is needed due to what appears to be a bug in how R.oo
#   finalizes Object:s assuming R.oo is/can be attached.  Until that
#   is resolved, we make sure R.oo is attached.
# 2012-04-15
# o Now no longer calling .Internal() readBin() and writeBin().
# o Now utilizing new .seekCon() instead of .Internal(seek(...)).
# o CLEANUP: readAllValues() for AbstractFileArray no longer uses
#   .Internal(seek(...)).
# o CLEANUP: Now as.vector() for AbstractFileArray uses the exact same
#   arguments as the base::as.vector() generic function.  This avoids
#   having to create a new one in R.huge, which R CMD check complaints
#   about.
# 2011-02-01
# o ROBUSTNESS: Now using argument 'nchars' (not 'nchar') when calling
#   readChar().
# 2009-09-22
# o BUG FIX: From R v2.10.0, writeValues() of AbstractFileArray would give
#   "Error: 4 arguments passed to .Internal(writeBin) which requires 5".
#   Updated so it works with all versions of R.
# 2009-05-19
# o ROBUSTNESS: Now open() of AbstractFileArray first tries to open the
#   file for reading and updating (as before).  If that fails, then it
#   tries to open the file for reading only.  It might be that the file
#   is only used for reading, so if the permission allows for that but
#   not updating, then open it.
# o EXCEPTION HANDLING: Methods that creates/modifies files will give
#   a clear error message if the file permissions does not allow it.
# 2006-08-21
# o BUG FIX: Argument 'dimOrder' of the AbstractFileArray constructor was
#   not recognized causing weird results if it was intended to be in a
#   non-increasing order, e.g. FileMatrix(..., byrow=TRUE).
# 2006-05-21
# o Added some more Rdoc comments.
# 2006-05-09
# o Added some more Rdoc comments.
# o Made a few methods protected.
# 2006-03-14
# o Added getComments() and setComments().
# o Added a default buffer of free bytes after the header comments.  This
#   will make it easier to update the comments.
# o Modified the file format of the header.
# 2006-03-04
# o Added readSeqsOfValues().
# 2006-02-27
# o Created from FileMatrix.R.
############################################################################
