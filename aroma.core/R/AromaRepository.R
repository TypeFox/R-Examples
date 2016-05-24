###########################################################################/**
# @RdocClass AromaRepository
#
# @title "The AromaRepository class"
#
# \description{
#  @classhierarchy
#
#  An AromaRepository object provides methods for downloading annotation data
#  from the Aroma repository.
# }
#
# @synopsis
#
# \arguments{
#   \item{urlPath}{The URL to the Aroma reposity.}
#   \item{verbose}{The @see "R.utils::Verbose" to be used during processing.}
#   \item{...}{Not used.}
# }
#
# \section{Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("AromaRepository", function(urlPath="http://www.aroma-project.org/data", verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'urlPath':
  urlPath <- Arguments$getCharacter(urlPath);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose, ...);

  extend(Object(), "AromaRepository",
    .urlPath = urlPath,
    .verbose = verbose
  );
})

setMethodS3("getUrlPath", "AromaRepository", function(this, ...) {
  this$.urlPath;
}, protected=TRUE)


setMethodS3("setVerbose", "AromaRepository", function(this, ...) {
  verbose <- Arguments$getVerbose(verbose, ...);
  this$.verbose <- verbose;
  invisible(this);
}, protected=TRUE)


setMethodS3("getVerbose", "AromaRepository", function(this, ...) {
  this$.verbose;
}, protected=TRUE)


setMethodS3("clearCache", "AromaRepository", function(this, ...) {
  dirs <- c("aroma.core", "AromaRepository", as.character(Sys.Date()));
  path <- R.cache::getCachePath(dirs);
  removeDirectory(path, recursive=TRUE, mustExist=FALSE);
}, protected=TRUE)




###########################################################################/**
# @RdocMethod listFiles
#
# @title "Retrieves the files available on the repository under a particular path"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{path}{The relative path to be listed.}
#   \item{orderBy}{A @character specifying how the returns files should be ordered.}
#   \item{...}{Additional arguments passed to @see "R.utils::downloadFile".}
#   \item{force}{If @FTRUE, cached results are redownloaded, otherwise not.}
# }
#
# \value{
#  Returns the relative pathnames of the files available.
# }
#
# \details{
#   Note that this method makes strong assumptions of the format of the
#   downloaded HTML index file.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("listFiles", "AromaRepository", function(this, path=NULL, full=TRUE, orderBy=c("name", "time"), ..., force=FALSE, verbose=getVerbose(this)) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'orderBy':
  orderBy <- match.arg(orderBy);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Get the URL to download
  urlPath <- getUrlPath(this);
  if (!is.null(path)) {
    urlPath <- file.path(urlPath, path);
  }

  verbose && enter(verbose, "Listing files");

  verbose && cat(verbose, "URL to download: ", urlPath);

  dirs <- c("aroma.core", "AromaRepository", as.character(Sys.Date()));
  key <- list(method="downloadListFiles", class=class(this), urlPath=urlPath, full=full, orderBy=orderBy);
  res <- loadCache(key=key, dirs=dirs);
  if (!force && !is.null(res)) {
    verbose && enter(verbose, "Available files:");
    verbose && print(verbose, res);
    verbose && exit(verbose);
    return(res);
  }

  # Download the URL HTML index file
  filename <- tempfile();
  on.exit(file.remove(filename));
  pathname <- tryCatch({
    suppressWarnings({
      downloadFile(url=urlPath, filename=filename, ...);
    });
  }, error = function(ex) {
    # Failed to download index file.  Assume directory does not exists.
    return(NULL);
  });
  if (is.null(pathname)) {
    verbose && exit(verbose);
    return(NULL);
  }

  # Parse the index file
  bfr <- readLines(pathname);

  # Extract the filenames
  pattern <- ".*<a href=\"([^\"/?][^\"]*)\">.*";
  bfr <- grep(pattern, bfr, value=TRUE);
  filenames <- gsub(pattern, "\\1", bfr);

  # Reorder?
  if (orderBy == "name") {
    o <- order(filenames);
  } else if (orderBy == "time") {
    patternT <- ".*>([0-9]*-[^-]*-[0-9]*)[ ]+([0-9]*:[0-9]*)[ ]*<.*";
    timestamps <- gsub(patternT, "\\1 \\2", bfr);
    timestamps <- strptime(timestamps, format="%d-%b-%Y %H:%M");
    o <- order(timestamps);
  }
  filenames <- filenames[o];

  if (full) {
    path <- gsub("/$", "", path);
    filenames <- file.path(path, filenames);
  }

  saveCache(filenames, key=key, dirs=dirs);

  verbose && enter(verbose, "Available files:");
  verbose && print(verbose, filenames);

  verbose && exit(verbose);

  filenames;
}, protected=TRUE) # listFiles()



###########################################################################/**
# @RdocMethod downloadFile
#
# @title "Download a particular file from the reposity"
#
# \description{
#  @get "title" by its relative pathname.
# }
#
# @synopsis
#
# \arguments{
#   \item{filename, path}{The filename and the relative path of the file
#         to be download.}
#   \item{gzipped}{If @TRUE, a gzipped file is downloaded and decompressed.}
#   \item{skip}{If @TRUE, an already downloaded file is skipped.}
#   \item{overwrite}{If @TRUE, an not skipping, an already downloaded file
#      is overwritten, otherwise an error is thrown.}
#   \item{...}{Additional arguments passed to @see "R.utils::downloadFile".}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns the pathname of the uncompressed downloaded file.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("downloadFile", "AromaRepository", function(this, filename, path=NULL, caseSensitive=FALSE, gzipped=TRUE, skip=TRUE, overwrite=FALSE, ..., verbose=getVerbose(this)) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'filename' & 'path':
  if (is.null(path)) {
    pathname <- filename;
  } else {
    pathname <- file.path(path, filename);
  }
  pathnameL <- Arguments$getWritablePathname(pathname, mustNotExist=!skip & !overwrite);

  # Argument 'caseSensitive':
  caseSensitive <- Arguments$getLogical(caseSensitive);

  # Argument 'gzipped':
  gzipped <- Arguments$getLogical(gzipped);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Downloading file");

  verbose && cat(verbose, "Local pathname: ", pathnameL);
  if (skip && isFile(pathnameL)) {
    verbose && cat(verbose, "Already downloaded: ", pathnameL);
    verbose && exit(verbose);
    return(pathnameL);
  }

  # If decompressed, check if already downloaded
  if (gzipped) {
    pathnameD <- sprintf("%s.gz", pathname);
  } else {
    pathnameD <- pathname;
  }
  pathnameDL <- Arguments$getWritablePathname(pathnameD,
                                      mustNotExist=!skip & !overwrite);

  verbose && cat(verbose, "File to download: ", pathnameD);

  # The filename and the relative path of the file to be downloaded
  filenameD <- basename(pathnameD);
  path <- dirname(pathnameD);

  # Get the list of files available for download
  pathnames <- listFiles(this, path=path, verbose=less(verbose,1));
#  verbose && cat(verbose, "Available files:");
#  verbose && print(verbose, pathnames);

  # Is the file available for download?
  if (!any(is.element(c(pathname, pathnameD), pathnames))) {
    msg <- paste("File not available for download: ", pathnameD, sep="");
    verbose && cat(verbose, msg);
    warning(msg);
    verbose && exit(verbose);
    return(NULL);
  }

  # Try to download the file
  urlPath <- getUrlPath(this);
  tryCatch({
    url <- file.path(urlPath, pathnameD);
    suppressWarnings({
      pathnameD <- downloadFile(url, filename=pathnameDL, skip=skip, overwrite=overwrite, ..., verbose=less(verbose,5));
    });
  }, error = function(ex) {
    # If gzipped file did not exists, try the regular one
    verbose && cat(verbose, "Failed to download compressed file. The reason was: ", ex$message);
    if (gzipped) {
      verbose && enter(verbose, "Trying to download non-compressed file");
      url <- file.path(urlPath, pathname);
      verbose && cat(verbose, "URL to download: ", url);
      pathname <- downloadFile(url, filename=pathnameL, skip=skip, overwrite=overwrite, ..., verbose=less(verbose,5));
      gzipped <<- FALSE;
      verbose && exit(verbose);
    } else {
      throw(ex);
    }
  })

  if (gzipped) {
    verbose && enter(verbose, "Decompressing file");
    gunzip(pathnameDL, overwrite=overwrite, remove=TRUE);
    verbose && exit(verbose);
  }

  # Sanity check
  stopifnot(isFile(pathname));

  verbose && exit(verbose);

  pathname;
}, protected=TRUE) # downloadFile()






setMethodS3("findAnnotationDataByChipType", "AromaRepository", function(this, chipType, tags=NULL, pattern=NULL, firstOnly=TRUE, ..., verbose=getVerbose(this)) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  orderByFullName <- function(pathnames, ..., verbose=FALSE) {
    # Nothing to do?
    if (length(pathnames) == 0) {
      return(integer(0));
    }

    verbose && enter(verbose, "Ordering in increasing lengths of fullnames");

    # Order located pathnames in increasing length of the fullnames
    # This is an AD HOC solution for selecting GenomeWideSNP_6 before
    # GenomeWideSNP_6,Full.

    # (a) Get filenames
    filenames <- basename(pathnames);

    # (b) Get fullnames by dropping filename extension
    fullnames <- gsub("[.][^.]*$", "", filenames);

    # (c) Order by length of fullnames
    o <- order(nchar(fullnames));

    verbose && cat(verbose, "Order:");
    verbose && print(verbose, o);

    verbose && exit(verbose);

    o;
  } # orderByFullNames()


  sortByFullName <- function(pathnames, ..., verbose=FALSE) {
    o <- orderByFullName(pathnames, verbose=verbose);
    pathnames <- pathnames[o];
    pathnames;
  } # sortByFullName()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(pattern)) {
    pattern <- Arguments$getRegularExpression(pattern);
  }

  # Argument 'firstOnly':
  firstOnly <- Arguments$getLogical(firstOnly);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Finding annotation data file by chip type");
  chipTypeS <- gsub(",.*", "", chipType);
  path <- file.path("annotationData", "chipTypes", chipTypeS);
  filenames <- listFiles(this, path=path, full=FALSE, verbose=less(verbose,1));

  # Drop directories
  filenames <- grep("/$", filenames, value=TRUE, invert=TRUE);

  # Sort files
  filenames <- sortByFullName(filenames);

  if (!is.null(pattern)) {
    filenames <- grep(pattern, filenames, value=TRUE);
  }

  if (firstOnly && length(filenames) > 1) {
    filenames <- filenames[1L];
  }

  # Return full pathnames
  pathnames <- file.path(path, filenames);

  verbose && exit(verbose);

  pathnames;
}, protected=TRUE) # findAnnotationDataByChipType()



###########################################################################/**
# @RdocMethod downloadChipTypeFile
#
# @title "Download a particular chip type annotation file"
#
# \description{
#  @get "title" by its chip type, tags and suffix.
# }
#
# @synopsis
#
# \arguments{
#   \item{chipType}{The chip type of the file to be downloaded.}
#   \item{tags}{Optional tags of the file to be downloaded.}
#   \item{suffix}{The filename suffix (including any preceeding period) of
#      the file to be downloaded.}
#   \item{ext}{The filename extension.}
#   \item{...}{Additional arguments passed to @seemethod "downloadFile".}
#   \item{skip}{If @TRUE, an already downloaded file is skipped.}
#   \item{overwrite}{If @TRUE, an not skipping, an already downloaded file
#      is overwritten, otherwise an error is thrown.}
#   \item{mustExist}{If @TRUE, an exception is thrown if no file matching
#      is available either locally or on the repository.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns the relative pathname of the uncompressed downloaded file.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("downloadChipTypeFile", "AromaRepository", function(this, chipType, tags=NULL, suffix=sprintf(".%s", ext), ext=NULL, ..., gunzip=TRUE, skip=TRUE, overwrite=FALSE, mustExist=TRUE, verbose=getVerbose(this)) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'chipType':
  chipType <- Arguments$getCharacter(chipType);

  # Argument 'tags':
  tags <- Arguments$getTags(tags);

  # Argument 'suffix':
  suffix <- Arguments$getCharacter(suffix);

  # Argument 'mustExist':
  mustExist <- Arguments$getLogical(mustExist);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Downloading chiptype file");

  chipTypeS <- gsub(",.*", "", chipType);

  path <- file.path("annotationData", "chipTypes", chipTypeS);
  verbose && cat(verbose, "Path: ", path);

  chipTypeF <- paste(c(chipType, tags), collapse=",");
  verbose && cat(verbose, "Full chip type: ", chipTypeF);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (a) Is file available on local file system?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Searching local file system");
  pattern <- sprintf("^%s(%s)$", chipTypeF, suffix);
  pathname <- findAnnotationDataByChipType(chipType, pattern=pattern, verbose=less(verbose,10));
  if (!is.null(pathname)) {
    verbose && cat(verbose, "Already downloaded: ", pathname);
    verbose && exit(verbose);
    return(pathname);
  }
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (b) Is file available on repository?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Searching repository");
  pattern <- sprintf("^%s(%s)(|.gz)$", chipTypeF, suffix);
  pathnameR <- findAnnotationDataByChipType(this, chipType, pattern=pattern, verbose=less(verbose,10));
  if (length(pathnameR) == 0) {
    msg <- sprintf("No such file available (with or without *.gz): %s/%s", path, pattern);
    verbose && cat(verbose, msg);
    if (mustExist) {
      throw(msg);
    }
    return(NULL);
    verbose && exit(verbose);
  }

  verbose && cat(verbose, "File found: ", pathnameR);
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (c) Download
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Downloading");
  pathname <- downloadFile(this, filename=pathnameR, ..., verbose=less(verbose,1));
  verbose && cat(verbose, "Downloaded file: ", pathname);

  # Gunzip
  if (gunzip) {
    if (regexpr("[.]gz$", pathname) != -1) {
      gunzip(pathname, remove=TRUE);
      pathname <- gsub("[.]gz$", "", pathname);
    }
  }

  verbose && exit(verbose);

  pathname;
}) # downloadChipTypeFile()



setMethodS3("downloadAll", "AromaRepository", function(this, ...) {
  suffixes <- c("acc", "acm", "acp", "acs", "(cdf|CDF)", "ufl", "ugp");
  suffixes <- sprintf("[.]%s", suffixes);
  pathnames <- lapply(suffixes, FUN=function(suffix) {
    downloadChipTypeFile(this, ..., suffix=suffix);
  });
  names(pathnames) <- suffixes;
  pathnames;
})


setMethodS3("downloadACC", "AromaRepository", function(this, ...) {
  downloadChipTypeFile(this, ..., ext="acc");
})

setMethodS3("downloadACM", "AromaRepository", function(this, ...) {
  downloadChipTypeFile(this, ..., ext="acm");
})

setMethodS3("downloadACP", "AromaRepository", function(this, ...) {
  downloadChipTypeFile(this, ..., ext="acp");
})

setMethodS3("downloadACS", "AromaRepository", function(this, ...) {
  downloadChipTypeFile(this, ..., ext="acs");
})

setMethodS3("downloadCDF", "AromaRepository", function(this, ...) {
  downloadChipTypeFile(this, ..., suffix="[.](cdf|CDF)");
})

setMethodS3("downloadUFL", "AromaRepository", function(this, ...) {
  downloadChipTypeFile(this, ..., ext="ufl");
})

setMethodS3("downloadUGP", "AromaRepository", function(this, ...) {
  downloadChipTypeFile(this, ..., ext="ugp");
})

setMethodS3("downloadUGC", "AromaRepository", function(this, ...) {
  downloadChipTypeFile(this, ..., ext="ugc");
})

setMethodS3("downloadTXT", "AromaRepository", function(this, ...) {
  downloadChipTypeFile(this, ..., ext="txt");
})

setMethodS3("downloadProbeSeqsTXT", "AromaRepository", function(this, ...) {
  downloadChipTypeFile(this, ..., suffix=",probeSeqs.txt");
})


######################################################################
# HISTORY:
# 2012-09-14
# o Now downloadFile() and listFiles() for AromaRepository no
#   longer gives a warning if the file does not exists.
# 2012-09-04
# o Added downloadUGC() for AromaRepository.
# 2012-08-31
# o Added findAnnotationDataByChipType() for AromaRepository.
# o Added argument 'mustExist=TRUE' to downloadChipTypeFile().
# o GENERALIZATION: Now downloadCDF() download both *.cdf and *.CDF.
# 2012-08-22
# o Added clearCache().
# o Added downloadAll().
# o Now all methods are non-static.
# o Now downloadFile() precheck which files are available, by
#   querying listFiles().
# o Added listFiles() for AromaReposity.
# o Added some Rdoc help.
# 2011-09-29
# o The purpose of this class is to simplify downloading of
#   data files needed in the Aroma Framework.
# o Created.
######################################################################
