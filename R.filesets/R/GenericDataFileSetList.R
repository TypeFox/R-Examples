###########################################################################/**
# @RdocClass GenericDataFileSetList
#
# @title "The GenericDataFileSetList class"
#
# \description{
#  @classhierarchy
#
#  A GenericDataFileSetList object represents a list of
#  @see "GenericDataFileSet"s.
# }
#
# @synopsis
#
# \arguments{
#   \item{dsList}{A single or a @list of @see "GenericDataFileSet":s.}
#   \item{tags}{A @character @vector of tags.}
#   \item{...}{Not used.}
#   \item{allowDuplicates}{If @FALSE, files with duplicated names are not
#     allowed and an exception is thrown, otherwise not.}
#   \item{.setClass}{A @character string specifying a name of the
#     class that each data set must be an instance of.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @examples "../incl/GenericDataFileSetList.Rex"
#
# @author
#*/###########################################################################
setConstructorS3("GenericDataFileSetList", function(dsList=list(), tags="*", ..., allowDuplicates=TRUE, .setClass="GenericDataFileSet") {
  # Argument '.setClass':
  .setClass <- Arguments$getCharacter(.setClass, length=c(1,1));

  # Argument 'dsList':
  if (!is.list(dsList)) {
    dsList <- list(dsList);
  }

  if (is.list(dsList)) {
    for (ds in dsList) {
      ds <- Arguments$getInstanceOf(ds, .setClass);
    }
  } else {
    throw("Argument 'dsList' is not a list: ", class(dsList)[1]);
  }

  # Arguments 'allowDuplicates':
  allowDuplicates <- Arguments$getLogical(allowDuplicates);

  this <- extend(Object(), c("GenericDataFileSetList", uses("FullNameInterface")),
    .dsList = dsList,
    .tags = NULL,
    .allowDuplicates = allowDuplicates
  );

  setTags(this, tags);

  assertDuplicates(this);

  this;
})

setMethodS3("assertDuplicates", "GenericDataFileSetList", function(this, ...) {
  allowDuplicates <- this$.allowDuplicates;

  if (!allowDuplicates) {
    names <- getNames(this);
    dups <- names[duplicated(names)];
    n <- length(dups);
    if (n > 0) {
      throw(sprintf("Detected %n files with duplicated names, which are not allowed (onDuplicates=FALSE): %s", n, hpaste(dups)));
    }
  }
}, protected=TRUE)



setMethodS3("as.GenericDataFileSetList", "GenericDataFileSetList", function(this, ...) {
  # Nothing to do
  this;
})



setMethodS3("as.character", "GenericDataFileSetList", function(x, ...) {
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

  # Unique names
  names <- getNames(this);
  n <- length(names);
  s <- c(s, sprintf("Names: %s [%d]", hpaste(names), n));

  # Number of file sets
  n <- nbrOfSets(this);
  s <- c(s, sprintf("Number of file sets: %d", n));

  for (kk in seq_len(n)) {
    ds <- getSet(this, kk);
    s <- c(s, sprintf("<File set #%d ('%s') of %d>", kk, getName(ds), n));
    s <- c(s, as.character(ds));
  }

  GenericSummary(s);
}, protected=TRUE)


setMethodS3("clone", "GenericDataFileSetList", function(this, ...) {
  res <- NextMethod("clone");
  dsList <- getSets(this);
  for (kk in seq_along(dsList)) {
    dsList[[kk]] <- clone(dsList[[kk]]);
  }
  res$.dsList <- dsList;
  res;
}, protected=TRUE)


setMethodS3("getAsteriskTags", "GenericDataFileSetList", function(this, ..., collapse=",") {
  # Get the tags of each data set
  dsList <- getSets(this);
  tags <- sapply(dsList, getTags, collapse=",");
  tags <- tags[nchar(tags, type="chars") > 0L];
  tags <- paste(tags, collapse=collapse);
  tags;
}, protected=TRUE)


setMethodS3("setTags", "GenericDataFileSetList", function(this, tags=NULL, ...) {
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));
    tags <- tags[nchar(tags, type="chars") > 0L];
  }

  this$.tags <- tags;

  invisible(this);
}, protected=TRUE)


setMethodS3("getTags", "GenericDataFileSetList", function(this, ...) {
  # AD HOC; custom tags are taken care of in getDefaultFullName().
  NextMethod("getTags", collapse=NULL, useCustomTags=FALSE);
}, protected=TRUE)


setMethodS3("getDefaultFullName", "GenericDataFileSetList", function(this, collapse="+", ...) {
  # Get the names of each data set
  dsList <- getSets(this);
  names <- sapply(dsList, getName);

  # By default, merge names
  name <- mergeByCommonTails(names, collapse=collapse);

  # Get the tags
  tags <- this$.tags;
  tags[tags == "*"] <- getAsteriskTags(this, collapse=",");
  tags <- tags[nchar(tags, type="chars") > 0L];

  fullname <- paste(c(name, tags), collapse=",");

  fullname;
}, protected=TRUE)


setMethodS3("as.list", "GenericDataFileSetList", function(x, ...) {
  # To please R CMD check
  this <- x;

  getSets(this, ...);
})


setMethodS3("getSets", "GenericDataFileSetList", function(this, idxs=NULL, ...) {
  res <- this$.dsList;
  if (is.null(idxs)) {
  } else {
    n <- length(res);
    idxs <- Arguments$getIndices(idxs, max=n);
    res <- res[idxs];
  }
  res;
})

setMethodS3("getSet", "GenericDataFileSetList", function(this, idx, ...) {
  if (length(idx) != 1L)
    throw("Argument 'idx' must be a single index.");
  res <- this$.dsList;
  n <- length(res);
  idx <- Arguments$getIndex(idx, max=n);
  res[[idx]];
})

setMethodS3("nbrOfSets", "GenericDataFileSetList", function(this, ...) {
  length(getSets(this));
})


setMethodS3("getFileListClass", "GenericDataFileSetList", function(this, ...) {
  classNames <- class(this);
  pattern <- "(.*)(List|Tuple)$";
  exts <- gsub(pattern, "\\2", classNames);
  keep <- is.element(exts, c("List", "Tuple"));
  classNames <- classNames[keep];
  exts <- exts[keep];
  classNames <- gsub("FileSet", "Set", classNames);
  classNames <- gsub("Set", "File", classNames);

  # Keep only existing class names
  keep <- sapply(classNames, FUN=exists, mode="function");
  classNames <- classNames[keep];

  if (length(classNames) == 0L) {
    throw("Failed to locate a file list class for this set list: ",
                                                      class(this)[1L]);
  }

  className <- classNames[1L];

  clazz <- Class$forName(className, envir=parent.frame());
  classNames <- c(getKnownSubclasses(clazz), className);
  clazz <- NULL;
  for (kk in seq_along(classNames)) {
     className <- classNames[kk];
     tryCatch({
       clazz <- Class$forName(className, envir=parent.frame());
     }, error = function(ex) {});
     if (!is.null(clazz)) {
       return(className);
     }
  } # for (kk ...)

  throw("Failed to locate a file list class for this set list: ",
                                                      class(this)[1]);
}, protected=TRUE)



setMethodS3("getFullNames", "GenericDataFileSetList", function(this, ...) {
  dsList <- getSets(this);
  names <- lapply(dsList, FUN=getFullNames);
  names <- unlist(names, use.names=FALSE);
  names <- unique(names);
  names <- sort(names);
  names;
})


setMethodS3("getNames", "GenericDataFileSetList", function(this, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  unionWithDuplicates <- function(x, y) {
    keep <- is.element(x, y);
    x <- x[keep];
    keep <- is.element(y, x);
    y <- y[keep];
    z <- if(length(x) > length(y)) x else y;
    sort(z);
  }

  dsList <- getSets(this);
  names <- lapply(dsList, FUN=getNames);

  # Note, we cannot use union() if there may be duplicates
  names <- Reduce(unionWithDuplicates, names);

  names <- unlist(names, use.names=FALSE);

  # Return unique names?
  allowDuplicates <- this$.allowDuplicates;
  unique <- (!allowDuplicates);
  if (unique) {
    names <- unique(names);
  }

  names <- sort(names);
  names;
})


setMethodS3("length", "GenericDataFileSetList", function(x, ...) {
  # To please R CMD check
  this <- x;
  names <- getNames(this, ...);
  length(names);
})

setMethodS3("nbrOfFiles", "GenericDataFileSetList", function(this, ...) {
  length(this, ...);
}, protected=TRUE)


setMethodS3("indexOf", "GenericDataFileSetList", function(this, ...) {
  # Since indexOf() for GenericDataFileSetList were identical to that
  # of GenericDataFileSet, we call the latter explicitly.
  indexOf.GenericDataFileSet(this, ...);
})



setMethodS3("extract", "GenericDataFileSetList", function(this, files, ..., drop=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'files':
  nbrOfFiles <- length(this);
  if (is.null(files)) {
    files <- seq_len(nbrOfFiles);
  } else if (is.logical(files)) {
    files <- Arguments$getVector(files, length=rep(nbrOfFiles, 2));
    files <- which(files);
  } else if (is.character(files)) {
    files <- indexOf(this, files, ...);
  } else if (is.numeric(files)) {
    files <- Arguments$getIntegers(files, disallow="NaN");

    # Exclude indices?
    if (any(files < 0, na.rm=TRUE)) {
      incl <- files[files > 0L];
      if (length(incl) == 0L) {
        incl <- seq_along(this);
      }
      excl <- na.omit(files[files < 0L]);
      files <- setdiff(incl, -excl);
      # Not needed anymore
      incl <- excl <- NULL;
    }
  }

  # Avoid return duplicates if not allowed
  allowDuplicates <- this$.allowDuplicates;
  if (!allowDuplicates) {
    dups <- files[duplicated(files)];
    n <- length(dups);
    if (n > 0L) {
      throw(sprintf("Argument 'files' contains %s duplicates, which is not allowed (allowDuplicates=FALSE): %s", n, hpaste(dups)));
    }
  }

  names <- getNames(this);
  files <- Arguments$getIndices(files, max=nbrOfFiles, disallow="NaN");
  files <- names[files];

  dsList <- getSets(this);
  setNames <- sapply(dsList, FUN=getFullName);
  dsList <- lapply(dsList, FUN=extract, files, ...);
  names(dsList) <- setNames;

  # Drop empty data sets?
  if (drop) {
    ns <- sapply(dsList, nbrOfFiles);
    dsList <- dsList[ns > 0L];
  }

  res <- clone(this);
  res$.dsList <- dsList;

  assertDuplicates(res);

  res;
}, protected=TRUE)



setMethodS3("as.data.frame", "GenericDataFileSetList", function(x, row.names=NULL, ..., names=NULL, onDuplicates=c("drop", "error")) {
  # To please R CMD check
  this <- x;

  # Argument 'names':
  if (is.null(names)) {
    names <- getNames(this);
  } else {
    names <- Arguments$getCharacters(names);
  }

  # Argument 'row.names':
  if (is.null(row.names)) {
    if (any(is.na(names))) {
      row.names <- seq_along(names);
    } else {
      row.names <- names;
    }
  } else {
    row.names <- Arguments$getCharacters(row.names,
                           length=rep(length(names), times=2L));
  }

  # Argument 'onDuplicates':
  onDuplicates <- match.arg(onDuplicates);


  # Extract subset of interests (missing files ok; duplicates too)
  if (is.element(onDuplicates, c("drop", "error"))) {
    onDuplicates2 <- onDuplicates;
  } else {
    onDuplicates2 <- "ignore";
    # Not supported
    throw("Unsupported value of 'onDuplicates': ", onDuplicates);
  }

  res <- extract(this, names, onMissing="NA", onDuplicates=onDuplicates2, ...);

  # Get list of files
  dsList <- getSets(res);
  dfList <- lapply(dsList, FUN=getFiles);
  # Not needed anymore
  res <- dsList <- NULL;

  # Sanity check
  ns <- sapply(dfList, FUN=length);
  if (length(unique(ns)) != 1L) {
    throw("Internal error. Cannot extract data frame. Non-balanced data sets.");
  }

  attr(dfList, "row.names") <- row.names;
  attr(dfList, "class") <- "data.frame";

  dfList;
}, protected=TRUE)


setMethodS3("getFileList", "GenericDataFileSetList", function(this, ..., dropMissing=TRUE) {
  dsList <- extract(this, ..., onDuplicated="error");
  dfList <- lapply(dsList, FUN=getFile, 1L);

  if (dropMissing) {
    keep <- sapply(dfList, FUN=isFile);
    dfList <- dfList[keep];
  }

  # Coerce to a file list
  className <- getFileListClass(this);
  clazz <- Class$forName(className, envir=parent.frame());
  dfList <- newInstance(clazz, dfList);

  dfList;
})


###########################################################################
# HISTORY:
# 2012-11-12
# o CLEANUP: Now nbrOfFiles() for GenericDataFileSetList calls length()
#   and not vice versa.  This will make it easier to phase out
#   nbrOfFiles() in the future.
# 2010-07-06
# o Now indexOf() for GenericDataFileSetList is identical to that of
#   GenericDataFileSet.
# 2010-02-07
# o BUG FIX: indexOf() of GenericDataFileSetList did not handle names with
#   regular expression symbols '+' and '*'.
# 2010-01-24
# o ROBUSTNESS: If argument 'files' is logical, then extract() of
#   GenericDataFileSetList now asserts that the length of 'files' matches
#   the number of available files.
# 2010-01-01
# o BUG FIX: getNames() of GenericDataFileSetList would drop duplicated
#   names if there where more than one data set.
# o Added constructor argument 'allowDuplicates=TRUE'.  If FALSE, an
#   error is thrown if a object with duplicates are tried to be created,
#   either via the constructor or via extract().
# 2009-12-30
# o Rename getFileClass() to getFileListClass().
# o Replaced getListOfSets() with getSets(), cf. getFiles().
# o Created new getFileList().
# o Turned into an Object().
# o Renamed argument '.fileSetClass' to '.setClass'.
# o Now GenericDataFileSetList is a FullNameInterface class.
# o Added Rdoc comments for GenericDataFileSetList.
# o Added protected as.data.frame() to GenericDataFileSetList.
# o Added protected extract().
# 2009-06-03
# o Added argument '.fileSetClass' to the constructor.
# 2009-05-12
# o Created.
###########################################################################
