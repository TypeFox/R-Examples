###########################################################################/**
# @RdocDefault findAnnotationData
#
# @title "Locates an annotation data file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{name}{Optional @character string.}
#   \item{tags}{Optional @character string.
#     Only used if argument \code{pattern} is not specified.}
#   \item{pattern}{A filename pattern to search for.
#     If @NULL, then defaults to the fullname as defined by
#     arguments \code{name} and \code{tags}.}
#   \item{private}{If @FALSE, files and directories starting
#     with a periods are ignored.}
#   \item{escapes}{A @character @vector specify symbols to be escaped
#     in argument \code{pattern}.}
#   \item{...}{Arguments passed to @see "R.utils::findFiles".}
#   \item{firstOnly}{If @TRUE, only the first matching pathname is returned.}
#   \item{paths}{A @character @vector of paths to search.
#     If @NULL, default paths are used.}
#   \item{set}{A @character string specifying what type of annotation
#     to search for.}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns @NULL, one or several matching pathnames.
# }
#
# @author
#
# @keyword internal
#*/###########################################################################
setMethodS3("findAnnotationData", "default", function(name=NULL, tags=NULL, set, pattern=NULL, private=FALSE, escapes=c("+"), ..., firstOnly=TRUE, paths=NULL, verbose=FALSE) {
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


  localFindFiles <- function(..., paths=paths) {
    isInPrivateDirectory <- function(pathname) {
      pathname <- strsplit(pathname, split="[/\\\\]")[[1]];
      pathname <- pathname[!(pathname %in% c(".", ".."))];
      any(regexpr("^[.]", pathname) != -1);
    } # isInPrivateDirectory()

    # By explicitly searching each path seperately as here, we
    # can guarantee that the root paths are search in the correct
    # order according to the aroma search conventions. /HB 2011-03-03
    pathnames <- c();
    for (kk in seq_along(paths)) {
      path <- paths[kk];
      verbose && enter(verbose, sprintf("Path #%d of %d", kk, length(paths)));

      verbose && cat(verbose, "Path: ", path);
      pathnamesKK <- findFiles(..., paths=path);
      verbose && print(verbose, pathnamesKK);

      if (length(pathnamesKK) == 0) {
        verbose && exit(verbose);
        next;
      }

      # AD HOC: Clean out files in "private" directories
      if (!private) {
        excl <- sapply(pathnamesKK, FUN=isInPrivateDirectory);
        pathnamesKK <- pathnamesKK[!excl];
        verbose && cat(verbose, "Dropping private directories:");
        verbose && print(verbose, pathnamesKK);
      }

      if (length(pathnamesKK) == 0) {
        verbose && exit(verbose);
        next;
      }

##      verbose && enter(verbose, "Reorder by fullname");
##      pathnamesKK <- sortByFullName(pathnamesKK, verbose=verbose);
##      verbose && print(verbose, pathnamesKK);
##      verbose && exit(verbose);

      pathnames <- c(pathnames, pathnamesKK);

      verbose && exit(verbose);
    } # for (kk ...)
    pathnames;
  } # localFindFiles()



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'name':
  if (!is.null(name)) {
    name <- Arguments$getCharacter(name, length=c(1,1));
  }

  # Argument 'tags':
  tags <- Arguments$getTags(tags);
  if (!is.null(tags) && is.null(name)) {
    throw("Argument 'tags' must be NULL if argument 'name' is: ", tags);
  }

  # Argument 'escapes':
  escapes <- Arguments$getCharacters(escapes);

  # Argument 'set':
  set <- Arguments$getCharacter(set, length=c(1,1));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Main
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(pattern)) {
    if (!is.null(name)) {
      fullname <- paste(c(name, tags), collapse=",");
      pattern <- fullname;
    }
  }

  # Escape pattern?
  if (length(escapes) > 0) {
    escapes <- unlist(strsplit(escapes, split="", fixed=TRUE), use.names=FALSE);
    escapes <- unique(escapes);
    escapePattern <- paste(escapes, collapse="");
    escapePattern <- sprintf("([%s])", escapePattern);
    pattern <- gsub(escapePattern, "[\\1]", pattern);
  }

  verbose && enter(verbose, "Searching for annotation data file(s)");

  verbose && cat(verbose, "Name: ", name);
  verbose && cat(verbose, "Tags: ", paste(tags, collapse=", "));
  verbose && cat(verbose, "Set: ", set);
  verbose && cat(verbose, "Filename pattern: ", pattern);
  verbose && cat(verbose, "Paths (from argument): ", paste(paths, collapse=", "));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Getting root paths
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying root (search) paths");
  verbose && cat(verbose, "Requested root paths:");
  verbose && print(verbose, paths);

  if (is.null(paths)) {
    verbose && enter(verbose, "Using default root paths");

    # Default root paths: annotationData/, annotationData,<tags>/
    rootPattern <- "^annotationData(|,.*)";

    ## verbose && cat(verbose, "Directory patterns: ", rootPattern);
    paths <- list.files(path=".", pattern=rootPattern);
    paths <- sort(paths);

    verbose && cat(verbose, "Identified root paths:");
    verbose && print(verbose, paths);

    verbose && enter(verbose, "Adding root paths of aroma.* packages as a backup");
    # Searching for aroma.* packages with install.packages() is
    # more generic, but also much slower. Maybe later. /HB 2011-03-03
    pkgs <- c("aroma.core", "aroma.affymetrix", "aroma.cn");
    pkgPaths <- sapply(pkgs, FUN=function(pkg) {
      system.file("annotationData", package=pkg);
    });
    pkgPaths <- pkgPaths[nchar(pkgPaths) > 0];
    verbose && cat(verbose, "Additional root paths:");
    verbose && print(verbose, pkgPaths);
    paths <- c(paths, pkgPaths);
    verbose && exit(verbose);


    # Nothing to do?
    if (length(paths) == 0) {
      verbose && cat(verbose, "None of the search paths exists. Skipping.");
      verbose && exit(verbose);
      verbose && exit(verbose);
      return(NULL);
    }

    verbose && exit(verbose);
  } else {
    # Split path strings by semicolons.
    paths <- unlist(strsplit(paths, split=";"));

    verbose && cat(verbose, "Parsed root paths:");
    verbose && print(verbose, paths);
  }

  # Expand any file system links
  paths <- sapply(paths, FUN=function(path) {
    Arguments$getReadablePath(path, mustExist=FALSE);
  });

  # Keep only existing search paths
  paths <- paths[sapply(paths, FUN=isDirectory)];
  verbose && cat(verbose, "Existing root paths:");
  verbose && print(verbose, paths);

  verbose && exit(verbose);

  # Nothing to do?
  if (length(paths) == 0) {
    verbose && cat(verbose, "None of the search paths exists. Skipping.");
    verbose && exit(verbose);
    return(NULL);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Getting subdir sets in root paths, i.e. <rootPath>/<set>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying all <rootPath>/<set>/ paths");
  verbose && cat(verbose, "Set: ", set);

  # Expand any file system links
  paths <- file.path(paths, set);
  paths <- sapply(paths, FUN=function(path) {
    Arguments$getReadablePath(path, mustExist=FALSE);
  });

  verbose && cat(verbose, "Possible <rootPath>/<set>/ paths:");
  verbose && print(verbose, paths);

  # Keep only existing search paths
  paths <- paths[sapply(paths, FUN=isDirectory)];
  verbose && cat(verbose, "Existing <rootPath>/<set>/ paths:");
  verbose && print(verbose, paths);

  verbose && exit(verbose);

  # Nothing to do?
  if (length(paths) == 0) {
    verbose && cat(verbose, "None of the search paths exists. Skipping.");
    verbose && exit(verbose);
    return(NULL);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Getting named sub sets, i.e. <rootPath>/<set>/<name>?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(name)) {
    verbose && enter(verbose, "Identifying all <rootPath>/<set>/<name>/ paths");
    verbose && cat(verbose, "Name: ", name);

    # Expand any file system links
    paths <- file.path(paths, name);
    paths <- sapply(paths, FUN=function(path) {
      Arguments$getReadablePath(path, mustExist=FALSE);
    });

    verbose && cat(verbose, "Possible <rootPath>/<set>/<name>/ paths:");
    verbose && print(verbose, paths);

    # Keep only existing search paths
    paths <- paths[sapply(paths, FUN=isDirectory)];
    verbose && cat(verbose, "Existing <rootPath>/<set>/<name>/ paths:");
    verbose && print(verbose, paths);

    verbose && exit(verbose);

    # Nothing to do?
    if (length(paths) == 0) {
      verbose && cat(verbose, "None of the search paths exists. Skipping.");
      verbose && exit(verbose);
      return(NULL);
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Now, search all paths
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && cat(verbose, "Final paths to search:");
  verbose && print(verbose, paths);

  # AD HOC: Searching in sibling root paths is only done in the correct
  # order when firstOnly == FALSE. /HB 2011-03-03
##  if (firstOnly) {
##    warning("Searching for annotation data in multiple root paths with 'firstOnly=TRUE' is *not* guaranteed to be done in the correct order. /HB 2011-03-03");
##  }

  # Search recursively for all files
  args <- list(...);
  args$pattern <- pattern;
  args$allFiles <- private;
  args$paths <- paths;
  args$recursive <- TRUE;
  args$firstOnly <- FALSE;
  verbose && cat(verbose, "Arguments to findFiles:");
  verbose && str(verbose, args);

  pathnames <- do.call(localFindFiles, args=args);

  verbose && cat(verbose, "All located pathname(s):");
  verbose && print(verbose, pathnames);

  verbose && enter(verbose, "Reorder by fullname");
  pathnames <- sortByFullName(pathnames, verbose=verbose);
  verbose && print(verbose, pathnames);
  verbose && exit(verbose);

  # Keep first match?
  if (firstOnly && length(pathnames) > 1) {
    pathnames <- pathnames[1];
    verbose && cat(verbose, "Returning the first one only");
    verbose && print(verbose, pathnames);
  }

  verbose && exit(verbose);

  pathnames;
}, protected=TRUE)  # findAnnotationData()

############################################################################
# HISTORY:
# 2012-08-29
# o ROBUSTNESS: Added argument 'escape' to findAnnotationData(), which
#   causes such symbols that exist in argument 'pattern' to be
#   automatically be escaped.  This was done in order for chip types
#   such as "Mapping250K_Nsp+Sty" to work by default.
# 2012-04-16
# o CLEANUP: findAnnotationData() no longer needs affxparser, because
#   findFiles() is in R.utils (>= 1.13.1).
# 2011-03-03
# o Now findAnnotationData() again guarantees that the returned pathnames
#   are ordered by the (length of the) fullnames.
# o GENERALIZATION: Now findAnnotationData() falls back to annotation data
#   available in any of the aroma.* packages.
# o GENERALIZATION: In addition to search <rootPath>/<set>/<name> paths,
#   findAnnotationData() can also search <rootPath>/<set>/ by not
#   specifying argument 'name' (or setting it to NULL).
# o SPEEDUP: Now findAnnotationData() returns NULL as soon as it knows
#   there are no root paths, subdir sets etc to search.
# 2011-02-19
# o GENERALIZATION: Extended the default root paths of findAnnotationData()
#   to be annotationData/ and annotationData,<tags>/
# 2009-02-10
# o Now findAnnotationData() always returns pathnames ordered by the length
#   of their fullnames. Before this was only done if 'firstOnly=TRUE'.
# 2008-05-21
# o Updated findAnnotationData() to only "import" affxparser.
# 2008-05-18
# o BUG FIX: findAnnotationDataByChipType(chipType="GenomeWideSNP_6",
#   pattern="^GenomeWideSNP_6.*[.]ugp$") would find file
#   'GenomeWideSNP_6,Full,na24.ugp' before 'GenomeWideSNP_6,na24.ugp'.
#   Now we return the one with the shortest full name.
# 2008-05-10
# o BUG FIX: When searching with 'firstOnly=FALSE', findAnnotationData()
#   was identifying files that are in "private" directory.  This is how
#   affxparser::findFiles() works.  Such files are now filtered out.
# 2008-05-09
# o Removed the option to specify the annotation data path by the option
#   'aroma.affymetrix.settings'.
# 2008-02-14
# o Added argument 'private=TRUE' to findAnnotationData().
# 2007-09-15
# o Created from findAnnotationDataByChipType.R.
############################################################################
