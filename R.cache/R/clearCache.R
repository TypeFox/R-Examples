#########################################################################/**
# @RdocDefault clearCache
#
# @title "Removes all files in a cache file directory"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{path}{A @character string specifying the directory to be cleared.
#     By default, the path is what is returned by @see "getCachePath"
#     which arguments \code{...} are also passed.}
#   \item{...}{Arguments passed to @see "getCachePath", especially
#     argument \code{dirs} to specify subdirectories.}
#   \item{recursive}{If @TRUE, subdirectories are also removed, otherwise
#     just the files in the specified directory.}
#   \item{prompt}{If @TRUE, the user will be prompted to confirm that
#     the directory will cleared before files are removed.}
# }
#
# \value{
#   Returns (invisibly) a @character @vector of pathnames of the
#   files removed.  If no files were removed, @NULL is returned.
# }
#
# @author
#
# \details{
#   If the specified directory does not exists, an exception is thrown.
# }
#
# @keyword "programming"
# @keyword "IO"
# @keyword "internal"
#*/#########################################################################
setMethodS3("clearCache", "default", function(path=getCachePath(...), ..., recursive=FALSE, prompt=TRUE && interactive()) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'recursive':
  recursive <- Arguments$getLogical(recursive);

  # Argument 'path':
  path <- Arguments$getReadablePath(path, mustExist=TRUE);

  # Argument 'prompt':
  prompt <- Arguments$getLogical(prompt);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Find files to be removed
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  allFiles <- listDirectory(path, allNames=TRUE, fullNames=TRUE, recursive=recursive);

  # Exclude '.' and '..' (just in case; listDirectory() shouldn't return them)
  excl <- grep("[.][.]*$", allFiles);
  if (length(excl) > 0L) allFiles <- allFiles[-excl];

  # Exclude 'README.txt'
  excl <- grep("README.txt$", allFiles);
  if (length(excl) > 0L) allFiles <- allFiles[-excl];

  nbrOfFiles <- length(allFiles);
  if (nbrOfFiles == 0L) {
    if (prompt) {
      message("Nothing to clear. Cache directory is empty: ", path, "\n", sep="");
    }
    return(invisible(NULL));
  }

  # Identify files and directories
  isdir <- file.info(allFiles)$isdir;
  dirs <- allFiles[isdir];
  files <- allFiles[!isdir];

  # Remove subdirectories before parent ones.
  o <- order(nchar(dirs), decreasing=TRUE);
  dirs <- dirs[o];

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Prompt user?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (prompt) {
    answer <- ".";
    while (!(answer %in% c("y", "n", ""))) {
      message(sprintf("Are you really sure you want to delete the %d files and %d directories in '%s'? [y/N]: ", length(files), length(dirs), path));
      answer <- tolower(readline());
    }
    if (answer != "y") {
      return(invisible(NULL));
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Remove files and directories
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (a) Remove files
  removed <- file.remove(files);
  filesR <- files[!removed];

  # (b) Remove subdirectories
  # Here we could use unlink(..., recursive=TRUE), but it is
  # SUPER-DUPER DANGEROUS to do that, because it may spawn off a
  # recursive deletion in a different place (in case there is a
  # bug or an inconsistency in list.files() etc), but also if we
  # forget to drop '.' and '..' from list.files(). /HB 2012-11-28
  removed <- sapply(dirs, FUN=function(dir) {
    filesT <- list.files(path=dir, all.files=TRUE);
    filesT <- setdiff(filesT, c(".", ".."));
    # Remove only empty directories
    if (length(filesT) > 0L) return(FALSE);
    removeDirectory(dir);
  });
  dirsR <- dirs[!removed];


  # Files and directories removed
  files <- sort(setdiff(files, filesR));
  dirs <- sort(setdiff(dirs, dirsR));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Report results?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (prompt) {
    msg <- sprintf("Cache cleared. Removed %d files and %d directories", length(files), length(dirs));
    if (length(filesR) + length(dirsR) > 0L) {
      msg <- sprintf("%s, but failed to remove another %d files and another %d directories", msg, length(filesR), length(dirsR));
    }
    message(sprintf("%s.\n", msg));
  }

  # Add a README.txt file, if missing.
  .addREADME();

  invisible(c(dirs, files));
})


############################################################################
# HISTORY:
# 2012-11-28
# o GENERALIZATION: Now clearCache(..., recursive=TRUE) removes all
#   cache files in subdirectories too.  The actual subdirectories are
#   not removed.
# 2012-11-27
# o BUG FIX: clearCache() would give error "object 'dirs' not found".
# 2011-05-19
# o Added Rdoc comments.
# 2011-04-02
# o BUG FIX: clearCache() would also report on subdirectories.
# 2005-12-09
# o BUG FIX: 'prompt=FALSE' would not clear cache.
# 2005-12-07
# o Created.
############################################################################
