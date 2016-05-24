###########################################################################/**
# @RdocClass Explorer
#
# @title "The Explorer class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{tags}{A @character @vector of tags to be added to the output path.}
#   \item{version}{An optional @character string.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Output directory structure}{
#   The \emph{main directory} of an Explorer report is
#    \code{reports/<name>/<subname>/}.
#   The \code{<name>} is typically the same as the name of the input
#   data set, and the \code{<subname>} is typically the tags of ditto.
#   This main directory is where main HTML document is stored.
#
#   For each chip type, real or "virtual" (combined), there is a
#   subdirectory with the same name as the chip type, i.e.
#    \code{reports/<name>/<subname>/<chiptype>/}.
#
#   For each chip type directory, there are set of subdirectories each
#   specifying a so called \emph{image layer}, e.g. an image layer
#   showing the raw data, another containing the estimates of a model
#   fit and so on.  Path format:
#    \code{reports/<name>/<subname>/<chiptype>/<image layer>/}.
#   In this directory all image files are stored, e.g. PNG files.
#
#   In some cases one do not want to all input tags to become part of the
#   subname, but instead for instance use those to name the image layer(s).
#   In such cases one has to override the default names.
# }
#
# @author
#
#*/###########################################################################
setConstructorS3("Explorer", function(tags="*", version="0", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'tags':
  tags <- Arguments$getTags(tags, collapse=NULL);

  # Argument 'version':
  version <- Arguments$getCharacter(version);

  extend(Object(), "Explorer",
    .version = version,
    .alias = NULL,
    .tags = tags,
    .arrays = NULL,
    .parallelSafe = FALSE
  )
})



setMethodS3("as.character", "Explorer", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, sprintf("Version: %s", getVersion(this)));
  s <- c(s, sprintf("Name: %s", getName(this)));
  s <- c(s, sprintf("Tags: %s", getTags(this, collapse=",")));
  s <- c(s, sprintf("Main path: %s", getMainPath(this)));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  GenericSummary(s);
}, protected=TRUE)


setMethodS3("getVersion", "Explorer", function(this, ...) {
  this$.version;
})


setMethodS3("getArraysOfInput", "Explorer", abstract=TRUE, protected=TRUE);


###########################################################################/**
# @RdocMethod getNames
#
# @title "Gets the names of the input samples"
#
# \description{
#  @get "title" for which the explorer is displaying results.
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
#*/###########################################################################
setMethodS3("getNames", "Explorer", function(this, ...) {
  names <- this$.arrays;

  if (is.null(names)) {
    names <- getArraysOfInput(this);
  }

  # Sanity check
  names <- Arguments$getCharacters(names);

  names;
})


###########################################################################/**
# @RdocMethod setArrays
#
# @title "Sets the arrays"
#
# \description{
#  @get "title" to be processed by the explorer.
# }
#
# @synopsis
#
# \arguments{
#   \item{arrays}{A @character (or @integer) @vector of arrays to be
#      considered. If @NULL, all arrays of the data set are considered.}
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
#*/###########################################################################
setMethodS3("setArrays", "Explorer", abstract=TRUE);



###########################################################################/**
# @RdocMethod nbrOfArrays
#
# @title "Gets the total number of arrays"
#
# \description{
#  @get "title" considered by the explorer.
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
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("nbrOfArrays", "Explorer", function(this, ...) {
  length(getNames(this));
})



###########################################################################/**
# @RdocMethod getAlias
#
# @title "Gets the alias of the output set"
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
#   Returns a @character, or @NULL if no alias is set.
# }
#
# @author
#
# \seealso{
#   @seemethod "setAlias".
#   @seemethod "getName".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getAlias", "Explorer", function(this, ...) {
  this$.alias;
}, protected=TRUE)



###########################################################################/**
# @RdocMethod setAlias
#
# @title "Sets the alias of the output set"
#
# \description{
#   @get "title".
#   If specified, the alias overrides the data set name, which is used by
#   default.
# }
#
# @synopsis
#
# \arguments{
#  \item{alias}{A @character string for the new alias of the output set.
#   The alias must consists of valid filename characters, and must not
#   contain commas, which are used to separate tags.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns itself invisibly.
# }
#
# @author
#
# \seealso{
#   @seemethod "getAlias".
#   @seemethod "getName".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("setAlias", "Explorer", function(this, alias=NULL, ...) {
  # Argument 'alias':
  if (!is.null(alias)) {
    alias <- Arguments$getFilename(alias);  # Valid filename?

    # Assert that no commas are used.
    if (regexpr("[,]", alias) != -1) {
      throw("Output-set aliases (names) must not contain commas: ", alias);
    }
  }

  this$.alias <- alias;

  invisible(this);
}, protected=TRUE)


###########################################################################/**
# @RdocMethod getName
#
# @title "Gets the name of the explorer"
#
# \description{
#  @get "title", which is the same as the name of the data set.
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
# \details{
#  If a name alias has not been set explicitly, the name of the data set will
#  used.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getName", "Explorer", function(this, ...) {
  name <- getAlias(this);
  if (is.null(name)) {
    name <- getNameOfInput(this);
  }
  name;
})



###########################################################################/**
# @RdocMethod getTags
#
# @title "Gets the tags of the explorer"
#
# \description{
#  @get "title", which are the tags of the data set plus additional tags.
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
#*/###########################################################################
setMethodS3("getTags", "Explorer", function(this, collapse=NULL, ...) {
  tags <- getTagsOfInput(this, ...);

  tags <- c(tags, this$.tags);

  # In case this$.tags is not already split
  tags <- strsplit(tags, split=",", fixed=TRUE);
  tags <- unlist(tags);

  tags <- locallyUnique(tags);

  # Update asterisk tags
  tags[tags == "*"] <- getAsteriskTags(this, collapse=",");

  tags <- Arguments$getTags(tags, collapse=NULL);

  tags <- locallyUnique(tags);

  # Collapsed or split?
  tags <- Arguments$getTags(tags, collapse=collapse);

  tags;
})

setMethodS3("getAsteriskTags", "Explorer", function(this, ...) {
  "";
}, protected=TRUE)


setMethodS3("getTagsOfInput", "Explorer", function(this, ...) {
  "";
}, protected=TRUE)

setMethodS3("getNameOfInput", "Explorer", abstract=TRUE, protected=TRUE);



setMethodS3("getFullName", "Explorer", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})



# tags <- "100K,CEU,testSet,ACC,-X,+300,RMA,A+B,w,FLN,SRMA,gauss,b=50000"
# Example: setReportPathPattern(ce, "^(.*),(SRMA,.*)(,CNC|)$");

setMethodS3("setReportPathPattern", "Explorer", function(this, pattern, ...) {
  # Argument 'pattern':
  pattern <- Arguments$getRegularExpression(pattern);
  this$.reportPathPattern <- pattern;
}, protected=TRUE)


setMethodS3("getReportPathPattern", "Explorer", function(this, ...) {
	this$.reportPathPattern;
}, protected=TRUE)

setMethodS3("splitByReportPathPattern", "Explorer", function(this, tags, ...) {
  # Argument 'tags':
  tags <- Arguments$getTags(tags, collapse=",");

  # Get subname and sampleLayerPrefix
	pattern <- getReportPathPattern(this);
  res <- list();
  if (is.null(pattern) || regexpr(pattern, tags) == -1) {
    res$subname <- tags;
  } else {
    res$subname <- gsub(pattern, "\\1", tags);
    res$sampleLayerPrefix <- gsub(pattern, "\\2", tags);
  }
  res;
}, protected=TRUE)


###########################################################################/**
# @RdocMethod getRootPath
#
# @title "Gets the root path of the output directory"
#
# \description{
#  @get "title" that is returned by @seemethod "getPath".
#  A root path is a directory in the current working directory.
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
#   @seemethod "getPath".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getRootPath", "Explorer", function(this, ...) {
  "reports";
})


setMethodS3("setSubname", "Explorer", function(this, value, ...) {
  oldValue <- this$.subname;
  this$.subname <- value;
  invisible(oldValue);
}, protected=TRUE)


setMethodS3("getSubname", "Explorer", function(this, ...) {
  # Preset?
  subname <- this$.subname;
  if (!is.null(subname))
    return(subname);

  # Infer from tags
  tags <- getTags(this, collapse=",");
  if (length(tags) == 0 || nchar(tags) == 0) {
    tags <- "raw";  # Default
  }

  subname <- splitByReportPathPattern(this, tags)$subname;
  if (is.null(subname))
    throw("ERROR: No subname could be inferred from tags: ", tags);

  subname;
}, protected=TRUE)


setMethodS3("getSampleLayerPrefix", "Explorer", function(this, ...) {
  # Infer from tags
  tags <- getTags(this, collapse=",");
  if (length(tags) == 0 || nchar(tags) == 0) {
    tags <- "raw";  # Default
  }
  prefix <- splitByReportPathPattern(this, tags)$sampleLayerPrefix;
  prefix;
}, protected=TRUE)



setMethodS3("getMainPath", "Explorer", function(this, ...) {
  # Create the (sub-)directory tree for the data set

  # Root path
  rootPath <- getRootPath(this);

  # Name
  name <- getName(this);

  # Subname
  subname <- getSubname(this);

  # The full path
  path <- filePath(rootPath, name, subname);
  path <- Arguments$getWritablePath(path);

  path;
}, protected=TRUE)


###########################################################################/**
# @RdocMethod getPath
#
# @title "Gets the path of the output directory"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{create}{If @TRUE, the path is created, otherwise not.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# \details{
#   Windows Shortcut links are recognized.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getPath", "Explorer", abstract=TRUE);


setMethodS3("getTemplatePath", "Explorer", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Locating template files for ChromosomeExplorer");
  # Search for template files
  rootPath <- getRootPath(this);
  path <- filePath(rootPath, "templates");
  path <- Arguments$getReadablePath(path, mustExist=FALSE);
  if (!isDirectory(path)) {
    path <- system.file("reports", "templates", package="aroma.core");
  }
  verbose && exit(verbose);

  path;
}, protected=TRUE)


setMethodS3("getIncludePath", "Explorer", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Locating include files for ChromosomeExplorer");
  # Search for include files
  path <- system.file("reports", "includes", package="aroma.core");
  verbose && exit(verbose);

  path;
}, protected=TRUE)


setMethodS3("addIncludes", "Explorer", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Setting up ", class(this)[1], " report files");

  destPath <- filePath(getRootPath(this), "includes");
  verbose && enter(verbose, "Copying template files");
  srcPath <- getIncludePath(this);
  verbose && cat(verbose, "Source path: ", srcPath);
  verbose && cat(verbose, "Destination path: ", destPath);

  pathnames <- copyDirectory(from=srcPath, to=destPath, copy.mode=FALSE,
                             recursive=TRUE, overwrite=force);
  verbose && exit(verbose);

  verbose && exit(verbose);
}, protected=TRUE)




setMethodS3("addIndexFile", "Explorer", function(this, filename=sprintf("%s.html", class(this)[1]), ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  srcPath <- getTemplatePath(this);
  srcPathname <- filePath(srcPath, "html", class(this)[1], filename);
  outPathname <- filePath(getMainPath(this), filename);

  if (force || !isFile(outPathname)) {
    verbose && enter(verbose, "Copying ", filename);
    verbose && cat(verbose, "Source pathname: ", srcPathname);
    verbose && cat(verbose, "Destination pathname: ", outPathname);
    if (!isFile(srcPathname))
      throw("File not found: ", srcPathname);

    copyFile(srcPathname, outPathname, overwrite=TRUE, copy.mode=FALSE);

    verbose && exit(verbose);
  }
}, protected=TRUE)



setMethodS3("updateSetupExplorerFile", "Explorer", function(this, data, ..., verbose=FALSE) {
  pkg <- "R.rsp";
  require(pkg, character.only=TRUE) || throw("Package not loaded: ", pkg);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'data':
  data <- Arguments$getInstanceOf(data, "environment");

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  outFile <- "setupExplorer.js";

  verbose && enter(verbose, "Updating ", outFile);

  # Get RSP-embedded source file
  mainPath <- getMainPath(this);
  setTuple <- getSetTuple(this);
  filename <- sprintf("%s.rsp", outFile);
  srcPath <- getTemplatePath(this);
  pathname <- filePath(srcPath, "rsp", class(this)[1], filename);
  verbose && cat(verbose, "Source: ", pathname);

  # Output destination
  outPath <- mainPath;
  verbose && cat(verbose, "Output path: ", outPath);
  outPath <- Arguments$getWritablePath(outPath);

  # Input data
  verbose && cat(verbose, "Input data:");
  verbose && str(verbose, as.list(data));

  verbose && enter(verbose, "Compiling RSP");
  js <- rfile(pathname, workdir=outPath, envir=data, postprocess=FALSE);
  verbose && exit(verbose);

  verbose && exit(verbose);

  invisible(js);
}, protected=TRUE) # updateSetupExplorerFile()



setMethodS3("setup", "Explorer", function(this, ..., force=FALSE) {
  # Setup includes/
  addIncludes(this, ..., force=force);

  # Setup HTML explorer page
  addIndexFile(this, ..., force=force);

  # Update Javascript files
  updateSetupExplorerFile(this, ...);
}, protected=TRUE)



###########################################################################/**
# @RdocMethod process
#
# @title "Generates image files, scripts and dynamic pages for the explorer"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{arrays}{A @vector of arrays specifying which arrays to
#    be considered.  If @NULL, all are processed.}
#   \item{...}{Not used.}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
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
#*/###########################################################################
setMethodS3("process", "Explorer", abstract=TRUE);


###########################################################################/**
# @RdocMethod display
#
# @title "Displays the explorer in the default browser"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
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
#*/###########################################################################
setMethodS3("display", "Explorer", function(this, filename=sprintf("%s.html", class(this)[1]), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Opening ", class(this)[1]);

  # The path to the explorer HTML document
  path <- getMainPath(this);
  pathname <- Arguments$getReadablePathname(filename, path=path, sbsolute=TRUE, mustExist=FALSE);

  # Just in case, is setup needed?
  if (!isFile(pathname)) {
    setup(this, verbose=less(verbose));
    if (!isFile(pathname))
      throw("Cannot open ", class(this)[1], ". No such file: ", pathname);
  }

  verbose && cat(verbose, "Pathname: ", pathname);

  # WORKAROUND: browseURL('foo/bar.html', browser=NULL), which in turn
  # calls shell.exec('foo/bar.html'), does not work on Windows, because
  # the OS expects backslashes.  [Should shell.exec() convert to
  # backslashes?]  By temporarily setting the working directory to that
  # of the file, this works around this issue.
  # Borrowed from R.rsp. /HB 2014-09-19
  if (isFile(pathname)) {
    path <- dirname(pathname);
    pathname <- basename(pathname);
    opwd <- getwd();
    on.exit(setwd(opwd));
    setwd(path);
  }

  res <- browseURL(pathname, ...);

  verbose && exit(verbose);

  invisible(res);
})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# DEPRECATED
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("getArrays", "Explorer", function(this, ...) {
  getNames(this, ...);
}, protected=TRUE, deprecated=TRUE)




##############################################################################
# HISTORY:
# 2014-01-17
# o CLEANUP: Now the Explorer class utilizes R.rsp::rfile() for compiling
#   RSP files instead of the to-be-deprecated rspToHtml().
# 2012-02-06
# o Added implementation of setup() to Explorer.
# o Added updateSetupExplorerFile() to Explorer.
# o Added getVersion() to Explorer.
# 2009-05-17
# o Added missing abstract method getArraysOfInput() to Explorer.
# o Moved the Explorer class and its support files under inst/ to aroma.core.
# 2008-06-05
# o Made getMainPath(), addIncludes(), addIndexFile() parallel safe.
# o Added getParallelSafe() and setParallelSafe().
# 2007-11-20
# o Now addIncludes() no longer passes '...' to copyDirectory().
# 2007-10-11
# o Now addIncludes() always copies missing files in the includes/ directory.
# 2007-03-24
# o BUG FIX: getPath() created the root path before trying to expand
#   Windows shortcuts.
# 2007-03-19
# o Created from ChromosomeExplorer.R.
##############################################################################
