###########################################################################/**
# @RdocClass ChromosomalModel
#
# @title "The ChromosomalModel class"
#
# \description{
#  @classhierarchy
#
#  This \emph{abstract} class represents a chromosomal model.
# }
#
# @synopsis
#
# \arguments{
#   \item{cesTuple}{A @see "AromaMicroarrayDataSetTuple".}
#   \item{tags}{A @character @vector of tags.}
#   \item{genome}{A @character string specifying what genome is process.}
#   \item{chromosomes}{(optional) A @vector specifying which chromosomes
#    to process.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Requirements}{
#   This class requires genome information annotation files for
#   every chip type.
# }
#
# @author
#*/###########################################################################
setConstructorS3("ChromosomalModel", function(cesTuple=NULL, tags="*", genome="Human", chromosomes=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'cesTuple':
  if (!is.null(cesTuple)) {
    # Coerce, if needed
    if (!inherits(cesTuple, "AromaMicroarrayDataSetTuple")) {
      cesTuple <- as.AromaMicroarrayDataSetTuple(cesTuple);
    }
  }

  # Argument 'tags':
  tags <- Arguments$getTags(tags, collapse=NULL);

  this <- extend(Object(), "ChromosomalModel",
    .alias = NULL,
    .cesTuple = cesTuple,
    .chromosomes = NULL,
    .tags = tags,
    .genome = genome
  );

  # Validate?
  if (!is.null(this$.cesTuple)) {
    # Assert that a genome annotation file exists
    gf <- getGenomeFile(this);
    this <- setChromosomes(this, chromosomes);
  }

  this;
}, abstract=TRUE)


setMethodS3("as.character", "ChromosomalModel", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, paste("Name:", getName(this)));
  s <- c(s, paste("Tags:", getTags(this, collapse=",")));
  s <- c(s, paste("Chip type (virtual):", getChipType(this)));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  tuple <- getSetTuple(this);
  chipTypes <- getChipTypes(tuple);
  nbrOfChipTypes <- length(chipTypes);
  s <- c(s, sprintf("Number of chip types: %d", nbrOfChipTypes));
  s <- c(s, sprintf("Chip types: %d", paste(chipTypes, collapse=", ")));

  s <- c(s, "List of data sets:");
  s <- c(s, as.character(tuple));

  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  GenericSummary(s);
}, protected=TRUE)


setMethodS3("getRootPath", "ChromosomalModel", function(this, ...) {
  tag <- getAsteriskTags(this)[1];
  sprintf("%sData", tolower(tag));
}, protected=TRUE)


setMethodS3("getParentPath", "ChromosomalModel", function(this, ...) {
  # Root path
  rootPath <- getRootPath(this);

  # Full name
  fullname <- getFullName(this);

  # The full path
  path <- filePath(rootPath, fullname);
  path <- Arguments$getWritablePath(path);

  path;
}, protected=TRUE)


setMethodS3("getPath", "ChromosomalModel", function(this, ...) {
  path <- getParentPath(this, ...);

  # Chip type
  chipType <- getChipType(this);

  # The full path
  path <- filePath(path, chipType);
  path <- Arguments$getWritablePath(path);

  path;
})

setMethodS3("getReportPath", "ChromosomalModel", function(this, ...) {
  rootPath <- "reports";

  # Data set name
  name <- getName(this);

  # Data set tags
  tags <- getTags(this, collapse=",");

  # Get chip type
  chipType <- getChipType(this);

  # Image set
  set <- getSetTag(this);

  # The report path
  path <- filePath(rootPath, name, tags, chipType, set);
  path <- Arguments$getWritablePath(path);

  path;
}, protected=TRUE)



setMethodS3("getSetTuple", "ChromosomalModel", function(this, ...) {
  this$.cesTuple;
}, protected=TRUE)



setMethodS3("getSets", "ChromosomalModel", function(this, ...) {
  tuple <- getSetTuple(this);
  getSets(tuple);
})


###########################################################################/**
# @RdocMethod nbrOfChipTypes
#
# @title "Gets the number of chip types"
#
# \description{
#  @get "title" used in the model.
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
setMethodS3("nbrOfChipTypes", "ChromosomalModel", function(this, ...) {
  tuple <- getSetTuple(this);
  nbrOfChipTypes(tuple, ...);
})



setMethodS3("getListOfUnitNamesFiles", "ChromosomalModel", function(this, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Retrieving unit names files");

  tuple <- getSetTuple(this);

  tryCatch({
    unfList <- getListOfUnitNamesFiles(tuple, ...);
  }, error = function(ex) {
    msg <- sprintf("Failed to located unit-names files for one of the chip types (%s). The error message was: %s", paste(getChipTypes(this), collapse=", "), ex$message);
    throw(msg);
  });

  verbose && exit(verbose);

  unfList;
}, private=TRUE)


setMethodS3("getListOfAromaUgpFiles", "ChromosomalModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Retrieving list of UGP files");

  tuple <- getSetTuple(this);
#  unfList <- getListOfUnitNamesFiles(this);

  ugpList <- NULL;
  tryCatch({
    verbose && enter(verbose, "Retrieving UGP files from unit names files");
#    ugpList <- lapply(unfList, FUN=getAromaUgpFile, verbose=less(verbose));
#   TODO: Why not do this?  /HB 2010-01-12
    ugpList <- lapply(tuple, FUN=getAromaUgpFile, verbose=less(verbose));
    verbose && exit(verbose);
  }, error = function(ex) {
    msg <- sprintf("Failed to located UGP files for one of the chip types (%s). Please note that DChip GenomeInformation files are no longer supported.  The error message was: %s", paste(getChipTypes(this), collapse=", "), ex$message);
    throw(msg);
  });

  verbose && exit(verbose);

  ugpList;
}, protected=TRUE)


setMethodS3("getListOfUnitTypesFiles", "ChromosomalModel", function(this, ...) {
  tuple <- getSetTuple(this);
  getListOfUnitTypesFiles(tuple, ...);
}, private=TRUE)



setMethodS3("getChipTypes", "ChromosomalModel", function(this, ...) {
  tuple <- getSetTuple(this);
  getChipTypes(tuple, ...);
})


###########################################################################/**
# @RdocMethod getChipType
#
# @title "Gets a label for all chip types merged"
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
#*/###########################################################################
setMethodS3("getChipType", "ChromosomalModel", function(this, ...) {
  getChipTypes(this, merge=TRUE, ...);
})



###########################################################################/**
# @RdocMethod getNames
#
# @title "Gets the names of the arrays"
#
# \description{
#  @get "title" available to the model.
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
setMethodS3("getNames", "ChromosomalModel", function(this, ...) {
  tuple <- getSetTuple(this);
  getNames(tuple, ...);
})


setMethodS3("getFullNames", "ChromosomalModel", function(this, ...) {
  tuple <- getSetTuple(this);
  getFullNames(tuple, ...);
})





###########################################################################/**
# @RdocMethod getTableOfArrays
#
# @title "Gets a table of arrays"
#
# \description{
#  @get "title" showing their availability across chip types.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a \eqn{NxK} @matrix of @integers where \eqn{N} is the total number
#  of arrays and \eqn{K} is the number of chip types in the model.  The row
#  names are the names of the arrays, and the column names are the chip types.
#  If data is available for array \eqn{n} and chip type \eqn{k}, cell
#  \eqn{(n,k)} has value \eqn{n}, otherwise @NA.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getTableOfArrays", "ChromosomalModel", function(this, ...) {
  tuple <- getSetTuple(this);
  getTableOfArrays(tuple, ...);
}, protected=TRUE, deprecated=TRUE)


setMethodS3("indexOf", "ChromosomalModel", function(this, patterns=NULL, ..., onMissing=c("error", "NA")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'onMissing':
  onMissing <- match.arg(onMissing);


  # If 'patterns' is numeric, then...
  if (is.numeric(patterns)) {
    names <- getNames(this);
    n <- length(names);
    res <- Arguments$getIndices(patterns, max=n);
    names(res) <- names[res];
    return(res);
  }

  # ...otherwise, reuse indexOf() for GenericDataFileSet in R.filesets.
  indexOf.GenericDataFileSet(this, patterns=patterns, ...,
                                                         onMissing=onMissing);
})




###########################################################################/**
# @RdocMethod nbrOfArrays
#
# @title "Gets the number of arrays"
#
# \description{
#  @get "title" used in the model.
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
setMethodS3("nbrOfArrays", "ChromosomalModel", function(this, ...) {
  length(getNames(this, ...));
})


setMethodS3("getName", "ChromosomalModel", function(this, collapse="+", ...) {
  name <- getAlias(this);

  if (is.null(name)) {
    tuple <- getSetTuple(this);
    name <- getName(tuple, ...);
  }

  name;
})


setMethodS3("getAsteriskTags", "ChromosomalModel", function(this, collapse=NULL, ...) {
  # Create a default asterisk tags for any class by extracting all
  # capital letters and pasting them together, e.g. AbcDefGhi => ADG.
  name <- class(this)[1];

  # Remove any 'Model' suffixes
  name <- gsub("Model$", "", name);

  name <- capitalize(name);

  # Vectorize
  name <- strsplit(name, split="")[[1]];

  # Identify upper case
  name <- name[(toupper(name) == name)];

  # Paste
  name <- paste(name, collapse="");

  tag <- name;
}, protected=TRUE)




setMethodS3("getTags", "ChromosomalModel", function(this, collapse=NULL, ...) {
  tuple <- getSetTuple(this);
  tags <- getTags(tuple, collapse=collapse, ...);

  # Add model tags
  tags <- c(tags, this$.tags);

  # In case this$.tags is not already split
  tags <- strsplit(tags, split=",", fixed=TRUE);
  tags <- unlist(tags);

  # Update default tags
  asteriskTags <- getAsteriskTags(this, collapse=",");
  if (length(asteriskTags) == 0)
    asteriskTags <- "";
  tags[tags == "*"] <- asteriskTags;

  tags <- Arguments$getTags(tags, collapse=NULL);

  # Get unique tags
  tags <- locallyUnique(tags);

  # Collapsed or split?
  tags <- Arguments$getTags(tags, collapse=collapse);

  tags;
})


setMethodS3("getFullName", "ChromosomalModel", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})



###########################################################################/**
# @RdocMethod getChromosomes
# @alias setChromosomes.ChromosomalModel
#
# @title "Gets the chromosomes to be processed"
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
#*/###########################################################################
setMethodS3("getChromosomes", "ChromosomalModel", function(this, ...) {
  chromosomes <- this$.chromosomes;
  if (!is.null(chromosomes)) {
    return(chromosomes);
  }

  # The default is to process all available chromosomes
  ugpList <- getListOfAromaUgpFiles(this);
  chromosomes <- lapply(ugpList, FUN=getChromosomes);
  chromosomes <- unlist(chromosomes, use.names=TRUE);
  chromosomes <- sort(unique(chromosomes));

  chromosomes;
})

setMethodS3("setChromosomes", "ChromosomalModel", function(this, chromosomes=NULL, ...) {
  # Argument 'chromosomes':
  if (!is.null(chromosomes)) {
    chromosomes <- Arguments$getVector(chromosomes);
    chromosomes <- sort(unique(chromosomes));

    # All available chromosomes
    ugpList <- getListOfAromaUgpFiles(this);
    chromosomesA <- lapply(ugpList, FUN=getChromosomes);
    chromosomesA <- unlist(chromosomesA, use.names=TRUE);
    chromosomesA <- sort(unique(chromosomesA));

    unknown <- setdiff(chromosomes, chromosomesA);
    if (length(unknown) > 0L) {
      throw(sprintf("Unknown chromosomes detected: %s [%d]", hpaste(unknown), length(unknown)));
    }
  }

  this$.chromosomes <- chromosomes;

  invisible(this);
})






setMethodS3("getGenome", "ChromosomalModel", function(this, ...) {
  this$.genome;
})


setMethodS3("getGenomeFile", "ChromosomalModel", function(...) {
  getAromaGenomeTextFile(...);
}, protected=TRUE)


setMethodS3("getAromaGenomeTextFile", "ChromosomalModel", function(this, genome=getGenome(this), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'genome':
  genome <- Arguments$getCharacter(genome);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Locating genome annotation file");
  verbose && cat(verbose, "Genome name: ", genome);

  gf <- AromaGenomeTextFile$byGenome(genome, ..., verbose=verbose);
  verbose && print(verbose, gf);
  verbose && exit(verbose);

  gf;
}, protected=TRUE)  # getAromaGenomeTextFile()


setMethodS3("setGenome", "ChromosomalModel", function(this, genome, tags=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'genome':
  genome <- Arguments$getCharacter(genome, length=c(1,1));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  oldGenome <- this$.genome;

  fullname <- paste(c(genome, tags), collapse=",");
  verbose && cat(verbose, "Fullname: ", fullname);

  # Verify that there is an existing genome file
  tryCatch({
    this$.genome <- fullname;
    gf <- getGenomeFile(this, verbose=less(verbose, 10));
  }, error = function(ex) {
    this$.genome <- oldGenome;
    throw(ex$message);
  })

  invisible(oldGenome);
})



setMethodS3("getGenomeData", "ChromosomalModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Reading genome chromosome annotation file");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get genome annotation data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving genome annotation file");
  gf <- getGenomeFile(this, verbose=less(verbose, 10));
  verbose && exit(verbose);

  verbose && enter(verbose, "Reading data file");
  pathname <- getPathname(gf);
  verbose && cat(verbose, "Pathname: ", pathname);
  data <- readTable(pathname, header=TRUE,
                            colClasses=c(nbrOfBases="integer"), row.names=1);
  verbose && exit(verbose);

  verbose && enter(verbose, "Translating chromosome names");
  chromosomes <- row.names(data);
  map <- c("X"=23, "Y"=24, "Z"=25);
  for (kk in seq_along(map)) {
    chromosomes <- gsub(names(map)[kk], map[kk], chromosomes, fixed=TRUE);
  }
  row.names(data) <- chromosomes;
  verbose && exit(verbose);

  verbose && exit(verbose);

  data;
}, protected=TRUE)


setMethodS3("fit", "ChromosomalModel", abstract=TRUE);


setMethodS3("getSetTag", "ChromosomalModel", function(this, ...) {
  tolower(getAsteriskTags(this)[1]);
}, private=TRUE)


setMethodS3("getOutputSet", "ChromosomalModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Retrieving output set");

  verbose && enter(verbose, "Scanning output path");
  # Locate all
  path <- getPath(this);
  verbose && cat(verbose, "Path: ", path);
  fs <- GenericDataFileSet$byPath(path, ...);
  verbose && cat(verbose, "Number of matching files located: ", length(fs));
  verbose && exit(verbose);

  verbose && enter(verbose, "Keep those with fullnames matching the input data set");
  fullnames <- getFullNames(fs);
  verbose && cat(verbose, "Full names of *all* files found:");
  verbose && str(verbose, fullnames);

  # Drop extranous files
  keepFullnames <- getFullNames(this);
  verbose && cat(verbose, "Full names to be kept:");
  verbose && str(verbose, keepFullnames);

  patterns <- sprintf("^%s", fullnames);
  keep <- rep(FALSE, times=length(fullnames));
  for (pattern in patterns) {
    keep <- keep | (regexpr(pattern, fullnames) != -1);
  }

  if (any(!keep)) {
    verbose && enter(verbose, "Extract subset of files");
    keep <- which(keep);
    verbose && cat(verbose, "Keeping indices:");
    verbose && str(verbose, keep);
    fs <- extract(fs, keep);
    verbose && exit(verbose);
  }

  verbose && exit(verbose);

  verbose && print(verbose, fs);

  verbose && exit(verbose);

  fs;
}, private=TRUE)


setMethodS3("getAlias", "ChromosomalModel", function(this, ...) {
  this$.alias;
}, protected=TRUE)


setMethodS3("getArrays", "ChromosomalModel", function(this, ...) {
  getNames(this, ...);
}, protected=TRUE, deprecated=TRUE)



##############################################################################
# HISTORY:
# 2013-10-03
# o Added argument 'chromosomes' to ChromosomalModel() and setChromosomes()
#   for the same class.  If not specified, the default is as before to infer
#   the set of chromosomes from the UGP files.
# 2011-03-03
# o Now getAromaGenomeTextFile() for ChromosomalModel utilizes byGenome()
#   for AromaGenomeTextFile to locate and return the AromaGenomeTextFile.
# o getGenomeFile() calls getAromaGenomeTextFile().
# 2011-02-28
# o UNDO: getArrays() was needed.
# 2011-02-19
# o CLEANUP: Removed getChipEffectFiles(), getListOfChipEffectSets(),
#   getArrays(), getAlias(), and setAlias() for ChromosomalModel.
# 2010-07-06
# o BUG FIX: indexOf() for ChromosomalModel would return NA if a search
#   pattern contained parenthesis '(' and ')'.  There was a similar issue
#   in indexOf() for GenericDataFileSet/List in R.filesets, which was
#   solved in R.filesets 0.8.3.  Now indexOf() for ChromosomalModel
#   utilizes ditto for GenericDataFileSet for its solution.
# 2010-03-02
# BUG FIX: Forgot argument 'verbose' of getOutputSet() of ChromosomalModel.
# 2010-02-19
# o Updated getGenomeFile() for ChromosomalModel such that it can be used
#   to locate other types of genome annotation files as well, files that
#   may be optional (without giving an error).
# 2010-02-18
# o Added getOutputSet() for ChromosomalModel.
# 2010-01-13
# o getListOfAromaUgpFiles() for ChromosomalModel no longer goes via
#   getListOfUnitNamesFiles().  This opens up the possibility to work with
#   data files without unit names files, e.g. smoothed CN data.
# 2009-11-18
# o CLEAN UP: Removed all Affymetrix specific classes/methods.
# 2009-11-16
# o CLEAN UP: The ChromosomalModel no longer checks 'combineAlleles'.
# o Now getChromosomes() of ChromosomalModel locates UGP files.
#   DChip GenomeInformation files are no longer supported for this.
# 2009-07-08
# o Added getListOfUnitTypesFiles() for ChromosomalModel.
# 2009-01-26
# o Removed get[]ListOfCdfs() from ChromosomalModel.
# o Removed deprectated get[]ListOfChipEffects() from ChromosomalModel.
# o Added getListOfAromaUgpFiles() to ChromosomalModel.
# o Added getListOfUnitNamesFiles() to ChromosomalModel.
# 2007-09-25
# o Extracted ChromosomalModel from CopyNumberSegmentationModel.  For
#   previous HISTORY, see that class.
##############################################################################
