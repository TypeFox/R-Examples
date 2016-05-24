###########################################################################/**
# @RdocClass AromaMicroarrayDataSetTuple
#
# @title "The AromaMicroarrayDataSetTuple class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "R.filesets::GenericDataFileSetList".}
#   \item{.setClass}{The name of the class of the input set.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#
#*/###########################################################################
setConstructorS3("AromaMicroarrayDataSetTuple", function(..., .setClass="AromaMicroarrayDataSet") {
  extend(GenericDataFileSetList(..., .setClass=.setClass), "AromaMicroarrayDataSetTuple");
})


setMethodS3("as.AromaMicroarrayDataSetTuple", "AromaMicroarrayDataSetTuple", function(this, ...) {
  # Nothing to do
  this;
})


setMethodS3("as.character", "AromaMicroarrayDataSetTuple", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, paste("Name:", getName(this)));
  s <- c(s, paste("Tags:", paste(getTags(this), collapse=",")));
  s <- c(s, paste("Chip types:", paste(getChipTypes(this), collapse=", ")));
  dsList <- getSets(this);
  for (ds in dsList) {
    s <- c(s, as.character(ds));
  }
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));

  GenericSummary(s);
}, private=TRUE)



setMethodS3("indexOf", "AromaMicroarrayDataSetTuple", function(this, arrays=NULL, ...) {
  # Argument 'arrays':
  if (is.numeric(arrays)) {
    n <- length(this);
    arrays <- Arguments$getIndices(arrays, max=n);
  } else {
    arrays <- NextMethod("indexOf", arrays, onMissing="error");
  }

  arrays;
}, protected=TRUE)


setMethodS3("getAsteriskTags", "AromaMicroarrayDataSetTuple", function(this, ...) {
  "";
}, protected=TRUE)


setMethodS3("getTags", "AromaMicroarrayDataSetTuple", function(this, collapse=NULL, ...) {
  # Get tags of chip-effect set
  dsList <- getSets(this);

  # Get data set tags
  tags <- lapply(dsList, FUN=getTags);

  # Keep common tags
  tags <- getCommonListElements(tags);
  tags <- tags[[1]];
  tags <- unlist(tags, use.names=FALSE);

  # Add optional tuple tags
  tags <- c(tags, this$.tags);

  # In case this$.tags is not already split
  tags <- strsplit(tags, split=",", fixed=TRUE);
  tags <- unlist(tags);

  # Update asterisk tags
  tags[tags == "*"] <- getAsteriskTags(this, collapse=",");

  # Remove empty tags
  tags <- Arguments$getTags(tags, collapse=NULL);

  # Remove duplicated tags
  tags <- locallyUnique(tags);

  # Collapsed or split?
  tags <- Arguments$getTags(tags, collapse=collapse);

  tags;
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
#   @seemethod "getChipTypes".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("nbrOfChipTypes", "AromaMicroarrayDataSetTuple", function(this, ...) {
  length(getChipTypes(this, ...));
})





setMethodS3("getFullNames", "AromaMicroarrayDataSetTuple", function(this, arrays=NULL, exclude=NULL, translate=TRUE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  getFullNameOfList <- function(cfList, ...) {
    # Get sample name
    names <- sapply(cfList, FUN=getName);
    names <- names[!is.na(names)];
    # Sanity check
    stopifnot(length(names) > 0);

    name <- names[1];

    # Get chip-effect tags *common* across chip types
    tags <- lapply(cfList, FUN=getTags, ...);
    tags <- lapply(tags, FUN=function(x) {
      # To avoid warning on na.omit(NULL)
      if (length(x) > 0) na.omit(x) else x;
    });
    tags <- lapply(tags, FUN=setdiff, exclude);
    tags <- getCommonListElements(tags);
    tags <- tags[[1]];
    tags <- unlist(tags, use.names=FALSE);
    tags <- locallyUnique(tags);

    fullname <- paste(c(name, tags), collapse=",");

    fullname;
  } # getFullNameOfList()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'arrays':
  arrays <- indexOf(this, arrays);

  # Argument 'exclude':
  exclude <- Arguments$getCharacters(exclude);


  fullnames <- c();
  for (ii in arrays) {
    cfList <- getFileList(this, ii, ...);
    # Call local function
    fullname <- getFullNameOfList(cfList, translate=translate);
    fullnames <- c(fullnames, fullname);
  }

  fullnames;
})



setMethodS3("getChipTypes", "AromaMicroarrayDataSetTuple", function(this, fullname=FALSE, merge=FALSE, collapse="+", ...) {
  dsList <- getSets(this);
  chipTypes <- sapply(dsList, FUN=getChipType, fullname=fullname);

  # Invariant for order
#  chipTypes <- sort(chipTypes);

  # Merge to a single string?
  if (merge) {
    chipTypes <- mergeByCommonTails(chipTypes, collapse=collapse);
  }

  chipTypes;
})


setMethodS3("getSets", "AromaMicroarrayDataSetTuple", function(this, ...) {
  res <- NextMethod("getSets");
  # Name sets by their chip types
  chipTypes <- sapply(res, FUN=getChipType);
  chipTypes <- gsub(",monocell", "", chipTypes);
  names(res) <- chipTypes;
  res;
})



setMethodS3("byPath", "AromaMicroarrayDataSetTuple", abstract=TRUE, static=TRUE, protected=TRUE)



##############################################################################
# HISTORY:
# 2011-11-19
# o CLEANUP: Now getFullNames() for AromaMicroarrayDataSetTuple no
#   longer produces a warning on "is.na() applied to non-(list or
#   vector) of type 'NULL'".
# 2009-12-31
# o getSets() of AromaMicroarrayDataSetTuple overrides the default method
#   by adding the chip types as the names of the returns list.
# 2009-12-30
# o Renamed indexOfArrays() to indexOf().
# o Dropped nbrOfArrays(); use nbrOfFiles() instead.
# o CLEAN UP: Dropped get- and setAlias().
# o Made into a GenericDataFileSetList.
# o Create from AromaMicroarrayDataSetTuple.R.
# 2008-06-03
# o Added argument 'fullnames=FALSE' to getChipTypes() for
#   AromaMicroarrayDataSetTuple.
# 2008-05-21
# o BUG FIX: getFullNames() was not passing down new 'translate' correctly.
# 2008-05-16
# o Abstract methods right now: byPath().
# o Made getChipTypes() platform independent (not querying CDF:s anymore).
# o Now getFullName() getTableOfArrays(), getArrays(), indexOfArrays(),
#   getArrayTuple(), and asMatrixOfFiles() passes down '...' (for the
#   purpose passing down argument 'translate').
# o Now getFullNames() passes down 'translate'.
# 2008-03-29
# o getTableOfArrays() of AromaMicroarrayDataSetTuple returned the incorrect
#   array indices for the 2nd chip type if different arrays in the two sets.
# 2008-03-11
# o Renamed getTuple() to getArrayTuple().
# 2007-03-29
# o Added asMatrixOfFiles().
# 2007-03-20
# o Now getArrays() returns a named list where the names are the result from
#   getFullNames().
# 2007-03-19
# o TODO: Handle replicated sample names. It is not clear how this should be
#   done.
# o Created from GladModel.R.
##############################################################################
