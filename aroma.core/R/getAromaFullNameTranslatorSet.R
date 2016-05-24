# NOW:
# (1) annotationData/dataSets/<dataSet>/<chipType>/
# (2) annotationData/dataSets/<dataSet>/
#
# PROBABLY BETTER/MORE SPECIFIC (THINK OO DESIGN):
# (1) <rootPath>/<dataSet>/<chipType>/
# (2) <rootPath>/<dataSet>/
# (3) annotationData/dataSets/<dataSet>/<chipType>/
# (4) annotationData/dataSets/<dataSet>/
# (5) annotationData/dataSets/ [?!?]
setMethodS3("getAromaFullNameTranslatorSet", "character", function(dataSet, ..., chipType=NULL, paths=c("annotationData/dataSets/"), verbose=FALSE) {
  # Find existing dataset directory
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'dataSet':
  dataSet <- Arguments$getCharacter(dataSet);

  # Arguments 'chipType':
  if (!is.null(dataSet)) {
    chipType <- Arguments$getCharacter(chipType);
  }

  # Arguments 'paths':
  paths <- Arguments$getCharacters(paths);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Locating fullname TabularTextFileSet:s");
  verbose && cat(verbose, "Data set: ", dataSet);
  verbose && cat(verbose, "Chip type: ", chipType);

  subdirs <- c(chipType, ".");
  res <- getFullNameTranslatorSet(dataSet, ..., subdirs=subdirs, 
                                    paths=paths, verbose=less(verbose, 5));
  verbose && print(verbose, res);

  verbose && exit(verbose);

  res;
}, protected=TRUE) # getAromaFullNameTranslatorSet()


setMethodS3("getAromaFullNameTranslatorSet", "AromaMicroarrayDataSet", function(this, ...) {
  getAromaFullNameTranslatorSet(getFullName(this), tags=NULL, chipType=getChipType(this, fullname=FALSE), ...);
}, protected=TRUE)

setMethodS3("getAromaFullNameTranslatorSet", "AromaUnitSignalBinarySet", function(this, ...) {
  getAromaFullNameTranslatorSet(getFullName(this), tags=NULL, chipType=getChipType(this, fullname=FALSE), ...);
}, protected=TRUE)




setMethodS3("getFullNameTranslatorSet", "character", function(dataSet, ..., firstOnly=FALSE, verbose=FALSE) {
  # Find existing dataset directory
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'dataSet':
  dataSet <- Arguments$getCharacter(dataSet);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Locating fullname TabularTextFileSet:s");
  verbose && cat(verbose, "Data set: ", dataSet);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Find all existing search paths
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  paths <- GenericDataFileSet$findByName(dataSet, ..., firstOnly=firstOnly,
                                mustExist=FALSE, verbose=less(verbose, 5));

  verbose && cat(verbose, "Existing search paths:");
  verbose && print(verbose, paths);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Scan for existing fullname translator files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocate empty result
  res <- TabularTextFileSet();

  for (path in paths) {
    ds <- TabularTextFileSet$byPath(path, pattern=",fullnames[.]txt$");
    append(res, ds);
  } # for (path ...)

  verbose && print(verbose, res);

  verbose && exit(verbose);

  res;
}, protected=TRUE) # getFullNameTranslatorSet()




# NOT USED. /HB 2010-05-26
setMethodS3("getFullNameTranslatorSet", "GenericDataFileSet", function(this, subdirs=c("*"), paths=getParent(getPath(this)), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'dirs':
  subdirs <- Arguments$getCharacters(subdirs);

  # Arguments 'paths':
  paths <- Arguments$getCharacters(paths);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Locating fullname TabularTextFileSet:s");
  dataSet <- getFullName(this);
  verbose && cat(verbose, "Argument 'dataSet': ", dataSet);
  verbose && cat(verbose, "Argument 'subdirs': ");
  verbose && print(verbose, subdirs);
  verbose && cat(verbose, "Argument 'paths':");
  verbose && print(verbose, paths);

  res <- getFullNameTranslatorSet(dataSet, subdirs=subdirs, paths=paths,
                                           ..., verbose=less(verbose, 1));

  verbose && exit(verbose);

  res;
}, protected=TRUE)


#############################################################################
# HISTORY:
# 2010-05-26
# o Added auxillary/genric getFullNameTranslatorSet().  We should ideally 
#   define a class and corresponding findByName() and byName() for this.
# o Added auxillary getAromaFullNameTranslatorSet() utilizing
#   getFullNameTranslatorSet().
# o Created.
#############################################################################
