# @author "HB"
setConstructorS3("AffymetrixTabularFile", function(...) {
  extend(TabularTextFile(...), c("AffymetrixTabularFile",
              uses("AromaPlatformInterface"), uses("FileCacheKeyInterface"))
  );
})


setMethodS3("translateColumnNames", "AffymetrixTabularFile", function(this, names, ...) {
  # Convert 'FOO_BAR_DOO' and 'FOO.BAR.DOO' to 'foo bar doo'?
  if (any(regexpr("[_.]", names) != -1)) {
    names <- tolower(gsub("[_.]", " ", names));
  }

  # Finally, convert 'Foo bar Doo' to 'fooBarDoo'
  names <- toCamelCase(names);

  names;
}, protected=TRUE)


setMethodS3("findByChipType", "AffymetrixTabularFile", function(static, chipType, tags=NULL, pattern=NULL, ...) {
  if (is.null(pattern)) {
    name <- paste(c(chipType, tags), collapse=",");
    pattern <- sprintf("^%s.*[.]...$", name);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search in annotationData/chipTypes/<chipType>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pathname <- findAnnotationDataByChipType(chipType, pattern, ...);
  pathname;
}, static=TRUE, protected=TRUE)



setMethodS3("byChipType", "AffymetrixTabularFile", function(static, chipType, tags=NULL, ...) {
  # Search for the genome information file
  pathname <- findByChipType(static, chipType, tags=tags, ...);
  if (is.null(pathname))
    throw("Failed to located Affymetrix tabular file: ", chipType);
  newInstance(static, pathname, ...);
}, static=TRUE)



############################################################################
# HISTORY:
# 2008-05-17
# o Now inherits from TabularTextFile.
# 2008-04-25
# o Now byChipType() passes '...' to the constructor, e.g. 'verbose'.
# 2007-09-16
# o Now AffymetrixTabularFile reports the translated/cleaned up column
#   names.  Column names on the file can be retrieved from getHeader().
# o Now colnames() returns column names in camelCase.
# 2007-09-14
# o Now inheriting from GenericTabularFile.
# 2007-09-10
# o Created from AffymetrixCsvGenomeInformation.R.
############################################################################
