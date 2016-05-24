###########################################################################/**
# @RdocClass AromaGenomeTextFile
#
# @title "The AromaGenomeTextFile class"
#
# \description{
#  @classhierarchy
#
#  An AromaGenomeTextFile represents a annotation tabular text file that
#  specifies the number of bases (nucleotides) per chromosome for a
#  particular genome/organism.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "R.filesets::TabularTextFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \details{
#   An AromaGenomeTextFile is a tab-delimited text file with a header
#   containing (at least) column names 'chromosome' and 'nbrOfBases'.
#   The 'chromosome' column specifies the chromosomes (character strings)
#   and the 'nbrOfBases' column specifies the lengths (integer) of the
#   chromosomes in number of bases (nucleotides).
#
#   The filename of an AromaGenomeTextFile should have format
#   "<genome>,chromosomes(,<tag>)*.txt",
#   and be located in annotationData/genomes/<genome>/, e.g.
#   annotationData/genomes/Human/Human,chromosomes,max,20090503.txt.
# }
#
# @examples "../incl/AromaGenomeTextFile.Rex"
#
# @author
#*/###########################################################################
setConstructorS3("AromaGenomeTextFile", function(...) {
  extend(TabularTextFile(...), c("AromaGenomeTextFile",
                                              uses("FileCacheKeyInterface")));
})


setMethodS3("readDataFrame", "AromaGenomeTextFile", function(this, ..., colClasses=c("*"="NULL", chromosome="character", nbrOfBases="integer", nbrOfGenes="integer")) {
  NextMethod("readDataFrame", colClasses=colClasses);
})



## setMethodS3("readGenomeData", "AromaGenomeTextFile", function(this, ..., verbose=FALSE) {
##   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##   # Validate arguments
##   # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
##   # Argument 'verbose':
##   verbose <- Arguments$getVerbose(verbose);
##   if (verbose) {
##     pushState(verbose);
##     on.exit(popState(verbose));
##   }
##
##
##   verbose && enter(verbose, "Reading genome chromosome annotation file");
##
##   data <- readDataFrame(this, ..., verbose=less(verbose, 10));
##
##   verbose && enter(verbose, "Translating chromosome names");
##   chromosomes <- row.names(data);
##   map <- c("X"=23, "Y"=24, "Z"=25);  # AD HOC; only for the human genome
##   for (kk in seq_along(map)) {
##     chromosomes <- gsub(names(map)[kk], map[kk], chromosomes, fixed=TRUE);
##   }
##   row.names(data) <- chromosomes;
##   verbose && exit(verbose);
##
##   verbose && exit(verbose);
##
##   data;
## })



setMethodS3("findByGenome", "AromaGenomeTextFile", function(static, genome, tags=NULL, pattern=sprintf("^%s,chromosomes(|,.*)*[.]txt$", genome, paste(tags, collapse=",")), paths=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'genome':
  genome <- Arguments$getCharacter(genome);

  # Argument 'tags':
  tags <- Arguments$getTags(tags, collapse=",");

  # Argument 'pattern':
  if (!is.null(pattern)) {
    pattern <- Arguments$getRegularExpression(pattern);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Locating genome annotation file");

  fullname <- paste(c(genome, tags), collapse=",");

  verbose && cat(verbose, "Genome name: ", genome);
  verbose && cat(verbose, "Genome tags: ", tags);
  verbose && cat(verbose, "Genome fullname: ", fullname);
  verbose && cat(verbose, "Pattern (fmtstr): ", pattern);
  patternT <- sprintf(pattern, fullname);
  verbose && cat(verbose, "Pattern (updated): ", patternT);
  verbose && cat(verbose, "Paths:");
  verbose && print(verbose, paths);

  pathname <- findAnnotationData(name=fullname, set="genomes",
                                 pattern=patternT, paths=paths, ...,
                                      verbose=less(verbose, 10));

  if (!is.null(pathname)) {
    verbose && cat(verbose, "Found file: ", pathname);
  }

  verbose && exit(verbose);

  pathname;
}, static=TRUE, protected=TRUE)


setMethodS3("byGenome", "AromaGenomeTextFile", function(static, genome, ..., mustExist=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'mustExist':
  mustExist <- Arguments$getLogical(mustExist);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Locating genome annotation file");

  # Locate genome file
  pathname <- findByGenome(static, genome=genome, ..., verbose=verbose);

  if (is.null(pathname)) {
    res <- NULL;
  } else {
    res <- newInstance(static, pathname);
  }

  if (mustExist && is.null(res)) {
    msg <- sprintf("Failed to locate a genome annotation data file for genome '%s'.", genome);
    verbose && cat(verbose, msg);
    throw(msg);
  }

  verbose && exit(verbose);

  res;
}, static=TRUE)


############################################################################
# HISTORY:
# 2014-02-03
# o BUG FIX: readDataFrame() for AromaGenomeTextFile explicitly passed
#   arguments '...' to NextMethod(), which would cause them to be
#   duplicated in certain cases.
# 2011-03-03
# o Replaced argument 'onMissing' of byGenome() with 'mustExist'.
# o Now findByGenome() for AromaGenomeTextFile follows the new aroma
#   search conventions.
# o Updated the default filename patterns used by findByGenome() for
#   AromaGenomeTextFile to "^%s,chromosomes(|,.*)*[.]txt$".  Before
#   additional tags were needed.
# 2010-08-22
# o Added Rdoc comments and an example.
# o Created.
############################################################################
