###########################################################################/**
# @set "class=ChipEffectFile"
# @RdocMethod fromDataFile
#
# @title "Retrive an existing CEL file, or create from CDF template if missing"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{static}{}
#   \item{df}{}
#   \item{filename}{The filename of the CEL file.}
#   \item{path}{The path to the directory where to find/create the CEL file.}
#   \item{name}{The name of the array to be stored in the CEL header.}
#   \item{cdf}{The template @see "AffymetrixCdfFile" used for creating
#              a CEL file from scratch.}
#   \item{...}{Passed to @see "affxparser::createCel".}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @see "ChipEffectFile".
# }
#
# @author "HB"
#
# \seealso{
#   \code{allocateFromCdf()} of @see "AffymetrixCelFile".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fromDataFile", "ChipEffectFile", function(static, df=NULL, filename=sprintf("%s,chipEffects.CEL", getFullName(df)), path, name=getName(df), cdf=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'df':
  if (!is.null(df)) {
    df <- Arguments$getInstanceOf(df, "AffymetrixCelFile");
  }

  # Argument 'cdf':
  if (is.null(cdf)) {
    if (is.null(df))
      throw("Either argument 'df' or 'cdf' must specified.");
  } else {
    cdf <- Arguments$getInstanceOf(cdf, "AffymetrixCdfFile");
  }

  # Argument 'filename' & 'path':
  # First, remove any replicated "chipEffects" tags.
  filename <- gsub("(,chipEffects)*,chipEffects", ",chipEffects", filename);
  pathname <- Arguments$getWritablePathname(filename, path=path);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # Rename lower-case *.cel to *.CEL, if that is the case.  Older versions
  # of the package generated lower-case CEL files. /HB 2007-08-09
  pathname <- AffymetrixFile$renameToUpperCaseExt(pathname);


  # Does the file have to be created?
  if (!isFile(pathname)) {
    verbose && enter(verbose, "Allocating empty chip-effect file");
    verbose && cat(verbose, "Pathname: ", pathname);

    # Get CDF for chip effects
    if (is.null(cdf)) {
      cdf <- createParamCdf(static, getCdf(df), verbose=less(verbose));
    }

    # Get CDF header
    cdfHeader <- getHeader(cdf);

    # Build a valid CEL header
    celHeader <- .cdfHeaderToCelHeader(cdfHeader, sampleName=name);

    # Add some extra information about what the CEL file is for
    params <- c(Descripion="This CEL file contains chip-effect estimates from the aroma.affymetrix package.");
    parameters <- gsub(" ", "_", params);
    names(parameters) <- names(params);
    parameters <- paste(names(parameters), parameters, sep=":");
    parameters <- paste(parameters, collapse=";");
    parameters <- paste(celHeader$parameters, parameters, "", sep=";");
    parameters <- gsub(";;", ";", parameters);
    parameters <- gsub(";$", "", parameters);
    celHeader$parameters <- parameters;

    # Write to a temporary file
    pathnameT <- pushTemporaryFile(pathname, verbose=verbose);

    # Create the CEL file
    .createCel(pathnameT, header=celHeader, ..., verbose=less(verbose));

##    # Fill with negative values
##    nbrOfProbes <- celHeader$total;
##    .updateCel(pathnameT, indices=1:nbrOfProbes, intensities=rep(-1,nbrOfProbes), verbose=less(verbose));

    # Rename temporary file
    popTemporaryFile(pathnameT, verbose=verbose);

    verbose && exit(verbose);
  }

  verbose && enter(verbose, "Setting up ", class(static)[1]);
  verbose && cat(verbose, "Pathname: ", pathname);
  res <- newInstance(static, pathname);

  # Inherit the CDF?
  if (!is.null(cdf))
    setCdf(res, cdf);

  verbose && exit(verbose);

  res;
}, static=TRUE, private=TRUE)



############################################################################
# HISTORY:
# 2010-05-12
# o ROBUSTNESS: When fromDataFile() of ChipEffectFile creates a file, it
#   is created first as a temporary file which is then renamed.  This
#   lowers the risk of generating incomplete chip-effect files.
# 2008-03-18
# o TO DO: Use new static allocateFromCdf() of AffymetrixCelFile.
# o Removed the backward compatibility patch from 2007-01-10 that made
#   fromDataFile() of ChipEffectFile to add missing tags. If anyone still
#   has such old chip effect files lying around, they have to either add
#   the tags manually or reanalyze the data if they want the fullnames
#   of the chip effect files to match the fullnames raw data files.
# o Extracted/moved from ChipEffectFile.R.
############################################################################
