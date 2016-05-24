###########################################################################/**
# @set "class=AffymetrixCelFile"
# @RdocMethod allocateFromCdf
#
# @title "Creates an empty CEL file from a template CDF"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{cdf}{A @see "AffymetrixCdfFile" to be used as a template.}
#   \item{name, tags}{The name and the tags of the file created.}
#   \item{path}{The directory where the file is created.}
#   \item{suffix}{Filename suffix.}
#   \item{...}{Arguments passed to @see "affxparser::createCel".}
#   \item{overwrite}{If @FALSE and the file already exists, then an
#      error is thrown.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an @see "AffymetrixCelFile".
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("allocateFromCdf", "AffymetrixCelFile", function(static, cdf, name, tags=NULL, path=".", suffix=".CEL", ..., overwrite=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'cdf':
  cdf <- Arguments$getInstanceOf(cdf, "AffymetrixCdfFile");

  # Argument 'name':
  name <- Arguments$getCharacter(name, nchar=c(1,256), length=c(1,1));
  name <- trim(name);
  name <- Arguments$getCharacter(name, nchar=c(1,256), length=c(1,1));

  # Argument 'tags':
  tags <- Arguments$getCharacters(tags);
  tags <- trim(tags);
  tags <- tags[nzchar(tags)];

  fullname <- paste(c(name, tags), collapse=",");
  parts <- unlist(strsplit(fullname, split=","));
  name <- parts[1];
  tags <- parts[-1];

  # Argument 'suffix':
  suffix <- Arguments$getCharacter(suffix, length=c(1,1));

  filename <- sprintf("%s%s", fullname, suffix);

  # Argument 'filename' & 'path':
  pathname <- Arguments$getWritablePathname(filename, path=path,
                                                mustNotExist=!overwrite);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);



  verbose && enter(verbose, "Creating chip-effect file");
  verbose && cat(verbose, "Pathname: ", pathname);
  verbose && cat(verbose, "Fullname: ", fullname);
  verbose && cat(verbose, "Name: ", name);
  verbose && cat(verbose, "Tags: ", paste(tags, collapse=","));

  # Get CDF header
  cdfHeader <- getHeader(cdf);

  # Build a valid CEL header
  celHeader <- .cdfHeaderToCelHeader(cdfHeader, sampleName=fullname);

  # Add some extra information about what the CEL file is for
  params <- c(Descripion="This CEL file was created by the aroma.affymetrix package.");
  parameters <- gsub(" ", "_", params);
  names(parameters) <- names(params);
  parameters <- paste(names(parameters), parameters, sep=":");
  parameters <- paste(parameters, collapse=";");
  parameters <- paste(celHeader$parameters, parameters, "", sep=";");
  parameters <- gsub(";;", ";", parameters);
  parameters <- gsub(";$", "", parameters);
  celHeader$parameters <- parameters;

  # Overwrite existing file?
  if (overwrite && isFile(pathname)) {
    verbose && enter(verbose, "Removing existing file (overwrite=TRUE)");
    file.remove(pathname);
    if (isFile(pathname))
      throw("Failed to remove existing file: ", pathname);
    verbose && exit(verbose);
  }

  # Create the CEL file
  .createCel(pathname, header=celHeader, ..., verbose=less(verbose));

##    # Fill with negative values
##    nbrOfProbes <- celHeader$total;
##    updateCel(pathname, indices=1:nbrOfProbes,
##         intensities=rep(-1,nbrOfProbes), verbose=less(verbose));

  verbose && enter(verbose, "Setting up ", class(static)[1]);
  verbose && cat(verbose, "Pathname: ", pathname);
  res <- newInstance(static, pathname);

  # Inherit the CDF?
  if (!is.null(cdf))
    setCdf(res, cdf);

  verbose && exit(verbose);

  res;
}, static=TRUE, protected=TRUE)


############################################################################
# HISTORY:
# 2008-03-18
# o Extracted from the static fromDataFile() of ChipEffectFile.
# o Created.
############################################################################
