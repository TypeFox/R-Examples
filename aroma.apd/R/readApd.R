#########################################################################/**
# @RdocDefault readApd
#
# @title "Reads an Affymetrix probe data (APD) file"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# \arguments{
#   \item{filename}{The filename of the APD file.}
#   \item{indices}{An optional @numeric @vector of cell (probe) indices
#     specifying what cells to read.  If @NULL, all are read.}
#   \item{readMap}{A @vector remapping cell indices to file indices.
#     If \code{"byMapType"}, the read map of type according to APD header
#     will be search for and read.  It is much faster to specify the
#     read map explicitly compared with searching for it each time.
#     If @NULL, no map is used.}
#   \item{name}{The name of the data field.
#     If @NULL, the APD header \code{name} is used.  If not specified,
#     it defaults to \code{"intensities"}.}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
#   \item{.checkArgs}{If @TRUE, arguments are validated, otherwise not.}
# }
#
# \value{
#   A named @list with the two elements \code{header} and
#   \code{data}.  The header is in turn a @list structure and
#   the second is a @numeric @vector holding the queried data.
# }
#
# \details{
#   To read one \emph{large} contiguous block of elements is faster than
#   to read individual elements one by one.  For this reason, internally
#   more elements than requested may be read and therefore allocation more
#   memory than necessary.  This means, in worst case \eqn{N} elements
#   may read allocation \eqn{N*8} bytes of \R memory, although only two
#   elements are queried.  However, to date even with the largest arrays
#   from Affymetrix this will still only require tens of megabytes of
#   \emph{temporary} memory.  For instance, Affymetrix Mapping 100K arrays
#   holds 2,560,000 probes requiring 20Mb of temporary memory.
# }
#
# \section{Remapping indices}{
#   Argument \code{readMap} can be used to remap indices.  For instance,
#   the indices of the probes can be reorder such that the probes within
#   a probeset is in a contiguous set of probe indices.  Then, given that
#   the values are stored in such an order, when reading complete probesets,
#   data will be access much faster from file than if the values were
#   scatter all over the file.
#
#   Example of speed improvements.  Reading all 40000 values in units
#   1001 to 2000 of an Affymetrix Mapping 100K Xba chip is more than
#   10-30 times faster with mapping compared to without.
# }
#
# \section{File format}{
#   The file format of an APD file is identical to the file format of an
#   @see "R.huge::FileVector".
# }
#
# @author
#
# \examples{\dontrun{#See ?createApd for an example.}}
#
# \seealso{
#   @see "createApd" and @see "updateApd".
#   See also @see "readApdHeader".
#   To create a cell-index read map from an CDF file, see
#   @see "affxparser::readCdfUnitsWriteMap".
# }
#
# @keyword "file"
# @keyword "IO"
#*/#########################################################################
setMethodS3("readApd", "default", function(filename, indices=NULL, readMap="byMapType", name=NULL, ..., verbose=FALSE, .checkArgs=TRUE) {
  header <- NULL;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'filename':
  if (.checkArgs) {
    filename <- Arguments$getReadablePathname(filename, mustExists=TRUE);
  }

  apd <- FileVector(filename);
  on.exit(close(apd));
  nbrOfProbes <- length(apd);

  if (.checkArgs) {
    # Argument 'indices':
    if (!is.null(indices)) {
      # An APD file has zero- or one-based indices
      indices <- Arguments$getIndices(indices, range=c(1, nbrOfProbes));
    }

    # Argument 'readMap':
    if (is.null(readMap)) {
      # No probe map specified.
    } else if (is.character(readMap)) {
    } else if (is.vector(readMap)) {
      readMap <- Arguments$getIndices(readMap, range=c(1, nbrOfProbes),
                                                     length=nbrOfProbes);
    }

    # Argument 'name':
    if (!is.null(name)) {
      name <- Arguments$getCharacter(name, nchar=c(1,256));
    }

    # Argument 'verbose':
    verbose <- Arguments$getVerbose(verbose);
  }


  if (identical(readMap, "byMapType")) {
    # Assume that all APD file have the same map type.
    header <- readApdHeader(apd);
    mapType <- header$mapType;
    if (!is.null(mapType)) {
      mapType <- trim(mapType);
      if (is.na(mapType) || nchar(mapType) == 0)
        mapType <- NULL;
    }

    readMap <- mapType;
  }

  if (is.character(readMap)) {
    mapFile <- findApdMap(mapType=readMap);
    if (is.null(mapFile))
      throw("No APD map file found for the given map type: ", readMap);

    readMap <- readApdMap(mapFile)$map;
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get APD header and data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- list();
  if (is.null(header))
    header <- readApdHeader(apd);
  res$header <- header;

  # Get the name from the APD header?
  if (is.null(name)) {
    name <- header$name;
    if (is.null(name)) {
      name <- "intensities";
    }
  }

  # Affymetrix probe data is stored in a file vector.
  if (is.null(indices)) {
    # Read all values
    values <- readAllValues(apd, verbose=verbose);
    # Re-map indices?
    if (!is.null(readMap))
      values <- values[readMap];
  } else {
    # Re-map ?
    if (!is.null(readMap)) {
      indices <- readMap[indices];
    }
    # To speed things up a lot, we read the data in one big patch and then
    # keep only those values we are interested in.  However, we only read
    # the necessary range of values, i.e. the first to the last.
    range <- range(indices);
    first <- range[1];
    n <- range[2]-first+1;
    verbose && cat(verbose, sprintf("Reading %d cells from %d to %d.\n", n, range[1], range[2]));
    values <- readContiguousValues(apd, indices=first, lengths=n, verbose=verbose, .checkArgs=FALSE);
    indices <- indices - (first - 1);
    values <- values[indices];
  }
  res[[name]] <- values;

  res;
})


############################################################################
# HISTORY:
# 2009-05-16
# o Updated readApd() to coerce argument 'indices' and 'readMap' to integer
#   indices.  Before it used to coerce to doubles (before updating R.utils).
# 2006-04-08
# o Remove internal rm().
# o When reading the the APD header, the 'apd' argument is passed.
# 2006-03-28
# o Removed argument 'indexOffset'; now affxparser uses one-based indices.
# o Renamed argument 'map' to 'readMap'.
# 2006-03-18
# o Added argument 'indexOffset' and made it one by default (as in R).
# 2006-03-14
# o If not specified, the probe map is now searched for (and read) using
#   the APD header.
# 2006-03-04
# o Removed all gc(). They slow down quite a bit.
# o Added support for index maps.
# 2006-02-27
# o Created by HB.
############################################################################
