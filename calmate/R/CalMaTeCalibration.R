###########################################################################/**
# @RdocClass CalMaTeCalibration
#
# @title "The CalMaTeCalibration class"
#
# \description{
#  @classhierarchy
#
#  This class represents the CalMaTe method [1], which 
#  corrects for SNP effects in allele-specific copy-number estimates
#  (ASCNs).
# }
# 
# @synopsis 
#
# \arguments{
#   \item{data}{A named @list with data set named \code{"total"} and
#     \code{"fracB"} where the former should be of class
#     @see "aroma.core::AromaUnitTotalCnBinarySet" and the latter of
#     class @see "aroma.core::AromaUnitFracBCnBinarySet".  The
#     two data sets must be for the same chip type, have the same
#     number of samples and the same sample names.}
#   \item{tags}{Tags added to the output data sets.}
#   \item{references}{An optional @numeric @vector specifying which samples
#     should be as reference samples for estimating the model parameters.
#     If @NULL, all samples are used.}
#   \item{flavor}{A @character string specifying which flavor of the
#     CalMaTe algorithm to use for fitting the model.
#     See @see "fitCalMaTeInternal" for details.}
#   \item{...}{Additional arguments passed to @see "calmateByTotalAndFracB".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
# 
# \section{Reference samples}{
#  In order to estimate the calibration parameters, the model assumes
#  that, for any given SNP, there are a majority of samples that are
#  diploid at that SNP.  Note that it does not have to be the same set
#  of samples for all SNPs.
#
#  By using argument \code{references}, it is possible so specify which
#  samples should be used when estimating the calibration parameters.
#  This is useful when for instance there are several tumor samples with
#  unknown properties as well as a set of normal samples that can be
#  assumed to be diploid.
#
#  Theoretical, a minimum of three reference samples are needed in order
#  for the model to be identifiable.  If less, an error is thrown.
#  However, in practice more reference samples should be used, that is,
#  in the order of at least 6-10 reference samples with a diverse set
#  of genotypes.
# }
#
# \section{Flavors}{
#   For backward compatibility, we try to keep all major versions of
#   the CalMaTe algorithm available.  Older versions can be used by
#   specifying argument \code{flavor}.
#   For more information about the different flavors, 
#   see @see "fitCalMaTeInternal".
# }
#
# \examples{\dontrun{
#   @include "../incl/CalMaTeCalibration.Rex"
# }}
#
# \references{
#  [1] @include "../incl/OrtizM_etal_2012.Rd" \cr 
# }
#
# \seealso{
#   Low-level versions of the CalMaTe method is available
#   via the @see "calmateByThetaAB" and 
#   @see "calmateByTotalAndFracB" methods.
#
#  For further information on the internal fit functions, see
#  @see "fitCalMaTeInternal".
# }
#
#*/###########################################################################
setConstructorS3("CalMaTeCalibration", function(data=NULL, tags="*", references=NULL, flavor=c("v2", "v1"), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(data)) {
    if (!is.list(data)) {
      throw("Argument 'data' is not a list: ", class(data)[1]);
    }
    reqNames <- c("total", "fracB");
    ok <- is.element(reqNames, names(data));
    if (!all(ok)) {
      throw(sprintf("Argument 'data' does not have all required elements (%s): %s", paste(reqNames, collapse=", "), paste(reqNames[!ok], collapse=", ")));
    }
    data <- data[reqNames];

    # Assert correct classes
    className <- "AromaUnitTotalCnBinarySet";
    ds <- data$total;
    if (!inherits(ds, className)) {
      throw(sprintf("The 'total' data set is not of class %s: %s", className, class(ds)[1]));
    }

    className <- "AromaUnitFracBCnBinarySet";
    ds <- data$fracB;
    if (!inherits(ds, className)) {
      throw(sprintf("The 'fracB' data set is not of class %s: %s", className, class(ds)[1]));
    }

    # Assert that the chip types are compatile
    if (getChipType(data$total) != getChipType(data$fracB)) {
      throw("The 'total' and 'fracB' data sets have different chip types: ", 
            getChipType(data$total), " != ", getChipType(data$fracB));
    }

    # Assert that the data sets have the same number data files
    nbrOfFiles <- nbrOfFiles(data$total);
    if (nbrOfFiles != nbrOfFiles(data$fracB)) {
      throw("The number of samples in 'total' and 'fracB' differ: ", 
            nbrOfFiles, " != ", nbrOfFiles(data$fracB));
    }

    # Assert that the data sets have the same samples
    if (!identical(getNames(data$total), getNames(data$fracB))) {
      throw("The samples in 'total' and 'fracB' have different names.");
    }
  }

  # Argument 'references':
  if (!is.null(references)) {
    nbrOfFiles <- nbrOfFiles(data$total);
    references <- Arguments$getIndices(references, max=nbrOfFiles);
    references <- unique(references);
    references <- sort(references);

    if (length(references) < 3) {
      throw("Argument 'references' specifies less than three reference samples: ", length(references));
    }
  }

  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Arguments '...'; optional arguments to calmateByTotalAndFracB()
  extraArgs <- list(...);
  if (length(extraArgs) > 0) {
    keys <- names(extraArgs);
    if (is.null(keys)) {
      throw("Optional arguments to CalMaTeCalibration passed via '...' must be named.");
    }

    nok <- which(nchar(keys) == 0);
    if (length(nok) > 0) {
      throw("All arguments to CalMaTeCalibration passed via '...' must be named.");
    }
  }

  this <- extend(Object(), c("CalMaTeCalibration", uses("ParametersInterface")),
    .data = data,
    .references = references,
    .flavor = flavor,
    .extraArgs = extraArgs
  );

  setTags(this, tags);

  this; 
})


setMethodS3("as.character", "CalMaTeCalibration", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);

  dsList <- getDataSets(this);
  s <- c(s, sprintf("Data sets (%d):", length(dsList)));
  for (kk in seq(along=dsList)) {
    ds <- dsList[[kk]];
    s <- c(s, sprintf("<%s>:", capitalize(names(dsList)[kk])));
    s <- c(s, as.character(ds));
  }
  refs <- getReferences(this);
  nbrOfFiles <- nbrOfFiles(this);
  s <- c(s, sprintf("Number of arrays: %d", nbrOfFiles));
  nbrOfRefs <- length(refs);
  if (nbrOfRefs == 0) {
    s <- c(s, "Number of references: <all arrays> (100%)");
  } else {
    s <- c(s, sprintf("Number of references: %d (%.2f%%)", nbrOfRefs, 100*nbrOfRefs/nbrOfFiles));
  }
 
  params <- getParametersAsString(this);
  nparams <- length(params);
  params <- paste(params, collapse=", ");
  s <- c(s, sprintf("Additional parameters: [%d] %s", nparams, params));

  class(s) <- "GenericSummary";
  s;
}, private=TRUE)



setMethodS3("getAsteriskTags", "CalMaTeCalibration", function(this, collapse=NULL, ...) {
  tags <- "CMTN";
  tags <- c(tags, this$.flavor);
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  }
  
  tags;
}, private=TRUE) 


setMethodS3("getName", "CalMaTeCalibration", function(this, ...) {
  dsList <- getDataSets(this);
  ds <- dsList$total;
  getName(ds);
}) 



setMethodS3("getTags", "CalMaTeCalibration", function(this, collapse=NULL, ...) {
  # "Pass down" tags from input data set
  dsList <- getDataSets(this);
  ds <- dsList$total;
  tags <- getTags(ds, collapse=collapse);

  # Get class-specific tags
  tags <- c(tags, this$.tags);

  # Update default tags
  tags[tags == "*"] <- getAsteriskTags(this, collapse=",");

  # Collapsed or split?
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  } else {
    tags <- unlist(strsplit(tags, split=","));
  }

  if (length(tags) == 0)
    tags <- NULL;

  tags;
})


setMethodS3("setTags", "CalMaTeCalibration", function(this, tags="*", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));
    tags <- tags[nchar(tags) > 0];
  }
  
  this$.tags <- tags;
})


setMethodS3("getFullName", "CalMaTeCalibration", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})


setMethodS3("getDataSets", "CalMaTeCalibration", function(this, ...) {
  this$.data;
})
 

setMethodS3("getRootPath", "CalMaTeCalibration", function(this, ...) {
  "totalAndFracBData";
})


setMethodS3("getPath", "CalMaTeCalibration", function(this, create=TRUE, ...) {
  # Create the (sub-)directory tree for the data set

  # Root path
  rootPath <- getRootPath(this);

  # Full name
  fullname <- getFullName(this);

  # Chip type    
  dsList <- getDataSets(this);
  ds <- dsList$total;
  chipType <- getChipType(ds, fullname=FALSE);

  # The full path
  path <- filePath(rootPath, fullname, chipType, expandLinks="any");

  # Verify that it is not the same as the input path
  inPath <- getPath(ds);
  if (getAbsolutePath(path) == getAbsolutePath(inPath)) {
    throw("The generated output data path equals the input data path: ", path, " == ", inPath);
  }

  # Create path?
  if (create) {
    if (!isDirectory(path)) {
      mkdirs(path);
      if (!isDirectory(path))
        throw("Failed to create output directory: ", path);
    }
  }

  path;
})


setMethodS3("nbrOfFiles", "CalMaTeCalibration", function(this, ...) {
  dsList <- getDataSets(this);
  ds <- dsList$total;
  nbrOfFiles(ds);
})


setMethodS3("getReferences", "CalMaTeCalibration", function(this, ...) {
  this$.references;
})

setMethodS3("getParameters", "CalMaTeCalibration", function(this, ...) {
  params <- list();
  params$references <- getReferences(this);
  params$flavor <- this$.flavor;
  params <- c(params, this$.extraArgs);
  params;
}, protected=TRUE);



setMethodS3("getOutputDataSets", "CalMaTeCalibration", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  res <- this$.outputDataSets;
  if (is.null(res)) {
    res <- allocateOutputDataSets(this, ..., verbose=less(verbose, 10));
    this$.outputDataSets <- res;
  }
  res;
})


setMethodS3("allocateOutputDataSets", "CalMaTeCalibration", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 


  verbose && enter(verbose, "Retrieve/allocation output data sets");

  dsList <- getDataSets(this);
  path <- getPath(this);

  res <- list();
  for (kk in seq(along=dsList)) {
    ds <- dsList[[kk]];
    verbose && enter(verbose, sprintf("Data set #%d ('%s') of %d", 
                                 kk, getName(ds), length(dsList)));

    for (ii in seq(ds)) {
      df <- getFile(ds, ii);
      verbose && enter(verbose, sprintf("Data file #%d ('%s') of %d", 
                                        ii, getName(df), nbrOfFiles(ds)));

      filename <- getFilename(df);
      pathname <- Arguments$getReadablePathname(filename, path=path, 
                                                          mustExist=FALSE);

      # Skip?
      if (isFile(pathname)) {
        verbose && cat(verbose, "Already exists. Skipping.");
        verbose && exit(verbose);
        next;
      }

      # Create temporary file (also checks for write permissions)
      pathnameT <- pushTemporaryFile(pathname);

      # Copy source file
      copyFile(getPathname(df), pathnameT);

      # Make it empty by filling it will missing values
      # AD HOC: We should really allocate from scratch here. /HB 2010-06-21
      dfT <- newInstance(df, pathnameT);
      dfT[,1] <- NA;

      # Renaming temporary file
      pathname <- popTemporaryFile(pathnameT);

      verbose && cat(verbose, "Copied: ", pathname);
      verbose && exit(verbose);
    } # for (ii ...)

    dsOut <- byPath(ds, path=path, ..., verbose=less(verbose, 10));

    # AD HOC: The above byPath() grabs all *.asb files. /HB 2010-06-20
    keep <- is.element(sapply(dsOut, FUN=getFilename), sapply(ds, FUN=getFilename));
    dsOut <- extract(dsOut, keep);

    res[[kk]] <- dsOut;
    rm(ds, dsOut);
    
    verbose && exit(verbose);
  } # for (kk ...)

  names(res) <- names(dsList);

  this$.outputDataSets <- res;

  verbose && exit(verbose);

  res;
}, protected=TRUE)




setMethodS3("findUnitsTodo", "CalMaTeCalibration", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 


  verbose && enter(verbose, "Finding units to do");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (a) Check for missing values in total CNs, because then we will 
  #     check both non-polymorphic loci and SNPs.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dsList <- getOutputDataSets(this);
  dsT <- dsList$total;
  verbose && print(verbose, dsT);

  fullnames <- getFullNames(dsT);
  verbose && str(verbose, fullnames);

  o <- order(fullnames, decreasing=TRUE);
  idx <- o[1];
  dfT <- getFile(dsT, idx);
  verbose && print(verbose, dfT);

  # Read all values
  values <- dfT[,1,drop=TRUE];

  nbrOfUnits <- length(values);
  verbose && cat(verbose, "Number of units: ", nbrOfUnits);

  # Identify all missing values
  nokT <- is.na(values);
  unitsT <- which(nokT);
  verbose && printf(verbose, "Number of TCNs missing values: %d (%.2f%%)\n", length(unitsT), 100*length(unitsT)/nbrOfUnits);



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (b) In case the BAFs were not stored but the TCNs were, which may
  #     happen during a user interrupt or power failure and because BAFs
  #     are stored after TCNs, we must make sure to check those as well.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # The last data set is always updated last
  dsB <- dsList$fracB;
  verbose && print(verbose, dsB);

  # The last file (in lexicographic ordering) is always updated last
  # (This should be the same as above, but just in case...)
  fullnames <- getFullNames(dsB);
  verbose && str(verbose, fullnames);

  o <- order(fullnames, decreasing=TRUE);
  idx <- o[1];
  dfB <- getFile(dsB, idx);
  verbose && print(verbose, dfB);

  # Read all values
  values <- dfB[,1,drop=TRUE];

  # Identify all missing values
  nokB <- is.na(values);
  unitsB <- which(nokB);
  verbose && printf(verbose, "Number of BAFs with missing values: %d (%.2f%%)\n", length(unitsB), 100*length(unitsB)/nbrOfUnits);


  # Identify units that *may* be unfitted.
  nok <- (nokT | nokB);
  units <- which(nok);


  # Potentially unfitted units?
  if (length(units) > 0) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # (c) Now, for BAFs we will get missing values for all non-polymorphic 
    #     loci.  The problem is that we don't know which they are (without 
    #     turning to chip type-specific annotation data files which cannot
    #     assume exist).  Instead, we will check the corresponding input
    #     BAF file to see whether those units also have missing values.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    dsList <- getDataSets(this);
    # The last data set is always updated last
    dsB <- dsList$fracB;
    verbose && print(verbose, dsB);
  
    idxs <- indexOf(dsB, getFullName(dfB));
    # Sanity check
    stopifnot(length(idxs) > 0);
  
    idx <- idxs[1];
    dfB <- getFile(dsB, idx);
    verbose && print(verbose, dfB);
  
    values <- dfB[,1,drop=TRUE];
  
    # Identify all missing values
    nokBin <- is.na(values);
    unitsB <- which(nokBin);
    verbose && printf(verbose, "Number of input BAFs with missing values: %d (%.2f%%)\n", length(unitsB), 100*length(unitsB)/nbrOfUnits);
  
    # Conclusions
    nokB <- (nokB & !nokBin);
    unitsB <- which(nokB);
    verbose && printf(verbose, "Number of BAFs with missing values in output but not in input: %d (%.2f%%)\n", length(unitsB), 100*length(unitsB)/nbrOfUnits);
  
    # Update units that are unfitted.
    nok <- (nokT | nokB);
    units <- which(nok);
  } # if (length(units) > 0)

  verbose && printf(verbose, "Number of units to do: %d (%.2f%%)\n", 
                     length(units), 100*length(units)/length(values));

  verbose && cat(verbose, "Units to do (with missing values):");
  verbose && str(verbose, units);

  verbose && exit(verbose);

  units;
})


setMethodS3("process", "CalMaTeCalibration", function(this, units="remaining", force=FALSE, ram=NULL, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dsList <- getDataSets(this);
  dsTCN <- dsList$total;
  dsBAF <- dsList$fracB;
  
  # Argument 'units':
  if (is.null(units)) {
    units <- "remaining";
  }
  df <- getFile(dsTCN, 1);
  nbrOfUnits <- nbrOfUnits(df);
  if (identical(units, "remaining")) {
    units <- seq(length=nbrOfUnits);
  } else {
    units <- Arguments$getIndices(units, max=nbrOfUnits);
  }

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'ram':
  ram <- getRam(aromaSettings, ram);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Argument '...':
  userArgs <- list(...);
  if (length(userArgs) > 0) {
    if (is.element("references", names(userArgs))) {
      throw("Argument 'references' must not be passed via process(), but instead to the CalMaTeCalibration constructor.");
    }

    warning("Do not pass extra arguments via process(), but instead via CalMaTeCalibration: ", paste(names(userArgs), collapse=", "));
  }



  verbose && enter(verbose, "CalMaTe calibration of ASCNs");

  params <- getParameters(this);

  nbrOfFiles <- nbrOfFiles(this);
  verbose && cat(verbose, "Number of arrays: ", nbrOfFiles);

  references <- params$references;
  nbrOfRefs <- length(references);
  if (nbrOfRefs == 0L) nbrOfRefs <- nbrOfFiles;
  verbose && cat(verbose, "Number of references: ", nbrOfRefs);

  verbose && cat(verbose, "Units:");
  verbose && str(verbose, units);

  verbose && cat(verbose, "All parameters:");
  verbose && str(verbose, params);


  # Skip already processed units or not?
  if (!force) {
    verbose && enter(verbose, "Skipping already processed units");
    unitsTodo <- findUnitsTodo(this, verbose=less(verbose, 25));
    verbose && cat(verbose, "Units remaining:");
    verbose && str(verbose, unitsTodo);

    units <- intersect(units, unitsTodo);
    verbose && cat(verbose, "Units:");
    verbose && str(verbose, units);

    rm(unitsTodo);
    verbose && exit(verbose);
  }

  nbrOfUnits <- length(units);
  verbose && printf(verbose, "Number of units to do: %d (%.2f%%)\n", 
                      length(units), 100*length(units)/nbrOfUnits(df));


  chipType <- getChipType(dsTCN, fullname=FALSE);
  verbose && cat(verbose, "Chip type: ", chipType);
  rm(dsList);


  sampleNames <- getNames(dsTCN);
  dimnames <- list(NULL, sampleNames, c("total", "fracB"));

  outPath <- getPath(this);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocate output data sets
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- getOutputDataSets(this, verbose=less(verbose, 5));
  if (nbrOfUnits == 0) {
    verbose && cat(verbose, "No more units to process. Skipping.");
    verbose && exit(verbose);
    return(res);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Process in chunks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calculating number of units to fit per chunk");
  verbose && cat(verbose, "RAM scale factor: ", ram);

  bytesPerChunk <- 100e6;       # 100Mb
  verbose && cat(verbose, "Bytes per chunk: ", bytesPerChunk);

  bytesPerUnitAndArray <- 2*8;  # Just a rough number; good enough?
  verbose && cat(verbose, "Bytes per unit and array: ", bytesPerUnitAndArray);

  bytesPerUnit <- nbrOfFiles * bytesPerUnitAndArray;
  verbose && cat(verbose, "Bytes per unit: ", bytesPerUnit);

  unitsPerChunk <- ram * bytesPerChunk / bytesPerUnit;
  unitsPerChunk <- as.integer(max(unitsPerChunk, 1));
  unitsPerChunk <- min(unitsPerChunk, nbrOfUnits);
  verbose && cat(verbose, "Number of units per chunk: ", unitsPerChunk);

  nbrOfChunks <- ceiling(nbrOfUnits / unitsPerChunk);
  verbose && cat(verbose, "Number of chunks: ", nbrOfChunks);
  verbose && exit(verbose);

  idxs <- 1:nbrOfUnits;
  head <- 1:unitsPerChunk;

  count <- 1L;
  while (length(idxs) > 0) {
    tTotal <- processTime();

    verbose && enter(verbose, sprintf("Chunk #%d of %d", count, nbrOfChunks));
    if (length(idxs) < unitsPerChunk) {
      head <- 1:length(idxs);
    }
    uu <- idxs[head];

    verbose && cat(verbose, "Units: ");
    verbose && str(verbose, units[uu]);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Reading (total,fracB) data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Reading (total,fracB) data");
    total <- extractMatrix(dsTCN, units=units[uu], verbose=less(verbose,25));
    verbose && str(verbose, total);
    fracB <- extractMatrix(dsBAF, units=units[uu], verbose=less(verbose,25));
    verbose && str(verbose, fracB);
    verbose && exit(verbose);

    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Calibration
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Combining into an (total,fracB) array");
    dim <- c(nrow(total), ncol(total), 2);
    data <- c(total, fracB);
    rm(total, fracB);
    data <- array(data, dim=dim, dimnames=dimnames);
    data <- aperm(data, perm=c(1,3,2));
    verbose && str(verbose, data);
    verbose && exit(verbose);

    verbose && enter(verbose, "Calibration");
    args <- params;
    # BACKWARD COMPATIBILITY: Override with user arguments. /HB 2012-02-05
    for (key in names(userArgs)) {
      args[[key]] <- userArgs[[key]];
    }
    args <- c(list(data), args);
    verbose && cat(verbose, "Arguments passed to calmateByTotalAndFracB():");
    verbose && str(verbose, args);
    args$verbose <- verbose;

    dataN <- do.call("calmateByTotalAndFracB", args=args);

    fit <- attr(dataN, "modelFit");
    verbose && str(verbose, fit);
    verbose && str(verbose, dataN);
    verbose && exit(verbose);

    rm(data);  # Not needed anymore
    gc <- gc();
    verbose && print(verbose, gc);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Storing model fit
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Storing model fit");
    verbose && exit(verbose);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Storing calibrated data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Storing calibrated data");
    for (kk in seq(along=res)) {
      ds <- res[[kk]];
      verbose && enter(verbose, sprintf("Data set #%d ('%s') of %d", 
                                           kk, getName(ds), length(res)));

      # Store in lexicograph ordering
      fullnames <- getFullNames(ds);
      idxsT <- order(fullnames, decreasing=FALSE);
      
      for (ii in idxsT) {
        df <- getFile(ds, ii);
        verbose && enter(verbose, sprintf("Data file #%d ('%s') of %d", 
                                        ii, getName(df), nbrOfFiles(ds)));

        signals <- dataN[,kk,ii];
        verbose && cat(verbose, "Signals:");
        verbose && str(verbose, signals);
        df[units[uu],1] <- signals;
        rm(signals);

        verbose && exit(verbose);
      } # for (ii ...)

      verbose && exit(verbose);
    } # for (kk ...)
    verbose && exit(verbose);

    # Next chunk
    idxs <- idxs[-head];
    count <- count + 1;

    # Garbage collection
    tGc <- processTime();
    gc <- gc();
    verbose && print(verbose, gc);

    verbose && exit(verbose);
  } # while(length(idxs) > 0)

  verbose && exit(verbose);

  invisible(res);
})


############################################################################
# HISTORY:
# 2013-01-05 [HB]
# o CLEANUP: Dropped getParametersAsString() since CalMaTeCalibration 
#   now implements ParametersInterface.
# 2012-02-19 [HB]
# o Made findUnitsTodo() for CalMaTeCalibration smarter. Before it would
#   detect all non-polymorphic loci as non-fitted.
# o Added argument 'flavor' to CalMaTeCalibration, which now also add
#   the "flavor" to the default set of tags.
# 2012-02-06 [HB]
# o Now as.character() also reports any optional parameters.
# o Added getParametersAsString().
# 2012-02-05 [HB]
# o DEPRECATED: Now argument 'references' of process() is obsolete and
#   gives an error.  Specify it via CalMaTeCalibration() instead.
# o Added getParameters() for CalMaTeCalibration.
# 2011-07-15 [HB]
# o DOCUMENTATION: Added a section 'Reference samples' to the help
#   of CalMaTeCalibration.
# o DEPRECATED: Argument 'references' of process() of CalMaTeCalibration
#   is deprecated.  Instead, pass it to CalMaTeCalibration().
# o CLEANUP: Tidied up the validation of argument 'references' and
#   improved the corresponding error messages.
# 2011-07-12 [MO]
# o Check that the number of references has to be at least 3.
#   If no references are given, the total number of files has to be
#   at least 3. (Used to be 6).
# 2011-03-12
# o BUG FIX: After recent update, allocateOutputDataSets() would only
#   work for existing data sets, not to create new ones.
# 2011-03-08
# o GENERALIZATION: Now allocateOutputDataSets() for CalMaTeCalibration
#   no longer requires write permissions if the data set already exists.
# o Added argument 'references' to the CalMaTeCalibration constructor.
# 2010-07-31
# o BUG FIX: process() for CalMaTeCalibration would only run one chunk.
# 2010-07-30
# o Renamed CalMaTeNormalization to CalMaTeCalibration.
# 2010-07-22
# o BUG FIX: Now process() for CalMaTeCalibration returns immediately 
#   if there are no units left.
# o BUG FIX: process(..., verbose=TRUE) would give "Error in sprintf("Chunk
#   #%d of %d", count, nbrOfChunks) :  invalid format '%d'; use format %f, 
#   %e, %g or %a for numeric objects".
# 2010-06-29
# o Added support for process(..., force=TRUE) in CalMaTeNormalization.
# 2010-06-23
# o Added support for argument 'references' and 'truncate'.
#   TODO: The asterisk tags should probably reflect these settings.
# o Code cleanup.
# 2010-06-21
# o Added minor Rdoc comments.
# o ROBUSTNESS: Added more assertions to the CalMaTe constructor.
# o ROBUSTNESS: Now process() stores results in lexicograph ordering to
#   assure that the lexicographicly last file is updated last, which is 
#   the file that findUnitsTodo() is querying.
# o Added getOutputDataSets().
# o Added findUnitsTodo().
# o Now allocateOutputDataSets() "clears" the files.
# 2010-06-20
# o First test shows that it seems to work.
# o Created.
############################################################################ 
