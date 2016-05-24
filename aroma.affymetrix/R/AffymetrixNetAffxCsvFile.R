# @author "HB"
setConstructorS3("AffymetrixNetAffxCsvFile", function(..., .verify=TRUE) {
  this <- extend(AffymetrixCsvFile(..., .verify=FALSE),
                  c("AffymetrixNetAffxCsvFile", uses("UnitNamesFile")),
    "cached:.unitNames" = NULL
  );

  if (.verify)
    verify(this, ...);
  this;
})


setMethodS3("getDefaultExtension", "AffymetrixNetAffxCsvFile", function(static, ...) {
  "annot.csv";
}, static=TRUE, protected=TRUE)


setMethodS3("findByChipType", "AffymetrixNetAffxCsvFile", function(static, chipType, tags=".*", pattern=sprintf("^%s%s([.]|_)%s$", chipType, tags, getDefaultExtension(static)), ...) {
  NextMethod("findByChipType", chipType=chipType, pattern=pattern);
}, static=TRUE, protected=TRUE)



setMethodS3("readUnitNames", "AffymetrixNetAffxCsvFile", function(this, colClasses=c("*"="NULL", "^probe[sS]etI[dD]$"="character"), con=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Reading unitName from file");

  data <- readDataFrame(this, colClasses=colClasses, ..., verbose=less(verbose));

  data <- data[[1]];
  attr(data, "importNames") <- colnames(data);

  verbose && exit(verbose);

  data;
}, protected=TRUE)


setMethodS3("getUnitNames", "AffymetrixNetAffxCsvFile", function(this, ..., force=FALSE) {
  unitNames <- this$.unitNames;

  if (force || is.null(unitNames)) {
    unitNames <- readUnitNames(this, ...);
  }

  unitNames;
})



setMethodS3("readDataUnitChromosomePosition", "AffymetrixNetAffxCsvFile", function(this, colClasses=c("*"="NULL", "^probe[sS]etI[dD]$"="character", "^chromosome$"="character", "^(physicalPosition|chromosomeStart|probeStartPosition)$"="character"), con=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Reading (unitName, fragmentLength) from file");

  data <- readDataFrame(this, colClasses=colClasses, ..., verbose=less(verbose));

  # Convert chromosome strings to integers
  cc <- grep("^chr", colnames(data))[1];
  if (length(cc) == 0 || is.na(cc)) {
    throw("Failed to locate chromosome column.");
  }

  map <- c(X=23, Y=24, M=25, MT=25, Z=25);
  values <- data[[cc]];
  for (kk in seq_along(map)) {
    idxs <- which(values == names(map)[kk]);
    values[idxs] <- map[kk];
  }
  values <- as.integer(values);
  data[[cc]] <- values;
  # Not needed anymore
  values <- idxs <- NULL;
  gc <- gc();


  # Convert positions to integers
  cc <- grep("(p|P)os", colnames(data));
  if (length(cc) == 0)
    cc <- grep("(s|S)tart", colnames(data));
  if (length(cc) == 0 || is.na(cc))
    throw("Failed to locate position column.");
  data[[cc]] <- as.integer(data[[cc]]);
  gc <- gc();

  attr(data, "importNames") <- colnames(data);
  colnames(data) <- c("unitName", "chromosome", "position");
  attr(data, "header") <- NULL;

  verbose && exit(verbose);

  data;
}, protected=TRUE);



setMethodS3("readDataUnitFragmentLength", "AffymetrixNetAffxCsvFile", function(this, colClasses=c("*"="NULL", "^probe[sS]etI[dD]$"="character", "^fragment.*Length.*"="character"), enzymes=1, con=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'enzymes':
  if (is.numeric(enzymes)) {
    enzymes <- Arguments$getIndices(enzymes, max=10);
  } else {
    enzymes <- Arguments$getCharacters(enzymes);
  }
  if (any(duplicated(enzymes))) {
    throw("Argument 'enzymes' contains duplicated values: ",
                                          paste(enzymes, collapse=", "));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Reading (unitName, fragmentLength+) from file");

  data <- readDataFrame(this, colClasses=colClasses, ..., verbose=less(verbose));

  # Extract fragment lengths
  verbose && enter(verbose, "Extracting fragment lengths from ([enzyme], lengths, start, stop)");

  cc <- grep("^fragment.*Length", colnames(data))[1];
  fln <- data[[cc]];
  data[[cc]] <- NULL;
  # Remove all white spaces
  fln <- gsub(" ", "", fln, fixed=TRUE);
  # Replace '///' with ';'
  fln <- gsub("///", ";", fln, fixed=TRUE);
  # Replace '//' with ','
  fln <- gsub("//", ",", fln, fixed=TRUE);
  if (isVisible(verbose, level=-50))
    verbose && str(verbose, fln[1:min(10,length(fln))], level=-50);
  gc <- gc();


  # Are enzyme names specified?
  verbose && enter(verbose, "Inferring if enzyme names are specified");
  hasNames <- NA;
  for (kk in seq_along(fln)) {
    unit <- fln[kk];
    if (nzchar(unit)) {
      hasNames <- (regexpr("^([abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ]+|---)", unit) != -1);
      break;
    }
  }
  if (is.na(hasNames))
    throw("INTERNAL ERROR: Failed to parse CSV file for fragment lengths.");
  verbose && cat(verbose, "Has enzyme names: ", hasNames);
  verbose && exit(verbose);


  # Split by enzymes first
  fln <- strsplit(fln, split=";", fixed=TRUE);
  if (isVisible(verbose, level=-50))
    verbose && str(verbose, fln[1:min(10,length(fln))], level=-50);
  gc <- gc();

  if (hasNames) {
    verbose && enter(verbose, "Identifying number of enzymes");
    nbrOfEnzymes <- sapply(fln, length);
    verbose && print(verbose, table(nbrOfEnzymes));
    nbrOfEnzymes <- max(nbrOfEnzymes);
    verbose && cat(verbose, "Max number of enzymes: ", nbrOfEnzymes);
    verbose && exit(verbose);

    verbose && enter(verbose, "Splitting into subunits and padding with NAs");
    keep <- 1:nbrOfEnzymes;
    fln <- lapply(fln, FUN=function(unit) {
       # Pad missing enzymes with trailing NAs
       unit <- .subset(unit, keep);
       # Split values
       strsplit(unit, split=",", fixed=TRUE);
    });
    if (isVisible(verbose, level=-50))
      verbose && str(verbose, fln[1:min(10,length(fln))], level=-50);
    verbose && exit(verbose);

    # Extract the name for each enzyme
    verbose && enter(verbose, "Extracting enzyme names");
    enzymeIdxs <- lapply(fln, FUN=function(unit) {
      sapply(unit, FUN=.subset, 1, USE.NAMES=FALSE);
    });
    enzymeIdxs <- unlist(enzymeIdxs, use.names=FALSE);
    # Replace '---' with NAs
    enzymeIdxs[enzymeIdxs %in% c("---")] <- NA;

    # Identify the unique enzymes
    allEnzymes <- na.omit(sort(unique(enzymeIdxs)));

    # Map 'enzymes' names to found names
    if (is.character(enzymes)) {
      verbose && enter(verbose, "Mapping requested enzyme names to names in file");
      enzymeNames <- enzymes;
      enzymes <- match(enzymeNames, allEnzymes);
      verbose && cat(verbose, paste(enzymeNames, enzymes, sep=" = "));
      if (any(is.na(enzymes))) {
        throw("Argument 'enzymes' specifies enzyme names that does not exists: ", paste(enzymeNames[is.na(enzymes)], collapse=", "));
      }
      verbose && exit(verbose);
    }

    # AD HOC: In the na24 builds, for Mapping250K_{Nsp|Sty} and
    # the fault GenomeWideSNP_5, the enzyme name field is there
    # but all are '---'.  In that case, assume first enzyme.
    if (length(allEnzymes) == 0) {
      verbose && cat(verbose, "No identified enzymes. Assuming a single enzyme.");
      enzymeIdxs <- rep(1, times=length(enzymeIdxs));
    } else {
      verbose && cat(verbose, "Identified enzymes: ",
                                    paste(allEnzymes, collapse=", "));
      # Map names to indices
      enzymeIdxs <- match(enzymeIdxs, allEnzymes);
    }

    # Sanity check
    if (length(enzymeIdxs) %% nbrOfEnzymes != 0) {
      throw("Internal error: The number of extracted (and NA padded) enzymes IDs (", length(enzymeIdxs), ") is not a multiple of the number of enzymes (", nbrOfEnzymes, ")");
    }

    # Put into an ExJ matrix
    enzymeIdxs <- matrix(enzymeIdxs, nrow=nbrOfEnzymes);
    if (isVisible(verbose, level=-50))
      verbose && str(verbose, enzymeIdxs, level=-50);
    verbose && exit(verbose);


    verbose && enter(verbose, "Identifying the location of the fragment lengths");
    offset <- 1;
    for (kk in seq_along(fln)) {
      unit <- fln[[kk]][[1]];
      if (length(unit) > 1) {
        if (length(unit) == 5) {
          # e.g. GenomeWideSNP_5 for na25
          offset <- 3;
        } else {
          offset <- 2;
        }
        break;
      }
    }

    verbose && cat(verbose, "Offset: ", offset);
    verbose && exit(verbose);

    # Extract the fragment length for each enzyme
    verbose && enter(verbose, "Extracting fragment lengths");
    fln <- lapply(fln, FUN=function(unit) {
      sapply(unit, FUN=.subset, offset, USE.NAMES=FALSE);
    });
    fln <- unlist(fln, use.names=FALSE);
    fln <- as.integer(fln);
    verbose && cat(verbose, "Summary of *all* fragment lengths:");
    verbose && summary(verbose, fln);

    # Sanity check
    if (length(fln) %% nbrOfEnzymes != 0) {
      throw("Internal error: The number of extracted (and NA padded) fragment lengths (", length(fln), ") is not a multiple of the number of enzymes (", nbrOfEnzymes, ")");
    }

    # Put into an ExJ matrix
    fln <- matrix(fln, nrow=nbrOfEnzymes);
    verbose && exit(verbose);


    verbose && enter(verbose, "Sorting data by enzyme");
    # Reorganize as an JxE matrix (transposed compared with 'fln'!)
    naValue <- as.integer(NA);
    fln2 <- matrix(naValue, nrow=ncol(fln), ncol=nbrOfEnzymes);
    for (ee in seq_along(allEnzymes)) {
      for (rr in seq_along(allEnzymes)) {
        # Identify all indices that have enzyme 'ee' in row 'rr'
        idxs <- which(enzymeIdxs[rr,] == ee);
        if (length(idxs) > 0) {
          values <- fln[rr,idxs];
          ok <- is.finite(values);
          fln2[idxs[ok],ee] <- values[ok];
          # Not needed anymore
          ok <- values <- NULL;
        }
      }
    } # for (ee ...)
    colnames(fln2) <- c(allEnzymes, rep(NA, ncol(fln2)-length(allEnzymes)));
    verbose && str(verbose, fln2);

    # Keep only requested enzymes
    verbose && summary(verbose, fln2);
    fln <- fln2[,enzymes,drop=FALSE];
#    verbose && summary(verbose, fln);
    # Not needed anymore
    enzymeIdxs <- fln2 <- NULL;
    verbose && exit(verbose);
  } else {
    nbrOfEnzymes <- length(enzymes);
    # Extract the fragment length for each enzyme
    verbose && enter(verbose, "Extracting fragment lengths");
    fln <- lapply(fln, FUN=function(unit) {
      # Keep only requested enzymes
      unit <- .subset(unit, enzymes);
      parts <- strsplit(unit, split=",", fixed=TRUE);
      sapply(parts, FUN=.subset, 1, USE.NAMES=FALSE);
    });
    verbose && exit(verbose);

    # Reorganize as an JxE integer matrix
    fln <- unlist(fln, use.names=FALSE);
    fln <- as.integer(fln);
    fln <- matrix(fln, ncol=nbrOfEnzymes, byrow=TRUE);
  }

  if (isVisible(verbose, level=-10))
    verbose && str(verbose, fln, level=-10);

  # Sanity check
  if (nrow(fln) != nrow(data)) {
    throw("Internal error. nrow(fln) != nrow(data): ",
                                       nrow(fln), " != ", nrow(data));
  }

  verbose && exit(verbose);

  names <- rep("fragmentLength", ncol(fln));
  if (ncol(fln) > 1)
    names[-1] <- sprintf("%s.%02d", names[-1], 2:ncol(fln));

  # Keep only enzymes of interest
#  names <- names[enzymes];

  data <- data.frame(unitName=data[[1]], fln);
  colnames(data) <- c("unitName", names);
  attr(data, "enzymeNames") <- colnames(fln);
  attr(data, "header") <- NULL;

  verbose && exit(verbose);

  data;
}, protected=TRUE);



############################################################################
# HISTORY:
# 2011-11-19
# o Added getDefaultExtension() for AffymetrixNetAffxCsvFile.
# 2010-02-15
# o BUG FIX: readDataUnitChromosomePosition() of AffymetrixNetAffxCsvFile
#   failed to map chromosome 'MT' to 25.  This bug was introduced
#   in aroma.affymetrix v1.0.7.
# 2009-05-19
# o Now readDataUnitChromosomePosition() for AffymetrixNetAffxCsvFile
#   also recognize 'M' for chromosome 25.
# 2008-09-15
# o Now it is possible to specify enzyme names in argument 'enzymes' to
#   readDataUnitFragmentLength().
# 2008-08-21
# o Now readDataUnitChromosomePosition() also recognizes "MT" as Chr25
#   (mitochondrial).
# 2008-07-20
# o Updated the following methods to preallocate matrixes with the correct
#   data type to avoid coercing later: readDataUnitFragmentLength().
# 2008-04-25
# o Added getUnitNames().
# 2008-04-23
# o NOTE: The CNV for GenomeWideSNP_5 for na25 contains one SNP
#   (SNP_A-1905053) with "duplicated" fragment lengths, i.e. NspI=617,
#   StyI=912, StyI=57.  This will cause the returned data to contain "three"
#   enzyme columns.
# 2007-12-08
# o BUG FIX: readDataUnitFragmentLength() of AffymetrixNetAffxCsvFile would
#   not handle cases where enzymes names a given but all are '---'.
# 2007-11-28
# o Updated readDataUnitFragmentLength() of AffymetrixNetAffxCsvFile to
#   also handle "enzyme data" columns that contain named or non-named
#   multiple enzymes.
# 2007-09-16
# o Now readDataUnitChromosomePosition() of AffymetrixNetAffxCsvFile also
#   recognizes column name 'probeStartPosition'.
# 2007-09-14
# o Added support to read fragment lengths for each enzyme.
# 2007-09-11
# o Added readDataUnitFragmentLength().
# 2007-09-10
# o Created from AffymetrixCsvGenomeInformation.R.
############################################################################
