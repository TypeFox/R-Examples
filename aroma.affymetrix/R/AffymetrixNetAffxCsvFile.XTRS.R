setMethodS3("readGeneAssignments", "AffymetrixNetAffxCsvFile", function(this, ..., unique=TRUE, parse=TRUE, fields=NULL, flatten=TRUE, na.rm=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'unique':
  unique <- Arguments$getLogical(unique);

  # Argument 'unique':
  parse <- Arguments$getLogical(parse);

  knownFields <- c("accession", "geneSymbol", "geneTitle", "cytoband", "entrezGeneId");
  if (!is.null(fields)) {
    fields <- match.arg(fields, choices=knownFields, several.ok=TRUE);
  }

  # Argument 'flatten':
  flatten <- Arguments$getLogical(flatten);

  # Argument 'na.rm':
  na.rm <- Arguments$getLogical(na.rm);
  na.rm <- (na.rm || flatten);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Reading gene assignments from ", class(this)[1]);
  verbose && print(verbose, this);

  # Force unique?
  if (flatten && !unique) {
    unique <- TRUE;
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading annotation data");
  colClasses <- c("*"="NULL", "(transcriptClusterId|probesetId|geneAssignment)"="character");
  colClasses <- c("*"="NULL", "(transcriptClusterId|probesetId|unigene)"="character");
  map <- readDataFrame(this, colClasses=colClasses, ...);
  verbose && cat(verbose, "Number of entries : ", nrow(map));
  verbose && cat(verbose, "Number of unique probesetIds: ", length(unique(map$probesetId)));
  verbose && cat(verbose, "Number of unique transcriptClusterIds: ", length(unique(map$transcriptClusterId)));
  verbose && exit(verbose);

  if (unique) {
    verbose && enter(verbose, "Dropping duplicated entries");

    n0 <- nrow(map);
    map <- unique(map);
    n <- nrow(map);
    if (n < n0) {
      verbose && printf(verbose, "Dropped %d (%.2f%%) out of %d duplicated entries\n", n-n0, 100*(n-n0)/n0, n0);
      verbose && cat(verbose, "Number of unique entries: ", nrow(map));
    } else {
      verbose && printf(verbose, "All %d entries are unique.\n", n);
    }

    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Parse
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (parse) {
    verbose && enter(verbose, "Parsing entries");
    
    rows <- seq_len(nrow(map));
    pairs <- map$geneAssignment;
    verbose && cat(verbose, "Number of entries: ", length(rows));

    pairs <- strsplit(pairs, split=" /// ", fixed=TRUE);
    if (unique) {
      verbose && enter(verbose, "Dropping duplicates");
      pairs <- lapply(pairs, FUN=unique);
      verbose && cat(verbose, "Number of entries: ", length(rows));
      verbose && exit(verbose);
    }

    ok <- !sapply(pairs, FUN=function(x) any(is.na(x)));
    nok <- sum(ok);
    verbose && printf(verbose, "Number of entries with annotations: %d (%.2f%%) out of %d\n", nok, 100*nok/length(ok), length(ok));
    if (na.rm) {
      verbose && enter(verbose, "Dropping entries without annotation");
      rows <- rows[ok];
      pairs <- pairs[ok];
      # Update
      ok <- !sapply(pairs, FUN=function(x) any(is.na(x)));
      # Sanity check
      stopifnot(all(ok));
      verbose && cat(verbose, "Number of entries: ", length(rows));
      verbose && cat(verbose, "Number of unique probesetIds: ", length(unique(map$probesetId[rows])));
      verbose && cat(verbose, "Number of unique transcriptClusterIds: ", length(unique(map$transcriptClusterId[rows])));
      verbose && exit(verbose);
    }


    # Valid number of annotation fields per entry    
    nFields <- c(2L, 5L);

    # Infer number of fields from data?
    if (length(nFields) > 1) {
      verbose && enter(verbose, "Inferring number of annotation fields");
      x <- pairs[ok][[1]];
      x <- strsplit(x, split=" // ", fixed=TRUE);
      ns <- sapply(x, FUN=length);
      n <- unique(ns);
      # Sanity check
      stopifnot(length(n) == 1);
      stopifnot(is.element(n, nFields));
      nFields <- n;
      verbose && cat(verbose, "Number of annotation fields: ", nFields);
      verbose && exit(verbose);
    }
    knownFields <- knownFields[1:nFields];
    verbose && cat(verbose, "Annotation fields: ", hpaste(knownFields, maxHead=Inf));


    verbose && enter(verbose, "Parsing annotation field");
    pairs[ok] <- lapply(pairs[ok], FUN=function(x) {
      x <- strsplit(x, split=" // ", fixed=TRUE);
      ns <- sapply(x, FUN=length);
      # Sanity check
      stopifnot(all(ns == nFields));
      x <- unlist(x, use.names=FALSE);
      dimNA(x) <- c(nFields,NA);
      t(x);
    });
    verbose && exit(verbose);


    # Remove any duplicates?
    if (unique) {
      verbose && enter(verbose, "Removing duplicates");
      pairs[ok] <- lapply(pairs[ok], FUN=unique);
      verbose && exit(verbose);
    }


    verbose && enter(verbose, "Adding column names");
    verbose && cat(verbose, "Column names: ", hpaste(knownFields, maxHead=Inf));
    # [1] HuGene-1_0-st-v1.na31.AFFX_README.NetAffx-CSV-Files.txt
    ok <- !sapply(pairs, FUN=function(x) any(is.na(x)));
    pairs[ok] <- lapply(pairs[ok], FUN=function(x) {
      colnames(x) <- knownFields;
      x;
    });
    verbose && exit(verbose);


    # Check argument 'fields' again
    if (is.null(fields)) {
      fields <- knownFields;
    }

    # Extract a subset of fields?
    if (!all(fields == knownFields)) {
      verbose && enter(verbose, "Extracting fields of interest");
      verbose && cat(verbose, "Fields: ", hpaste(fields, maxHead=Inf));

      ok <- !sapply(pairs, FUN=function(x) any(is.na(x)));
      pairs[ok] <- lapply(pairs[ok], FUN=function(x) {
        x[,fields,drop=FALSE];
      });

      if (unique) {
        verbose && enter(verbose, "Dropping duplicated entries");
        pairs[ok] <- lapply(pairs[ok], FUN=unique);
        verbose && exit(verbose);
      }

      verbose && exit(verbose);
    }

    # Sanity check
    stopifnot(length(pairs) == length(rows));


    # Flatten?
    if (flatten) {
      verbose && enter(verbose, "Flattens list to table");
      verbose && cat(verbose, "Number of entries: ", length(rows));

      # Sanity check
      stopifnot(all(ok));

      verbose && enter(verbose, "Identifying blocks of unique sizes");
      ns <- sapply(pairs, FUN=NROW, simplify=TRUE, USE.NAMES=FALSE);
      uns <- sort(unique(ns));
      verbose && print(verbose, uns);
      verbose && exit(verbose);

      verbose && enter(verbose, "Stacking by block size");
      ids <- map$transcriptClusterId;
      ids <- map$probesetId;
      ids <- ids[rows];

      verbose && cat(verbose, "Number of entries: ", length(rows));
      verbose && cat(verbose, "Number of ids: ", length(ids));
      verbose && cat(verbose, "Number of unique ids: ", length(unique(ids)));

      unitNames <- idxs <- c();
      for (ii in seq_along(uns)) {
        n <- uns[ii];
        verbose && enter(verbose, sprintf("Size %d (n=%d) of %d", ii, n, length(uns)));
        keep <- which(ns == n);
        verbose && str(verbose, keep);

        nII <- length(keep);
        unitNamesII <- ids[keep];
        unitNamesII <- rep(unitNamesII, each=n);
        idxsII <- rep(seq_len(n), times=nII);

        unitNames <- c(unitNames, unitNamesII);
        idxs <- c(idxs, idxsII);

        verbose && exit(verbose);
      } # for (ii ...)
  
      # Sanity checks
      stopifnot(length(unitNames) == length(idxs));
      stopifnot(length(idxs) == sum(ns));
      verbose && exit(verbose);

  
      verbose && enter(verbose, "Stacking unit entries");
      dataList <- pairs;
      verbose && cat(verbose, "Number of entries: ", length(dataList));
      chunkSizes <- c(3000, rep(5, times=100));
      while (length(dataList) > 1) {
        len <- length(dataList);
        chunkSize <- chunkSizes[1];
        starts <- seq(from=1, to=len, by=chunkSize);
        ends <- starts + chunkSize - 1L;
        ends[length(ends)] <- len;
        dataList2 <- vector("list", length(starts));
        for (kk in seq_along(starts)) {
          verbose && enter(verbose, sprintf("Chunk #%d of %d", kk, length(starts)));
          idxs <- starts[kk]:ends[kk];
          dataKK <- Reduce(rbind, dataList[idxs]);
          dataList2[[kk]] <- dataKK;
          verbose && exit(verbose);
        } # for (kk ...)
        dataList <- dataList2;
        chunkSizes <- chunkSizes[-1];
      } # for (chunkSize ...)

      data <- dataList[[1]];
      rownames(data) <- NULL;
      verbose && exit(verbose);
  
      verbose && enter(verbose, "Building final table");
      # Sanity check
      stopifnot(length(unitNames) == nrow(data));
      data <- cbind(unitName=unitNames, index=idxs, data);
      rownames(data) <- NULL;
      verbose && exit(verbose);

      verbose && enter(verbose, "Coercing to data frame");
      data <- as.data.frame(data, stringsAsFactors=FALSE);
      verbose && exit(verbose);
  
      pairs <- data;
  
      verbose && exit(verbose);
    } # if (flatten)

    res <- pairs;

    verbose && exit(verbose);
  } # if (parse)

  verbose && exit(verbose);

  res;
})  # readGeneAssignments()


setMethodS3("getHeaderAttributes", "AffymetrixNetAffxCsvFile", function(this, ...) {
  hdr <- getHeader(this, ...);
  comments <- hdr$comments;
  attrs <- grep("^#%", comments, value=TRUE);
  attrs <- gsub("^#%", "", attrs);
  pattern <- "([^=]*)=(.*)";
  keys <- gsub(pattern, "\\1", attrs);
  values <- gsub(pattern, "\\2", attrs);
  values <- trim(values);
  names(values) <- keys;
  attrs <- as.list(values);
  attrs;
}) # getHeaderAttributes()


setMethodS3("getGenomeBuild", "AffymetrixNetAffxCsvFile", function(this, ...) {
  attrs <- getHeaderAttributes(this, ...);
  attrs <- attrs[grep("^genome-version", names(attrs))];
  if (length(attrs) == 0) {
    return(NULL);
  }

  keys <- names(attrs);
  if (is.element("genome-version", keys)) {
    res <- attrs[["genome-version"]];
  } else if (is.element("genome-version-ucsc", keys)) {
    res <- attrs[["genome-version-ucsc"]];
  } else {
    res <- attrs[[1]];
  }
  res; 
})

setMethodS3("getNetAffxBuild", "AffymetrixNetAffxCsvFile", function(this, ...) {
  attrs <- getHeaderAttributes(this, ...);
  res <- attrs[["netaffx-annotation-netaffx-build"]];
  res;
})


setMethodS3("getNetAffxDate", "AffymetrixNetAffxCsvFile", function(this, ...) {
  attrs <- getHeaderAttributes(this, ...);
  res <- attrs[["netaffx-annotation-date"]];
  res <- as.Date(res);
  res;
})


##############################################################################
# HISTORY:
# 2011-09-11
# o Added getGenomeBuild(), getNetAffxBuild() and getNetAffxDate().
# o Added getHeaderAttributes() for AffymetrixNetAffxCsvFile.
# 2011-04-18
# o Updated readGeneAssignments().
# 2011-04-07
# o Now readGeneAssignments() for AffymetrixNetAffxCsvFile handles both
#   *.probeset.csv (2 fields) and *.transcript.csv (5 fields), at least
#   for the HuGene-1_0-st-v1 chip type.
# 2011-04-06
# o Added readGeneAssignments() for AffymetrixNetAffxCsvFile.
# o Created.
##############################################################################
