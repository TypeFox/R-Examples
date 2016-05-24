# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BEGIN: Affymetrix specific
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# @author "HB"
setMethodS3("allocateFromCdf", "AromaCellSequenceFile", function(static, cdf, path=NULL, tags="*", ...) {
  # Argument 'cdf':
  cdf <- Arguments$getInstanceOf(cdf, "AffymetrixCdfFile");

  # Argument 'tags':
  tags <- Arguments$getTags(tags, collapse=NULL);

  # Generate filename: <chipType>(,tags)*.<ext>
  chipType <- getChipType(cdf);

  # Move chip type tags to the other tags, if asterisk is used.
  parts <- strsplit(chipType, split=",", fixed=TRUE);
  parts <- unlist(parts, use.names=FALSE);
  tags[tags == "*"] <- paste(parts[-1], collapse=",");
  chipType <- parts[1];

  # Exclude 'monocell' tags (AD HOC)
  chipType <- gsub(",monocell", "", chipType);

  # Output path
  if (is.null(path)) {
    path <- file.path("annotationData", "chipTypes", chipType);
  }
  path <- Arguments$getWritablePath(path);

  # Get platform
  platform <- getPlatform(cdf);

  # Number of cells
  nbrOfCells <- nbrOfCells(cdf);

  fullname <- paste(c(chipType, tags), collapse=",");
  ext <- getFilenameExtension(static);
  filename <- sprintf("%s.%s", fullname, ext);

  # Create microarray tabular binary file
  allocate(static, filename=filename, path=path, nbrOfCells=nbrOfCells,
                               platform=platform, chipType=chipType, ...);
}, static=TRUE)



setMethodS3("importFromAffymetrixProbeTabFile", "AromaCellSequenceFile", function(this, srcFile, rows=NULL, ..., onDuplicates=c("warning", "error", "ignore"), keepSequenceLengths=NULL, ram=NULL, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'srcFile':
  if (inherits(srcFile, "AffymetrixProbeTabFile")) {
    # Validate chip type
    chipType <- getChipType(this, fullname=FALSE);
    chipTypeSrc <- getChipType(srcFile, fullname=FALSE);
    if (!identical(chipTypeSrc, chipType)) {
      throw("Argument 'srcFile' has a different chip type: ",
                                            chipTypeSrc, " != ", chipType);
    }
  } else {
    throw("Argument 'srcFile' is not an AffymetrixProbeTabFile: ",
                                                        class(srcFile)[1]);
  }

  # Argument 'rows':
  nbrOfCells <- nbrOfCells(this);
  if (is.null(rows)) {
    rows <- 1:nbrOfCells;
  } else {
    rows <- Arguments$getIndices(rows, max=nbrOfCells);
    rows <- sort(unique(rows));
  }

  # Argument 'onDuplicates':
  onDuplicates <- match.arg(onDuplicates);

  # Argument 'keepSequenceLengths':
  if (!is.null(keepSequenceLengths)) {
    keepSequenceLengths <- Arguments$getIndices(keepSequenceLengths);
  }

  # Argument 'ram':
  ram <- getRam(aromaSettings, ram);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Assert that the CDF can be found
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chipType <- getChipType(this);
  nbrOfCells <- nbrOfCells(this);
  cdf <- NULL;
  tryCatch({
    cdf <- AffymetrixCdfFile$byChipType(chipType, nbrOfCells=nbrOfCells);
  }, error = function(ex) {
    cdf <<- AffymetrixCdfFile$byChipType(chipType, tags=".*", nbrOfCells=nbrOfCells);
  })

  # Figure out the number of columns on the array
  nbrOfColumns <- nbrOfColumns(cdf);


  verbose && enter(verbose, "Importing probe sequences");
  verbose && cat(verbose, "Pathname: ", getPathname(srcFile));
  verbose && cat(verbose, "Chip type: ", getChipType(srcFile));
  verbose && cat(verbose, "Number of cell columns: ", nbrOfColumns);
  verbose && cat(verbose, "Rows:");
  verbose && str(verbose, rows);

  colClasses <- c("probe(X|Y)Pos"="integer", "probeSequence"="character", "targetStrandedness"="character");

  count <- 0;
  CHUNK.SIZE <- as.integer(ram*1.5e6);
  rowOffset <- as.integer(1);
  idxs <- 1:CHUNK.SIZE;
  while (length(rows) > 0) {
    if (length(rows) < CHUNK.SIZE) {
      idxs <- 1:length(rows);
    }

    rowsChunk <- rows[idxs];
    verbose && printf(verbose, "Row: %d-%d\n", min(rowsChunk), max(rowsChunk));
    df <- readDataFrame(srcFile, colClasses=colClasses,
                        rows=rowsChunk, ..., verbose=less(verbose, 25));
    if (nrow(df) == 0)
      break;
    verbose && cat(verbose, "Data read:");
    verbose && str(verbose, df);

    if (nrow(df) < length(idxs)) {
      idxs <- 1:nrow(df);
    }

    cells <- nbrOfColumns*df[,"probeYPos"] + df[,"probeXPos"] + as.integer(1);

    # Check for duplicated cell indices
    dups <- duplicated(cells);
    hasDuplicates <- any(dups);
    if (hasDuplicates) {
      setOfDups <- cells[dups];
      n <- length(setOfDups);
      msg <- paste("Identified ", n, " duplicated cell indices: ",
                             paste(head(setOfDups), collapse=", "), sep="");
      verbose && cat(verbose, msg);
      verbose && print(verbose, head(df[dups,]));
      if (onDuplicates == "error") {
        throw(msg);
      } else if (onDuplicates == "warning") {
        warning(msg);
      } else if (onDuplicates == "ignore") {
      }
    }
    # Not needed anymore
    dups <- NULL;

    # Clean up to save memory
    df[["probeXPos"]] <- df[["probeYPos"]] <- NULL;
    seqs <- df[,"probeSequence"];
    strands <- df[,"targetStrandedness"];
    # Not needed anymore
    df <- NULL;
    gc <- gc();
    verbose && str(verbose, cells);
    verbose && str(verbose, seqs);
    verbose && str(verbose, strands);

    verbose && cat(verbose, "Unique probe-sequence lengths:");
    verbose && str(verbose, seqs);
    n <- nchar(seqs);
    verbose && print(verbose, table(n));

    # Identify probe sequences to be kept?
    if (!is.null(keepSequenceLengths)) {
      keep <- is.element(n, keepSequenceLengths);
      keep <- which(keep);
      if (length(keep) != length(seqs)) {
        verbose && enter(verbose, "Dropping probe sequence with odd lengths");
        msg <- paste("Dropped ", length(seqs)-length(keep),
                     " sequences, because they are of unwanted lengths: ",
                     paste(keepSequenceLengths, collapse=", "), sep="");
        verbose && cat(verbose, msg);
        cells <- cells[keep];
        seqs <- seqs[keep];
        strands <- strands[keep];
        warning(msg);
        verbose && exit(verbose);
      }
      # Not needed anymore
      keep <- NULL;
    }

    updateSequences(this, cells=cells, seqs=seqs, verbose=less(verbose, 25));
    updateTargetStrands(this, cells=cells, strands=strands, verbose=less(verbose, 25));
    count <- count + length(cells);

    # Next set of rows
    rows <- rows[-idxs];
  } # while(TRUE)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update footer with import information
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Updating file footer:");
  footer <- readFooter(this);
  keys <- paste("srcFile", c("", 2:99), sep="");
  key <- keys[(!keys %in% names(footer))][1];
  footer[[key]] <- list(
    filename = getFilename(srcFile),
    filesize = getFileSize(srcFile),
    checksum = getChecksum(srcFile)
  );
  writeFooter(this, footer);
  verbose && exit(verbose);

  verbose && exit(verbose);

  invisible(as.integer(count));
})



setMethodS3("inferMmFromPm", "AromaCellSequenceFile", function(this, cdf, units=NULL, ..., safe=FALSE, ram=NULL, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Contract: X == dnaComplement(dnaComplement(X)) is always true
  dnaComplement <- function(X, ...) {
    # A (0x01) -> T (0x04)
    # C (0x02) -> G (0x03)
    # G (0x03) -> C (0x02)
    # T (0x04) -> A (0x01)
    from <- as.raw(1:4);
    to <- as.raw(4:1);

    Xc <- X;
    for (kk in 1:4) {
      idxs <- which(X == from[kk]);
      Xc[idxs] <- to[kk];
    }

    Xc;
  } # dnaComplement()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'cdf':
  cdf <- Arguments$getInstanceOf(cdf, "AffymetrixCdfFile");

  chipType <- getChipType(this, fullname=FALSE);
  if (getChipType(cdf, fullname=FALSE) != chipType) {
    throw("Argument 'cdf' is of a different chip type than this: ",
                       getChipType(cdf, fullname=FALSE), " != ", chipType);
  }

  # Argument 'safe':
  safe <- Arguments$getLogical(safe);
  if (safe) {
    throw("Argument 'safe' was TRUE. There is no implementation available that infers MM cell indices safely from the CDF. The non-safe version assumes that stratifyBy=\"pmmm\" will return MMs in every 2nd index.");
  }

  # Argument 'ram':
  ram <- getRam(aromaSettings, ram);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  if (is.null(units)) {
    units <- seq_len(nbrOfUnits(cdf));
  } else {
    units <- Arguments$getIndices(units, max=nbrOfUnits(cdf));
  }

  CHUNK.SIZE <- as.integer(ram*100e3);
  nbrOfChunks <- ceiling(length(units) / CHUNK.SIZE);
  chunk <- 1;
  while (length(units) > 0) {
    verbose && enter(verbose, sprintf("Chunk #%d of %d", chunk, nbrOfChunks));

    uu <- 1:min(length(units), CHUNK.SIZE);
    unitsChunk <- units[uu];
    verbose && cat(verbose, "Units in this chunk:");
    verbose && str(verbose, unitsChunk);

    # Read cell indices from CDF
    verbose && enter(verbose, "Querying CDF for cell indices");
    cells <- getCellIndices(cdf, units=unitsChunk, stratifyBy="pmmm",
                                               unlist=TRUE, useNames=FALSE);
    # Not needed anymore
    unitsChunk <- NULL;
    # Sanity check
    if (length(cells) %% 2 != 0) {
      throw("Expected an even number of cell indices: ", length(cells));
    }
    verbose && exit(verbose);

    # Create (PM,MM) pairs
    cells <- matrix(cells, nrow=2);
    rownames(cells) <- c("pm", "mm");
    verbose && cat(verbose, "(PM,MM) cell indices:");
    verbose && str(verbose, cells);

    # Read all PM sequences
    verbose && cat(verbose, "PM cell indices:");
    verbose && str(verbose, cells["pm",]);

    seqs <- readSequenceMatrix(this, cells=cells["pm",], what="raw",
                                                verbose=less(verbose, 5));
    verbose && cat(verbose, "PM sequences:");
    verbose && str(verbose, seqs);

    # Keep only (PM,MM) pairs for which we know the PM sequence
    keep <- (seqs[,1] != as.raw(0));
    keep <- which(keep);
    verbose && printf(verbose, "Keeping %d of %d (%.1f%%) non-missing PM sequences\n", length(keep), nrow(seqs), 100*length(keep)/nrow(seqs));

    seqs <- seqs[keep,, drop=FALSE];
    cells <- cells[,keep, drop=FALSE];

    strands <- readTargetStrands(this, cells=cells["pm",], what="raw",
                                                verbose=less(verbose, 5));
    verbose && cat(verbose, "PM strands:");
    verbose && str(verbose, strands);

    # Keep only MM cell indices
    cells <- cells["mm",];
    # Not needed anymore
    keep <- NULL;
    verbose && cat(verbose, "MM cell indices:");
    verbose && str(verbose, cells);

    # Update target strands for MM cells
    verbose && cat(verbose, "MM strands (same as PM strands):");
    verbose && str(verbose, strands);
    updateTargetStrands(this, cells=cells, strands=strands, verbose=less(verbose, 5));
    # Not needed anymore
    strands <- NULL;

    # Complement 13th base to get MM sequences
    seqs[,13] <- dnaComplement(seqs[,13]);

    # Build MM sequences
    verbose && cat(verbose, "MM sequences:");
    verbose && str(verbose, seqs);

    updateSequenceMatrix(this, cells=cells, seqs=seqs, verbose=less(verbose, 5));
    # Not needed anymore
    cells <- seqs <- NULL;

    # Next chunk of units
    units <- units[-uu];
    # Not needed anymore
    uu <- NULL;

    gc <- gc();
    verbose && print(verbose, gc);

    chunk <- chunk + 1;
    verbose && exit(verbose);
  } # while()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update footer with MM infer information
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Updating file footer:");
  footer <- readFooter(this);
  footer$updates <- list(
    srcForMMs = list(
      methodName="inferMmFromPm",
      cdf=list(
        filename = getFilename(cdf),
        filesize = getFileSize(cdf),
        checksum = getChecksum(cdf)
      )
    )
  );
  writeFooter(this, footer);
  verbose && exit(verbose);
})



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# END: Affymetrix specific
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


############################################################################
# HISTORY:
# 2011-02-28
# o STANDARDIZATION: Now the default output path for all allocateFromCdf()
#   is annotationData/chipTypes/<chipType>/.  Before it was the same
#   directory as the CDF, which may for instance have been in a deeper
#   subdirectory, or recently also in a sibling root path.
# 2010-03-14 [HB]
# o BUG FIX: allocateFromCdf() of AromaCellPositionFile would drop all
#   but the first tag.
# 2009-09-27
# o Added argument 'keepSequenceLengths' to importFromAffymetrixProbeTabFile()
#   for AromaCellSequenceFile so that one can drop sequences of incorrect
#   lengths, cf. HuGene-1_0-st-v1.probe.tab.
# o Added argument 'onDuplicates' to importFromAffymetrixProbeTabFile() for
#   AromaCellSequenceFile.  If "error" ("warning"), an exception (warning)
#   is generated whenever duplicated cell indices are detected in the probe
#   tab file.  If "ignore", they are ignored, which means that the last
#   duplicated probe sequence will be what is finally imported.
# 2008-07-10
# o Renamed inferMmFromPmSequences() to inferMmFromPm(), because now it also
#   infers target strandedness.
# 2008-07-09
# o Created from AromaProbeSequenceTextFile.AFFX.R.
############################################################################
