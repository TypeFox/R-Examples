# @RdocClass "AromaCellSequenceFile"
#
# @title "A binary file holding cell (probe/feature) sequences of equal lengths"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \details{
#   Note that this class does \emph{not} assume a rectangular chip layout.
#   In other words, there is no concept of mapping a \emph{spatial}
#   location on the array to a cell index and vice versa.
#   The reason for this to be able to use this class also for
#   non-rectangular chip types.
# }
#
# @author

setConstructorS3("AromaCellSequenceFile", function(...) {
  extend(AromaCellTabularBinaryFile(...), "AromaCellSequenceFile");
})

setMethodS3("getFilenameExtension", "AromaCellSequenceFile", function(static, ...) {
  "acs";
}, static=TRUE)

setMethodS3("getDefaultExtension", "AromaCellSequenceFile", function(static, ...) {
  "acs";
}, static=TRUE, protected=TRUE);


setMethodS3("getExtensionPattern", "AromaCellSequenceFile", function(static, ...) {
  "[.](acs)$";
}, static=TRUE, protected=TRUE)


setMethodS3("getDefaultColumnNames", "AromaCellSequenceFile", function(this, ...) {
  c(sprintf("b%02d", seq_len(getProbeLength(this))), "targetStrand");
}, protected=TRUE)


setMethodS3("getProbeLength", "AromaCellSequenceFile", function(this, ...) {
  as.integer(nbrOfColumns(this, ...) - 1L);
})


setMethodS3("readSequenceMatrix", "AromaCellSequenceFile", function(this, cells=NULL, positions=seq_len(getProbeLength(this)), drop=FALSE, what=c("character", "raw"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'cells':
  nbrOfCells <- nbrOfCells(this);
  if (!is.null(cells)) {
    cells <- Arguments$getIndices(cells, max=nbrOfCells);
    nbrOfCells <- length(cells);
  }

  # Argument 'positions':
  positions <- Arguments$getIndices(positions, max=getProbeLength(this));

  # Argument 'what':
  what <- match.arg(what);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Reading sequence matrix");

  # Read data
  verbose && enter(verbose, "Reading data frame");
  res <- readDataFrame(this, rows=cells, columns=positions, verbose=less(verbose, 5));
  verbose && exit(verbose);

  # The raw to character map
  map <- as.raw(0:4);
  names(map) <- c(NA, "A", "C", "G", "T");

  # Flatten (data frame)
  dim <- dim(res);
  res <- unlist(res, use.names=FALSE);

  # Coerce to character strings?
  if (what == "character") {
    verbose && enter(verbose, "Coerce to a character matrix");
    res <- as.integer(res);
    res <- res + 1L;
    res <- names(map)[res];
    verbose && exit(verbose);
  }

  # Coerce to matrix
  dim(res) <- dim;

  # Drop singleton dimensions?
  if (drop) {
    res <- drop(res);
  }

  # Add 'map' attribute
  attr(res, "map") <- map;

  verbose && exit(verbose);

  res;
})


setMethodS3("readPairSequenceMatrix", "AromaCellSequenceFile", function(this, ...) {
  readNeighborSequenceMatrix(this, nbrOfNeighbors=2L, ...);
})


setMethodS3("readNeighborSequenceMatrix", "AromaCellSequenceFile", function(this, ..., nbrOfNeighbors, drop=FALSE, what=c("character", "raw", "integer", "double"), useNames=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'what':
  what <- match.arg(what);
  if (what == "character") {
    if (!useNames) {
      throw("Argument 'useNames' must be TRUE if 'what' is \"character\": FALSE");
    }
  }

  # Argument 'nbrOfNeighbors':
  nbrOfNeighbors <- Arguments$getInteger(nbrOfNeighbors, range=c(1, getProbeLength(this)));

  indexWhat <- what;
  if (what == "raw") {
    if (nbrOfNeighbors > 4) {
      throw("Cannot group nucleotides in groups larger than 4 when 'what==\"raw\"': ", nbrOfNeighbors);
    }
  } else if (what == "integer") {
    if (nbrOfNeighbors > 15) {
      throw("Cannot group nucleotides in groups larger than 15 when 'what==\"integer\"': ", nbrOfNeighbors);
    }
  } else if (what == "double") {
  } else if (what == "character") {
    if (nbrOfNeighbors <= 4) {
      indexWhat <- "raw";
    } else if (nbrOfNeighbors <= 15) {
      indexWhat <- "integer";
    } else if (nbrOfNeighbors <= 30) {
      indexWhat <- "double";
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read sequences
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading sequence matrix");
  seqWhat <- "raw";
  if (nbrOfNeighbors == 1)
    seqWhat <- what;
  seqs <- readSequenceMatrix(this, ..., what=seqWhat, verbose=less(verbose, 5));
  verbose && exit(verbose);

  # Nothing more to do?
  if (nbrOfNeighbors == 1) {
    return(seqs);
  }

  if (useNames) {
    verbose && enter(verbose, "Generating neighbor names");
    nbrOfNames <- 4^nbrOfNeighbors;
    verbose && cat(verbose, "Number of names: ", nbrOfNames);
    if (nbrOfNames > 20e6) {
      throw("Safety stop. Too many names (use a smaller 'nbrOfNeighbors'): ", nbrOfNames);
    }

    # Identify nucleotides
    map <- attr(seqs, "map");
    names <- names(map);
    names <- names[!is.na(names)];

    neighborNames <- names;
    for (kk in seq(from=2, to=nbrOfNeighbors)) {
      neighborNames <- outer(neighborNames, names, FUN=paste, sep="");
    }
    neighborNames <- as.vector(neighborNames);
    neighborNames <- sort(neighborNames);
    verbose && cat(verbose, "Names:");
    verbose && str(verbose, neighborNames);
    verbose && exit(verbose);
  } # if (useNames)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Build neighbored sequences
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Allocating return matrix");
  zeroValue <- 0; storage.mode(zeroValue) <- indexWhat;
  res <- matrix(zeroValue, nrow=nrow(seqs), ncol=ncol(seqs)-nbrOfNeighbors+1);
  verbose && str(verbose, res);
  verbose && exit(verbose);

  verbose && enter(verbose, "Identifying non-missing sequences");
  # Only need to calculate non-missing sequences
  rr <- which(seqs[,1] != as.raw(0));
  verbose && exit(verbose);


  if (length(rr) > 0) {
    verbose && enter(verbose, "Mapping sequences to neighbor-group sequences");

    calcWhat <- indexWhat;
    if (indexWhat == "raw")
      calcWhat <- "integer";
    zeroValue <- 1;  storage.mode(zeroValue) <- calcWhat;
    oneValue <- 1;  storage.mode(oneValue) <- calcWhat;

    basis <- rev(4^(1:nbrOfNeighbors-1L));
    storage.mode(basis) <- calcWhat;
    verbose && cat(verbose, "Basis:");
    verbose && print(verbose, basis);

    for (cc in seq_len(ncol(res))) {
      verbose && enter(verbose, sprintf("Position #%d of %d", cc, ncol(res)));

      values <- seqs[rr,cc];
      storage.mode(values) <- calcWhat;
      values <- values - oneValue;
      values <- basis[1] * values;
      neighbors <- values;
      # Not needed anymore
      values <- NULL;

      for (tt in 2:nbrOfNeighbors) {
        values <- seqs[rr,cc+tt-1];
        storage.mode(values) <- calcWhat;
        if (tt < nbrOfNeighbors) {
          values <- values - oneValue;
          values <- basis[tt] * values;
        }
        neighbors <- neighbors + values;
        # Not needed anymore
        values <- NULL;
      }

      storage.mode(neighbors) <- indexWhat;

      res[rr,cc] <- neighbors;
      # Not needed anymore
      neighbors <- NULL;

      verbose && exit(verbose);
    } # for (cc ...)

    verbose && exit(verbose);
  }

  if (useNames) {
    map <- 0:length(neighborNames);
    storage.mode(map) <- indexWhat;
    names(map) <- c(NA, neighborNames);
  }

  # Coerce to character strings?
  if (what == "character") {
    dim <- dim(res);
    storage.mode(res) <- calcWhat;
    res <- res + oneValue;
    res <- names(map)[res];
    dim(res) <- dim;
  }

  # Drop singleton dimensions?
  if (drop) {
    res <- drop(res);
  }

  if (useNames) {
    attr(res, "map") <- map;
  }

  res;
}, protected=TRUE)



setMethodS3("readSequences", "AromaCellSequenceFile", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Reading sequences as strings");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read raw sequence matrix
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading raw sequence matrix");
  res <- readSequenceMatrix(this, ..., what="raw");
  map <- attr(res, "map");
  verbose && str(verbose, res);
  verbose && exit(verbose);

  nbrOfCells <- nrow(res);
  nbrOfPositions <- ncol(res);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup and allocation
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Allocating sequence strings");
  # Identify non-missing sequences
  idxs <- (res[,1] != as.raw(0));
  idxs <- which(idxs);

  # Keep only those
  res <- res[idxs,,drop=FALSE];

  # Allocate return vector with missing values set
  naValue <- as.character(NA);
  seqs <- rep(naValue, times=nbrOfCells);
  seqs[idxs] <- "";  # Redo non-missing
  verbose && str(verbose, seqs);
  verbose && exit(verbose);

  # Nothing more to do?
  if (nrow(res) == 0) {
    verbose && cat(verbose, "Nothing more to do.");
    verbose && exit(verbose);
    return(seqs);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Remap
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Remapping raw values to character values");
  verbose && cat(verbose, "(raw, character) map:");
  verbose && print(verbose, map);
  values <- paste(names(map)[-1], collapse="");
  values <- charToRaw(values);
  values <- c(as.raw(0), values);

  verbose && cat(verbose, "Raw (from, to) map:");
  names(values) <- names(map);
  verbose && print(verbose, values);

  for (kk in seq_along(map)) {
    verbose && enter(verbose, sprintf("Value #%d of %d", kk, length(map)));
    verbose && cat(verbose, "Translation: 0x", map[kk], " -> 0x", values[kk]);
    idxsT <- (res == map[kk]);
    idxsT <- which(idxsT);
    verbose && cat(verbose, "Number of occurances: ", length(idxsT));
    res[idxsT] <- values[kk];
    verbose && exit(verbose);
  }
  # Not needed anymore
  idxsT <- NULL;
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Build sequence strings
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Building sequence strings");
  while (ncol(res) > 0) {
    bases <- rawToChar(res[,1], multiple=TRUE);
    seqs[idxs] <- paste(seqs[idxs], bases, sep="");
    res <- res[,-1,drop=FALSE];
  }
  verbose && exit(verbose);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && exit(verbose);

  seqs;
})


setMethodS3("readTargetStrands", "AromaCellSequenceFile", function(this, cells=NULL, what=c("character", "raw"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'cells':
  nbrOfCells <- nbrOfCells(this);
  if (!is.null(cells)) {
    cells <- Arguments$getIndices(cells, max=nbrOfCells);
    nbrOfCells <- length(cells);
  }

  # Argument 'what':
  what <- match.arg(what);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  verbose && enter(verbose, "Reading target strands");

  # Read data
  verbose && enter(verbose, "Reading data frame");
  column <- which("targetStrand" == getColumnNames(this));
  res <- readDataFrame(this, rows=cells, columns=column, drop=TRUE, verbose=less(verbose, 5));
  verbose && str(verbose, res);
  verbose && exit(verbose);

  # The raw to character map
  map <- as.raw(0:2);
  names(map) <- c(NA, "+", "-");

  # Coerce to character strings?
  if (what == "character") {
    verbose && enter(verbose, "Coerce to a character matrix");
    res <- as.integer(res);
    res <- res + 1L;
    res <- names(map)[res];
    verbose && exit(verbose);
  }

  # Add 'map' attribute
  attr(res, "map") <- map;

  verbose && exit(verbose);

  res;
})



setMethodS3("updateTargetStrands", "AromaCellSequenceFile", function(this, cells=NULL, strands, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'cells':
  nbrOfCells <- nbrOfCells(this);
  if (!is.null(cells)) {
    cells <- Arguments$getIndices(cells, max=nbrOfCells);
    nbrOfCells <- length(cells);
  }

  # Argument 'seqs':
  nbrOfStrands <- length(strands);
  if (nbrOfStrands != nbrOfCells) {
    throw("The number target strands in argument 'strands' does not match the number of cells specified: ", nbrOfStrands, " != ", nbrOfCells);
  }


  what <- mode(strands);
  if (what == "character") {
    fKeys <- c("+", "f", "forward", "sense");
    rKeys <- c("-", "r", "reverse", "antisense");
    knownKeys <- c(NA, fKeys, rKeys);
    strands <- tolower(strands);
    if (any(!strands %in% knownKeys)) {
      missing <- strands[(!strands %in% knownKeys)];
      throw("Argument 'strands' contains unknown values: ",
                                     paste(head(missing), collapse=", "));
    }
  } else if (what == "raw") {
  } else {
    throw("Argument 'strands' is of unknown type: ", mode(strands));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Optimize
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Remove duplicated 'cells'
  keep <- which(!duplicated(cells));
  cells <- cells[keep];
  strands <- strands[keep];
  # Not needed anymore
  keep <- NULL;

  # Order by 'cells'
  srt <- sort(cells, method="quick", index.return=TRUE);
  o <- srt$ix;
  cells <- srt$x;
  # Not needed anymore
  srt <- NULL;
  strands <- strands[o];
  # Not needed anymore
  o <- NULL;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Coerce to raw sequence matrix
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  map <- as.raw(0:2);
  names <- c(NA, "+", "-");
  names(map) <- names;

  if (what == "character") {
    # Coerce to a raw vector
    rawStrands <- rep(as.raw(0), times=length(strands));
    rawStrands[(strands %in% fKeys)] <- map["+"];
    rawStrands[(strands %in% rKeys)] <- map["-"];
    strands <- rawStrands;
    # Not needed anymore
    rawStrands <- NULL;
  }

  lastColumn <- nbrOfColumns(this);
  this[cells,lastColumn] <- strands;
})



setMethodS3("updateSequenceMatrix", "AromaCellSequenceFile", function(this, cells=NULL, positions=seq_len(getProbeLength(this)), seqs, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'cells':
  nbrOfCells <- nbrOfCells(this);
  if (!is.null(cells)) {
    cells <- Arguments$getIndices(cells, max=nbrOfCells);
    nbrOfCells <- length(cells);
  }

  # Argument 'positions':
  nbrOfPositions <- getProbeLength(this);
  if (is.null(positions)) {
    positions <- seq_len(nbrOfPositions);
  } else {
    positions <- Arguments$getIndices(positions, max=nbrOfPositions);
  }

  # Argument 'seqs':
  nbrOfSeqs <- nrow(seqs);
  if (nbrOfSeqs != nbrOfCells) {
    throw("The number sequences in argument 'seqs' does not match the number of cells specified: ", nbrOfSeqs, " != ", nbrOfCells);
  }
  if (ncol(seqs) != nbrOfPositions) {
    throw("The number nucleotides in argument 'seqs' does not match the number of positions specified: ", ncol(seqs), " != ", nbrOfPositions);
  }


  what <- mode(seqs);
  if (what == "character") {
  } else if (what == "raw") {
  } else {
    throw("Argument 'seqs' is of unknown type: ", mode(seqs));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Optimize
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Remove duplicated 'cells'
  keep <- which(!duplicated(cells));
  cells <- cells[keep];
  seqs <- seqs[keep,,drop=FALSE];
  # Not needed anymore
  keep <- NULL;

  # Order by 'cells'
  srt <- sort(cells, method="quick", index.return=TRUE);
  o <- srt$ix;
  cells <- srt$x;
  # Not needed anymore
  srt <- NULL;
  seqs <- seqs[o,,drop=FALSE];
  # Not needed anymore
  o <- NULL;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Coerce to raw sequence matrix
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  map <- as.raw(0:4);
  names <- c(NA, "A", "C", "G", "T");
  names(map) <- names;

  if (what == "character") {
    dim <- dim(seqs);

    # Coerce to a raw vector
    seqs <- paste(seqs, collapse="");
    seqs <- charToRaw(seqs);

    # Remap
    for (kk in seq_along(map)) {
      # Source value
      value <- names(map)[kk];

      # Identify nucleotides with this value
      if (is.na(value)) {
        idxs <- is.na(seqs);
      } else {
        idxs <- (seqs == value);
      }
      idxs <- which(idxs);

      # Update their values
      seqs[idxs] <- map[kk];
    }
    # Not needed anymore
    idxs <- NULL;

    dim(seqs) <- dim;
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update data file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (kk in seq_len(ncol(seqs))) {
    pp <- positions[kk];
    this[cells,pp] <- seqs[,kk];
  }
  # Not needed anymore
  seqs <- NULL;

  invisible(cells);
})



setMethodS3("updateSequences", "AromaCellSequenceFile", function(this, ..., seqs, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Argument 'seqs':
  if (!is.character(seqs)) {
    throw("Argument 'seqs' must be a character vector: ", mode(seqs));
  }
  probeLength <- nchar(seqs);
  probeLength <- unique(probeLength);
  if (length(probeLength) != 1) {
    throw("Cannot write probe sequences. Sequences of varying lengths detected: ", paste(head(probeLength), collapse=", "));
  }


  verbose && enter(verbose, "Updating file with sequence strings");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Coerce to a raw sequence matrix
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Coercing to raw sequence matrix");
  nbrOfSequences <- length(seqs);
  seqs <- paste(seqs, collapse="");
  seqs <- charToRaw(seqs);

  # Remap
  from <- charToRaw(" ACGT");
  to <- as.raw(0:4);
  for (kk in seq_len(5L)) {
    idxs <- which(seqs == from[kk]);
    seqs[idxs] <- to[kk];
  }

  seqs <- matrix(seqs, nrow=nbrOfSequences, ncol=probeLength, byrow=TRUE);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Update file with raw sequence matrix");
  verbose && str(verbose, seqs);
  res <- updateSequenceMatrix(this, ..., seqs=seqs, verbose=less(verbose, 5));
  verbose && exit(verbose);

  verbose && exit(verbose);

  invisible(res);
})




setMethodS3("isMissing", "AromaCellSequenceFile", function(this, ...) {
  res <- readSequenceMatrix(this, ..., positions=1, what="raw", drop=TRUE);
  res <- (res == as.raw(0));
  res;
}, protected=TRUE)


setMethodS3("countBases", "AromaCellSequenceFile", function(this, bases=c("A", "C", "G", "T"), drop=FALSE, mode=c("integer", "raw"), ...) {
  # Argument 'mode':
  mode <- match.arg(mode);


  # Mode-specific values
  if (mode == "raw") {
    zeroValue <- as.raw(0);
    naValue <- as.raw(255);
  } else {
    zeroValue <- 0L;
    naValue <- as.integer(NA);
  }

  # Tabular nucleotides
  counts <- countBasesInternal(this, mode=mode, ...);

  # Identify missing sequences
  isMissing <- which(counts[,1] != zeroValue);

  # Keep only bases of interest
  counts <- counts[,bases,drop=FALSE];

  # Set missing values
  counts[isMissing,] <- naValue;

  # Drop singleton dimensions?
  if (drop) {
    if (mode == "raw") {
      # To be verified. /HB 2008-07-20
      dim <- dim(counts);
      dim <- dim[dim > 1];
      dim(counts) <- dim;
    } else {
      counts <- drop(counts);
    }
  }

  counts;
})

setMethodS3("countBasesInternal", "AromaCellSequenceFile", function(this, cells=NULL, positions=seq_len(getProbeLength(this)), mode=c("integer", "raw"), ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'cells':
  nbrOfCells <- nbrOfCells(this);
  if (!is.null(cells)) {
    cells <- Arguments$getIndices(cells, max=nbrOfCells);
    nbrOfCells <- length(cells);
  }

  # Argument 'positions':
  positions <- Arguments$getIndices(positions, max=getProbeLength(this));

  # Argument 'mode':
  mode <- match.arg(mode);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  verbose && enter(verbose, "Counting occurances of each nucleotide");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chipType <- getChipType(this);
  key <- list(method="countBases", class=class(this)[1],
              chipType=chipType, fullname=getFullName(this),
              cells=cells, positions=positions, mode=mode);
  if (getOption(aromaSettings, "devel/useCacheKeyInterface", FALSE)) {
    key <- getCacheKey(this, method="countBases", chipType=chipType,
                              cells=cells, positions=positions, mode=mode);
  }
  dirs <- c("aroma.affymetrix", chipType);
  if (!force) {
    counts <- loadCache(key=key, dirs=dirs);
    if (!is.null(counts)) {
      verbose && cat(verbose, "Cached results found.");
      return(counts);
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Optimize reading order?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(cells)) {
    cells <- seq_len(nbrOfCells);
    reorder <- FALSE;
  } else {
    srt <- sort(cells, method="quick", index.return=TRUE);
    o <- srt$ix;
    cells <- srt$x;
    # Not needed anymore
    srt <- NULL;
    reorder <- TRUE;
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sum over positions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (mode == "raw") {
    zeroValue <- as.raw(0);
  } else {
    zeroValue <- 0L;
  }
  counts <- matrix(zeroValue, nrow=nbrOfCells, ncol=5);

  map <- NULL;
  for (kk in seq_along(positions)) {
    pp <- positions[kk];
    verbose && enter(verbose, sprintf("Position #%d (%d) of %d", kk, pp, length(positions)));
    seqs <- readSequenceMatrix(this, cells=cells, positions=pp, what="raw", drop=TRUE, verbose=less(verbose, 20));

    if (is.null(map)) {
      map <- attr(seqs, "map");
      colnames(counts) <- names(map);
    }

    # Add to counts
    verbose && enter(verbose, "Summing counts");
    for (bb in seq_len(ncol(counts))) {
      verbose && enter(verbose, sprintf("Nucleotide #%d of %d", bb, ncol(counts)));
      idxs <- (seqs == map[bb]);
      idxs <- which(idxs);

      verbose && cat(verbose, "Increment:");
      countsBB <- counts[idxs,bb];
      if (mode == "raw") {
        countsBB <- as.integer(countsBB);
        countsBB <- countsBB + 1L;
        countsBB <- as.raw(countsBB);
      } else {
        countsBB <- countsBB + 1L;
      }
      counts[idxs,bb] <- countsBB;
      # Not needed anymore
      idxs <- countsBB <- NULL;
      verbose && exit(verbose);
    } # for (bb ...)
    # Not needed anymore
    seqs <- NULL;
    verbose && exit(verbose);

    verbose && exit(verbose);
  } # for (kk ...)

  # Reorder?
  if (reorder) {
    o <- order(o);
    counts <- counts[o,,drop=FALSE];
    # Not needed anymore
    o <- NULL;
    gc <- gc();
    verbose && print(verbose, gc);
  }

  # Cache results
  saveCache(key=key, dirs=dirs, counts);

  verbose && exit(verbose);

  counts;
}, protected=TRUE)


setMethodS3("allocate", "AromaCellSequenceFile", function(static, ..., nbrOfCells, platform, chipType, footer=list()) {
  # Argument 'nbrOfCells':
  nbrOfCells <- Arguments$getInteger(nbrOfCells, range=c(1, 1000e6));

  # Argument 'platform':
  platform <- Arguments$getCharacter(platform, length=c(1,1));

  # Argument 'chipType':
  chipType <- Arguments$getCharacter(chipType, length=c(1,1));

  # Argument 'footer':
  if (is.null(footer)) {
  } else if (!is.list(footer)) {
    throw("Argument 'footer' must be NULL or a list: ", class(footer)[1]);
  }

  footer <- c(
    list(
      createdOn=format(Sys.time(), "%Y%m%d %H:%M:%S", usetz=TRUE),
      platform=platform,
      chipType=chipType
    ),
    footer
  );

  probeLengths <- 25L;
  nbrOfColumns <- probeLengths + 1L;

  NextMethod("allocate", nbrOfRows=nbrOfCells, types=rep("raw", times=nbrOfColumns), sizes=rep(1L, times=nbrOfColumns), footer=footer);
}, static=TRUE, protected=TRUE)



############################################################################
# HISTORY:
# 2008-08-04
# o Added support for argument 'positions' to countBases().
# o BUG FIX: readSequences() translated raw values to incorrect nucleotides.
# 2008-07-20
# o Now countBases() returns "raw" counts, if argument 'mode="raw"'.
# 2008-07-10
# o Added read- and updateTargetStrands().
# o Now readNeighborSequenceMatrix() takes what="integer" and "double" too.
# o Update updateSequences().
# o Made the decision that for sequences in character mode, missing
#   sequences/nucleotides are always represented as NA.  A sequence is
#   defined to be missing if the 1st nucleotide is missing.
#   For "raw" sequences, the value zero is representing "missing".
# 2008-07-09
# o Added updateSequenceMatrix().
# o Created from AromaUgpFile.R.
############################################################################
