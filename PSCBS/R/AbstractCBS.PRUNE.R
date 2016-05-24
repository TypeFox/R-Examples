setMethodS3("seqOfSegmentsByDP", "AbstractCBS", function(fit, by, shift=+100, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'by':
  by <- Arguments$getCharacters(by);
  data <- getLocusData(fit);
  fields <- colnames(data);
  missing <- fields[!is.element(by, fields)];
  if (length(missing) > 0) {
    throw("Argument 'by' specifies one or more non-existing locus data fields: ", paste(missing, collapse=", "));
  }

  # Argument 'shift':
  shift <- Arguments$getNumeric(shift);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Identifying optimal sets of segments via dynamic programming");

  # Assert that known segments was used
  knownSegments <- fit$params$knownSegments;
  if (nrow(knownSegments) == 0L) {
    chromosome <- getChromosomes(fit);
    knownSegments <- data.frame(chromosome=chromosome, start=-Inf, end=+Inf);
  }

  chrOffsets <- getChromosomeOffsets(fit);
  # Sanity check
  stopifnot(all(is.finite(chrOffsets)));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Shift every other non-empty region
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Shifting TCN levels for every second segment");

  segPrefix <- getSegmentTrackPrefixes(fit)[1L];
  segKeys <- toCamelCase(paste(segPrefix, c("start", "end")));
  segRowsKey <- toCamelCase(paste(segPrefix, "seg rows"));

  verbose && enter(verbose, "Split up into non-empty independent regions");
  chromosomes <- getChromosomes(fit);

  fitList <- list();
  for (jj in seq(along=chromosomes)) {
    chr <- chromosomes[jj];
    verbose && enter(verbose, sprintf("Chromosome #%d ('%s') of %d", jj, chr, length(chromosomes)));

    # Subset segmentation results on this chromosome
    fitJJ <- extractChromosome(fit, chr);
    verbose && cat(verbose, "Number of loci on chromosome: ", nbrOfLoci(fitJJ));

    # Nothing to do and nothing to add?
    if (nbrOfLoci(fitJJ) == 0L) {
      verbose && exit(verbose);
      next;
    }

    # Find known segments on this chromosome
    knownSegmentsJJ <- subset(knownSegments, chromosome == chr);
    verbose && cat(verbose, "Known segments on chromosome:");
    verbose && print(verbose, knownSegmentsJJ);

    # Nothing to do?
    if (nrow(knownSegmentsJJ) == 0L) {
      fitList <- append(fitList, list(fitJJ));
      verbose && exit(verbose);
      next;
    }

    # Get the segments on this chromosome
    segsJJ <- getSegments(fitJJ);

    # Extract the individual known segments on this chromosome
    fitListJJ <- list();
    for (kk in seq(length=nrow(knownSegmentsJJ))) {
      verbose && enter(verbose, sprintf("Known segment #%d of %d", kk, nrow(knownSegmentsJJ)));
      seg <- knownSegmentsJJ[kk,];
      verbose && print(verbose, seg);

      start <- seg$start;
      end <- seg$end;
      idxStart <- min(which(segsJJ[[segKeys[1]]] >= start));
      idxEnd <- max(which(segsJJ[[segKeys[2]]] <= end));
      idxs <- idxStart:idxEnd;

      # Extract the particular known segment
      fitKK <- extractSegments(fitJJ, idxs);

      # Only add if it has loci
      if (nbrOfLoci(fitKK) > 0L) {
        fitListJJ <- append(fitListJJ, list(fitKK));
      }

      fitKK <- NULL; # Not needed anymore
      verbose && exit(verbose);
    } # for (kk ...)

    # Append
    fitList <- append(fitList, fitListJJ);
    fitListJJ <- NULL; # Not needed anymore

    verbose && exit(verbose);
  } # for (jj ...)

  nbrOfRegions <- length(fitList);
  verbose && cat(verbose, "Number of independent non-empty regions: ", nbrOfRegions);
  verbose && exit(verbose);

  verbose && enter(verbose, "Shift every other region");
  for (jj in seq(from=1L, to=nbrOfRegions, by=2L)) {
    fitJJ <- fitList[[jj]];
    fitJJ <- shiftTCN(fitJJ, shift=shift);
    fitList[[jj]] <- fitJJ;
  } # for (jj ...)
  verbose && exit(verbose);

  verbose && enter(verbose, "Merge");
  fitT <- Reduce(function(a, b) append(a,b, addSplit=FALSE), fitList);
  # Sanity check
##  stopifnot(nbrOfSegments(fitT) == nbrOfSegments(fit)); # Not true anymore
  verbose && exit(verbose);
  fitList <- NULL; # Not needed anymore

  segsT <- getSegments(fitT);
  verbose && print(verbose, tail(segsT));
  verbose && exit(verbose);

  fit <- fitT;

  # Not needed anymore
  fitT <- knownSegments <- NULL;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract signals for DP
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Extracting signals for dynamic programming");

  fit <- tileChromosomes(fit);

  # Locus-level signals
  data <- getLocusData(fit);
  Y <- as.matrix(data[,by,drop=FALSE]);
  verbose && print(verbose, summary(Y));

  # "DP" change-point indices (excluding the two outer/boundary ones)
  segRows <- fit[[segRowsKey]];
  segIdxs <- seq(length=nrow(segRows)-1L);
  cpIdxs <- segRows$endRow[segIdxs];

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Dynamic-programming segmention pruning
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Dynamic programming");

  verbose && cat(verbose, "Number of \"DP\" change points: ", length(cpIdxs));
  verbose && str(verbose, cpIdxs);

  res <- seqOfSegmentsByDP(Y, candidatechangepoints=cpIdxs, ...);
  verbose && str(verbose, res);

  # Sanity checks
  jumpList <- res$jump;
  lastJump <- jumpList[[length(jumpList)]];
  stopifnot(identical(cpIdxs, as.integer(lastJump)));

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Adjustments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Excluding cases where known segments no longer correct");

  verbose && cat(verbose, "Number of independent non-empty regions: ", nbrOfRegions);

  # Drop the K first
  nbrOfCPs <- nbrOfRegions - 1L;
  if (nbrOfCPs > 0L) {
    K <- nbrOfCPs - 1L;  # Don't count the flat segmentation
    jumpList <- jumpList[-(1:K)];
  }

  verbose && str(verbose, jumpList);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract the possible sets of segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose <- less(verbose, 20);
  verbose && enter(verbose, "Converting to physical segments");

  segs <- getSegments(fit, splitters=TRUE, addGaps=FALSE);
  segs <- segs[,c("chromosome", segKeys)];
  jumpIdxList <- lapply(jumpList, FUN=function(idxs) {
    match(idxs, table=cpIdxs);
  });
  # Sanity check
  stopifnot(identical(seq(along=cpIdxs), jumpIdxList[[length(jumpIdxList)]]));

  chrs <- segs$chromosome;
  starts <- segs[[segKeys[1]]];
  ends <- segs[[segKeys[2]]];
  nsegs <- nrow(segs);
  segList <- vector("list", length=length(jumpIdxList));
  for (kk in seq(along=segList)) {
    verbose && enter(verbose, sprintf("Sequence #%d of %d", kk, length(segList)));
    idxs <- jumpIdxList[[kk]];
    verbose && cat(verbose, "Change point indices:");
    verbose && str(verbose, idxs);

    chrsKK <- chrs[idxs];
    chr <- chrsKK[1];
#    stopifnot(all(chrsKK == chr));
    chrsKK <- c(chrsKK, chrsKK[length(chrsKK)]);
    startsKK <- starts[c(1L, idxs+1L)];
    endsKK <- ends[c(idxs, nsegs)];
    verbose && cat(verbose, "Chromosomes:");
    verbose && str(verbose, chrsKK);
    verbose && cat(verbose, "Starts:");
    verbose && str(verbose, startsKK);
    verbose && cat(verbose, "Ends:");
    verbose && str(verbose, endsKK);


    offsetsKK <- chrOffsets[chrsKK];
    startsKK <- startsKK - offsetsKK;
    endsKK <- endsKK - offsetsKK;
    segsKK <- data.frame(chromosome=chrsKK, start=startsKK, end=endsKK);
    verbose && print(verbose, tail(segsKK));
    segList[[kk]] <- segsKK;

    verbose && exit(verbose);
  } # for (kk ...)

  verbose && exit(verbose);
  verbose <- more(verbose, 20);


  # Sanity check
  lastSegs <- segList[[length(segList)]];
#  stopifnot(identical(lastSegs, segs));
  verbose && str(verbose, segList);

  nbrOfCPsSeq <- sapply(jumpList, FUN=length);
  verbose && cat(verbose, "Sequence of number of \"DP\" change points:");
  verbose && print(verbose, nbrOfCPsSeq);

  nbrOfSegsSeq <- sapply(segList, FUN=nrow);
  verbose && cat(verbose, "Sequence of number of segments:");
  verbose && print(verbose, nbrOfSegsSeq);

  nbrOfChangePointsSeq <- nbrOfSegsSeq - nbrOfRegions;
  verbose && cat(verbose, "Sequence of number of \"discovered\" change points:");
  verbose && print(verbose, nbrOfChangePointsSeq);

  stopifnot(nbrOfSegsSeq == nbrOfCPsSeq + 1L);
  stopifnot(nbrOfSegsSeq[1L] == nbrOfRegions);

  segList <- lapply(segList, FUN=function(seg) {
    attr(seg, "nbrOfChangePoints") <- nrow(seg) - nbrOfRegions;
    seg;
  });

  K <- (nbrOfRegions-1L);
  modelFit <- list(
    nbrOfSegments=nbrOfSegsSeq,
    nbrOfChangePoints=nbrOfSegsSeq - nbrOfRegions,
    nbrOfRegions=nbrOfRegions,
    idxBest=res$kbest - K + 1L,
    rse=res$rse[-(1:K)],
    V=res$V[-(1:K),-(1:K),drop=FALSE],
    seqOfSegmentsByDP=res
  );
  attr(segList, "modelFit") <- modelFit;

  verbose && exit(verbose);

  segList;
}, protected=TRUE) # seqOfSegmentsByDP()



setMethodS3("seqOfSegmentsByDP", "matrix", function(
### Segmentation of a multi-dimensional signal with dynamic programming.
Y,
### A n*p signal to be segmented
candidatechangepoints = c(1:(n-1)),
### A vector of candidate positions for change-points (default=[1:n-1])
kmax = n-1,
### Maximum number of change-points to test (default : n-1)
threshold = 0.5,
### Stopping criteria. Typically chosen to be in the interval
### (0 0.5]. The smaller the threshold, the higher the tendency to keep more
### changepoints. The criteria is based on the method found in 'Picard et al
### (2005)', "A statistical approach for array CGH data analysis" (BMC
### Bioinformatics). Default=0.5
...
){
  n <- dim(Y)[1];
  p <- dim(Y)[2];
  kmaxmin <- floor(n/4);

  if (kmax == 0 || kmax > length(candidatechangepoints)) {
    kmax <- min(length(candidatechangepoints), kmaxmin);
  }

  if (kmax > kmaxmin) {
    sprintf('warning : not enough points to optimize the number of the change-points up to %s\n', kmax);
    kmax <- kmaxmin;
    sprintf('Set the maximum number of change-points to %s\n', kmax);
  }

  # Compute boundaries of the smallest intervals considered
  b <- sort(union(0, union(n, candidatechangepoints)));
  k <- length(b) - 1; # k is the number of such intervals
  # Compute the k*k matrix J such that J[i,j] for i<=j is the RSE when intervales i to j are merged
  J <- matrix(numeric(k*k), ncol=k);

  # How should NAs be handled?!?
  Yz <- Y;
  Yz[is.na(Yz)] <- 0;
  s <- rbind(rep(0, times=p), colCumsums(Yz));
  v <- c(0, cumsum(rowSums(Y^2, na.rm=TRUE)));
  for (i in 1:k) {
    for (j in i:k) {
      Istart <- b[i] + 1;
      Iend <- b[j+1];
      J[i,j] <- v[Iend+1] - v[Istart] - sum((s[Iend+1,]-s[Istart,])^2, na.rm=TRUE)/(Iend-Istart+1);
    } # for (j ...)
  } # for (i ...)

  # Dynamic programming
  V <- matrix(numeric((kmax+1)*k), ncol=k); # V[i,j] is the best RSE for segmenting intervals 1 to j with at most i-1 change points
  jump <- matrix(numeric(kmax*k), ncol=k);
  # With no change points, V[i,j] is juste the precomputed RSE for intervals 1 to j
  V[1,] <- J[1,];
  # Then we apply the recursive formula
  for (ki in 1:kmax) {
    for (j in (ki+1):k) {
      val <- min(V[ki,ki:(j-1)] + t(J[(ki+1):j, j]));
      ind <- which.min(V[ki,ki:(j-1)] + t(J[(ki+1):j, j]));
      V[ki+1,j] <- val;
      jump[ki,j] <- ind + ki-1;
    } # for (j ...)
  } # for (ki ...)

  # Optimal segmentation
  res.jump <- list();
  for (ki in 1:kmax) {
    res.jump[[ki]] <- numeric(ki);
    res.jump[[ki]][ki] <- jump[ki,k];
    if (ki != 1) {
      for (i in seq(from=(ki-1), to=1, by=-1)) {
        res.jump[[ki]][i] <- jump[i, res.jump[[ki]][i+1]];
      }
    }
  } # for (ki ...)

  # Convert back the index of the interval to the last position before the jump
  rightlimit <- b[2:length(b)];
  for (ki in 1:kmax) {
    res.jump[[ki]] <- rightlimit[res.jump[[ki]]];
  } # for (ki ...)

  # RSE as a function of number of change-points
  res.rse <- V[,k];
  # Optimal number of change points
  options(warn=-1);
  J <- log(res.rse);
  Km <- length(J);
  Jtild <- (J[Km]-J)/(J[Km]-J[1])*(Km-1)+1; # Normalize
  res.kbest <- max(which(diff(diff(Jtild)) > threshold)) + 1;
  #if((res.kbest) == -Inf) { res.kbest <- 1 };
  return(list(jump=res.jump, rse=res.rse, kbest=res.kbest, V=V));
### \item{res.jump{i}}{a i*1 vector of change-point positions for the i-th lambda value (i depends on lambda). i varies between 1 and kmax}
### \item{res.rse}{a (kmax+1)-dimensional vector of residual squared error}
### \item{res.kbest}{the number of selected change-points}
}) # seqOfSegmentsByDP()



###########################################################################/**
# @set class=AbstractCBS
# @RdocMethod pruneByDP
#
# @title "Prunes the CN profile using dynamical programming"
#
# \description{
#  @get "title" by specifying the target number of segments or alternative
#  how of many change points to drop.
# }
#
# @synopsis
#
# \arguments{
#  \item{nbrOfSegments}{An @integer specifying the number of segments after
#     pruning. If negative, the it specifies the number of change points
#     to drop.}
#  \item{...}{Optional arguments passed to @seemethod "seqOfSegmentsByDP".}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a pruned object of the same class.
# }
#
# \examples{\dontrun{
#  # Drop two segments
#  fitP <- pruneByDP(fit, nbrOfSegments=-2);
# }}
#
# @author "HB, PN"
#
# \references{
#   [1] ... \cr
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("pruneByDP", "AbstractCBS", function(fit, nbrOfSegments, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Some pre-extraction
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  knownSegments <- fit$params$knownSegments;
  if (nrow(knownSegments) == 0L) {
    chromosome <- getChromosomes(fit);
    knownSegments <- data.frame(chromosome=chromosome, start=-Inf, end=+Inf);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'nbrOfSegments':
  nbrOfSegments <- Arguments$getInteger(nbrOfSegments);
  # Specifying number of change points *to drop*?
  if (nbrOfSegments < 0L) {
    nbrOfCPsToDrop <- -nbrOfSegments;
    nbrOfSegments <- nbrOfSegments(fit, splitters=FALSE) - nbrOfCPsToDrop;
  }

  if (nbrOfSegments < nrow(knownSegments)) {
    throw("Argument 'nbrOfSegments' is less than number of \"known\" segments: ", nbrOfSegments, " < ", nrow(knownSegments));
  }
  if (nbrOfSegments > nbrOfSegments(fit, splitters=FALSE)) {
    throw("Argument 'nbrOfSegments' is greater than the number of \"found\" segments: ", nbrOfSegments, " > ", nbrOfSegments(fit, splitters=FALSE));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Prune segments by dynamical programming");

  # Nothing to do?
  if (nbrOfSegments == nbrOfSegments(fit, splitters=FALSE)) {
    verbose && cat(verbose, "No need for pruning, because the number of \"target\" segments is the same as the number of \"found\" segments: ", nbrOfSegments, " == ", nbrOfSegments(fit, splitters=FALSE));
    return(fit);
    verbose && exit(verbose);
  }

  verbose && cat(verbose, "Target number of segments: ", nbrOfSegments);

  verbose && enter(verbose, "Dynamic programming");
  segList <- seqOfSegmentsByDP(fit, ...);

  # Select the one with expected number of segments
  nbrOfCPs <- sapply(segList, FUN=nrow);
  verbose && printf(verbose, "Range of number of CPs among solutions: [%d,%d]\n", min(nbrOfCPs), max(nbrOfCPs));

  keep <- which(nbrOfCPs == nbrOfSegments);
  stopifnot(length(keep) == 1L);

  knownSegments <- segList[[keep]];
  verbose && printf(verbose, "Solution with %d segments:\n", nbrOfSegments);
  verbose && print(verbose, knownSegments);
  verbose && exit(verbose);

  verbose && enter(verbose, "Rebuilding segmentation results");
  fitDP <- resegment(fit, knownSegments=knownSegments, undoTCN=+Inf, undoDH=+Inf, verbose=less(verbose, 10));
  verbose && print(verbose, fitDP);
  verbose && exit(verbose);

  verbose && exit(verbose);

  fitDP;
}, protected=TRUE) # pruneByDP()




############################################################################
# HISTORY:
# 2012-09-21
# o Now seqOfSegmentsByDP() return a 'modelFit' attribute.
# 2012-09-20
# o BUG FIX: seqOfSegmentsByDP() for AbstractCBS would not handle empty
#   segments, which could occur if 'knownSegments' for instance included
#   centromere gaps.
# 2012-09-14
# o Added pruneByDP(..., nbrOfSegments) for AbstractCBS.
# o Added seqOfSegmentsByDP() which is a matrix version of dpseg() from
#   GFLseg v0.1.6 by Morgane Pierre-Jean and Pierre Neuvial, which in turn
#   ported it from the Matlab GPL source code by the original authors ???.
# o Generalized seqOfSegmentsByDP() to work for AbstractCBS.
# 2012-09-13
# o Added seqOfSegmentsByDP() for PairedPSCBS.
# o Created.
############################################################################
