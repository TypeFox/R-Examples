###########################################################################/**
# @set class=PairedPSCBS
# @RdocMethod bootstrapTCNandDHByRegion
#
# @title "Estimate confidence intervals of TCN and DH segment levels"
#
# \description{
#  @get "title" using bootstrap.
# }
#
# @synopsis
#
# \arguments{
#   \item{B}{A postive @integer specifying the number of bootstrap samples.}
#   \item{boot}{Alternatively, to generating \code{B} bootstrap samples,
#      this specifies a pre-generated set of bootstrap samples as
#      returned by \code{bootstrapSegmentsAndChangepoints()}.}
#   \item{...}{Additional arguments passed to \code{bootstrapSegmentsAndChangepoints()}.}
#   \item{probs}{The default quantiles to be estimated.}
#   \item{statsFcn}{A (optional) @function that estimates confidence
#      intervals given locus-level data.
#      If @NULL, the @see "stats::quantile" function is used.}
#   \item{what}{A @character @vector specifying what to bootstrap.}
#   \item{force}{If @TRUE, already existing estimates are ignored,
#      otherwise not.}
#   \item{verbose}{See @see "R.utils::Verbose".}
#   \item{.debug}{(internal) If @TRUE, additional sanity checks are
#      performed internally.}
# }
#
# \value{
#   Returns a @see "PairedPSCBS" object.
# }
#
# @author "HB"
#*/###########################################################################
setMethodS3("bootstrapTCNandDHByRegion", "PairedPSCBS", function(fit, B=1000L, boot=NULL, ..., probs=c(0.025, 0.050, 0.95, 0.975), statsFcn=NULL, what=c("segment", "changepoint"), force=FALSE, verbose=FALSE, .debug=FALSE) {
  # Settings for sanity checks
  tol <- getOption("PSCBS/sanityChecks/tolerance", 0.0005);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  summarizeSamples <- function(X, statsFcn, stats=NULL, what=c("segment", "changepoint"), ..., verbose=FALSE) {
    # Argument 'X':
    stopifnot(is.array(X));
    dim <- dim(X);
    stopifnot(length(dim) == 3L);

    # Argument 'statsFcn':
    stopifnot(is.function(statsFcn));
    statsT <- statsFcn(1);
    stopifnot(!is.null(names(statsT)));
    nbrOfStats <- length(statsT);
    statsNames <- names(statsT);
    statsT <- NULL; # Not needed anymore

    # Argument 'stats':
    if (!is.null(stats)) {
      stopifnot(is.data.frame(stats));
    }

    # Argument 'what':
    what <- match.arg(what);
    whatC <- capitalize(what);

    # Argument 'verbose':
    verbose <- Arguments$getVerbose(verbose);
    if (verbose) {
      pushState(verbose);
      on.exit(popState(verbose));
    }


    dimnames <- dimnames(X);
    fields <- dimnames[[3L]];

    verbose && enter(verbose, sprintf("Summarizing bootstrapped %s (%s) data", what, paste(sQuote(fields), collapse=", ")));

    # Allocate JxQxF matrix S
    dim[2L] <- nbrOfStats;
    dimnames[[2L]] <- statsNames;
    S <- array(NA_real_, dim=dim, dimnames=dimnames);
    verbose && str(verbose, S);

    for (kk in seq(along=fields)) {
      field <- fields[kk];
      verbose && enter(verbose, sprintf("Field #%d ('%s') of %d", kk, field, length(fields)));

      Xkk <- X[,,kk,drop=FALSE];  # An JxB matrix
      dim(Xkk) <- dim(Xkk)[-3L];
      # Sanity check
      stopifnot(is.matrix(Xkk));
      stopifnot(nrow(Xkk) == dim(X)[1L]);
      stopifnot(ncol(Xkk) == B);

      for (jj in seq(length=dim(X)[1L])) {
        verbose && enter(verbose, sprintf("%s #%d of %d", whatC, jj, dim(X)[1L]));
        Xkkjj <- Xkk[jj,,drop=TRUE]; # A vector of length B
        S[jj,,kk] <- statsFcn(Xkkjj);
        verbose && exit(verbose);
      } # for (jj ...)

      Xkk <- NULL; # Not needed anymore

      verbose && exit(verbose);
    } # for (jj ...)

    # Not needed anymore
    X <- NULL;

    verbose && cat(verbose, "Bootstrap statistics");
    verbose && str(verbose, S);

    # Reshape JxQx4 array to Jx(4*Q) matrix
    T <- wrap(S, map=list(1,NA), sep="_");
    colnames(T) <- gsub("(.*)_(.*)", "\\2_\\1", colnames(T));

    # Append as new columns to the summary table
    stats <- cbind(stats, T);

    # Drop previously estimated values
    dups <- duplicated(colnames(stats), fromLast=TRUE);
    if (any(dups)) {
      stats <- stats[,!dups, drop=FALSE];
    }

    # Not needed anymore
    T <- dups <- NULL;


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Statistical sanity checks
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (what == "segment" && B >= 100L) {
      verbose && enter(verbose, "Statistical sanity checks (iff B >= 100)");

      stopifnot(is.array(S));

      # Find extreme quantiles
      probs <- dimnames(S)[[2L]];
      verbose && printf(verbose, "Available summaries: %s\n", paste(probs, collapse=", "));
      probs <- grep("%", probs, fixed=TRUE, value=TRUE);
      S <- S[,probs,,drop=FALSE];
      probs <- gsub("%", "", probs, fixed=TRUE);
      probs <- as.double(probs) / 100;
      verbose && printf(verbose, "Available quantiles: %s\n", paste(probs, collapse=", "));
      verbose && str(verbose, S);
      # Sanity check
      stopifnot(all(is.finite(probs)));

      # Is it possible to check?
      if (any(probs < 0.10) && any(probs > 0.90)) {
        tryCatch({
          fields <- dimnames(S)[[3L]];
          for (kk in seq(along=fields)) {
            field <- fields[kk];
            verbose && enter(verbose, sprintf("Field #%d ('%s') of %d", kk, field, length(fields)));

            # Bootstrap statistics
            Skk <- S[,,kk, drop=FALSE];
            dim(Skk) <- dim(Skk)[-3L];

            # Sanity checks
            stopifnot(is.matrix(Skk));

            range <- Skk[,c(1L,ncol(Skk)),drop=FALSE];

            # Segmentation means
            key <- sprintf("%sMean", field);
            segMean <- segs[[key]];

            # Segmentation counts
            cfield <- sprintf("%sNbrOfLoci", ifelse(field == "tcn", "tcn", "dh"));
            counts <- segs[,cfield,drop=TRUE];

            if (verbose) {
              for (rr in seq_len(length(segMean))) {
                printf(verbose, "Seg %3d. mean=%g, range=[%g,%g], n=%d\n", rr, segMean[rr], range[rr,1L], range[rr,2L], counts[rr]);
              } # for (rr ...)
            }

            # Compare only segments with enough data points
            keep <- (counts > 1L);
            range <- range[keep,,drop=FALSE];
            segMean <- segMean[keep];

            # Sanity checks
            stopifnot(all(range[,2L] + tol >= range[,1L], na.rm=TRUE));
            stopifnot(all(segMean + tol >= range[,1L], na.rm=TRUE));
            stopifnot(all(segMean - tol <= range[,2L], na.rm=TRUE));

            verbose && exit(verbose);
          } # for (kk ...)
        }, error = function(ex) {
          # If an error, display the data, then throw the exception
          verbose && cat(verbose, "Tolerance (option 'PSCBS/sanityChecks/tolerance'): ", tol);
          verbose && print(verbose, segs);
          throw(ex);
        })
      } else {
        verbose && cat(verbose, "Skipping. Not enough quantiles: ",
                                 paste(dimnames(S)[[2L]], collapse=", "));
      }

      verbose && exit(verbose);
    } # if (B >= 100L)


    verbose && exit(verbose);

    stats;
  } # summarizeSamples()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'B':
  B <- Arguments$getInteger(B, range=c(1,Inf));

  # Argument 'probs':
  probs <- Arguments$getNumerics(probs, range=c(0,1));
  # Always estimate the default quantiles
  probs0 <- eval(formals(bootstrapTCNandDHByRegion.PairedPSCBS)$probs);
  probs <- unique(sort(c(probs, probs0)));

  # Argument 'statsFcn':
  if (is.null(statsFcn)) {
    statsFcn <- function(x) quantile(x, probs=probs, na.rm=TRUE);
  }

  # Argument 'what':
  what <- unique(match.arg(what, several.ok=TRUE));

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Argument '.debug':
  .debug <- Arguments$getLogical(.debug);



  verbose && enter(verbose, "Resample (TCN,DH) signals and re-estimate summaries for ", paste(what, collapse=" & "));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract existing estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.element("segment", what)) {
    segs <- getSegments(fit);
  }
  if (is.element("changepoint", what)) {
    cps <- getChangePoints(fit);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Already done?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  stats <- statsFcn(1);
  stopifnot(!is.null(names(stats)));
  nbrOfStats <- length(stats);
  statsNames <- names(stats);

  if (is.element("segment", what)) {
    tcnStatsNames <- sprintf("tcn_%s", names(stats));
    dhStatsNames <- sprintf("dh_%s", names(stats));
    c1StatsNames <- sprintf("c1_%s", names(stats));
    c2StatsNames <- sprintf("c2_%s", names(stats));
    allStatsNames <- c(tcnStatsNames, dhStatsNames, c1StatsNames, c2StatsNames);
    isDone <- is.element(allStatsNames, names(segs));
    names(isDone) <- allStatsNames;
    verbose && cat(verbose, "Already done?");
    verbose && print(verbose, isDone);

    # Not needed anymore
    allStatsNames <- tcnStatsNames <- dhStatsNames <-
                     c1StatsNames <- c2StatsNames <- NULL;

    if (!force && all(isDone)) {
      verbose && cat(verbose, "Already done. Skipping.");
      verbose && exit(verbose);
      return(fit);
    }
  }


  # The object to be returned
  fitB <- fit;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Bootstrap (TCN,DH,C1,C2) segment mean levels
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(boot)) {
    boot <- bootstrapSegmentsAndChangepoints(fit, B=B, ...,
                         force=force, .debug=.debug, verbose=verbose);
  } else {
    B <- dim(boot$segments)[2L];
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Summarizing segment (TCN,DH,C1,C2) mean levels
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.element("segment", what)) {
    segs <- summarizeSamples(boot$segments, statsFcn=statsFcn, stats=segs, what="segment", verbose=verbose);
    # Record statistics
    fitB$output <- segs;
    segs <- NULL; # Not needed anymore
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Summarizing change point (alpha, radius, manhattan, d1, d2) data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.element("changepoint", what)) {
    cps <- summarizeSamples(boot$changepoints, statsFcn=statsFcn, stats=cps, what="changepoint", verbose=verbose);
    # Record statistics
    fitB$changepoints <- cps;
    cps <- NULL; # Not needed anymore
  }

  # Not needed anymore
  fit <- boot <- NULL;

  verbose && exit(verbose);

  fitB;
}, private=TRUE) # bootstrapTCNandDHByRegion()





#   \item{by}{A @character specifying whether DH should be calculated from
#      normalized ('betaTN') or non-normalized ('betaT') tumor BAFs.}
#   \item{seed}{(optional) A random seed.}
#
#
# \value{
#   Returns a named @list containing two @arrays of bootstrap samples.
#   These arrays also contains the original observation as the first
#   element before the actual bootstrap samples.
# }
setMethodS3("bootstrapSegmentsAndChangepoints", "PairedPSCBS", function(fit, B=1000L, by=c("betaTN", "betaT"), seed=NULL, force=FALSE, cache=FALSE, verbose=FALSE, .debug=FALSE, ...) {
  # Settings for sanity checks
  tol <- getOption("PSCBS/sanityChecks/tolerance", 0.0005);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'B':
  B <- Arguments$getInteger(B, range=c(1,Inf));

  # Argument 'by':
  by <- match.arg(by);

  # Argument 'seed':
  if (!is.null(seed)) {
    seed <- Arguments$getInteger(seed);
  }

  # Argument '.cache':
  cache <- Arguments$getLogical(cache);

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Argument '.debug':
  .debug <- Arguments$getLogical(.debug);



  verbose && enter(verbose, "Bootstrapping (TCN,DH,C1,C2) segment mean levels");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  key <- list(method="bootstrapSegmentsAndChangepoints", class=class(fit)[1L],
              fit=fit, B=B, by=by, seed=seed);
  dirs <- c("PSCBS", "bootstrap");
  boot <- loadCache(key=key, dirs=dirs);
  if (!force && !is.null(boot)) {
    verbose && cat(verbose, "Found cached results. Skipping.");
    verbose && exit(verbose);
    return(boot);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Set the random seed
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(seed)) {
    randomSeed("set", seed=seed, kind="L'Ecuyer-CMRG")
    on.exit(randomSeed("reset"), add=TRUE)
    verbose && printf(verbose, "Random seed temporarily set (seed=%d)\n", seed)
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data and estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  data <- getLocusData(fit);
  tcnSegRows <- fit$tcnSegRows;
  dhSegRows <- fit$dhSegRows;
  segs <- getSegments(fit);
  params <- fit$params;

  # Sanity checks
  stopifnot(all(!is.na(data$chromosome) & !is.na(data$x)));

  # Sanity checks
  if (!params$joinSegments) {
    throw("Cannot bootstrap TCN and DH by segments unless PSCNs are segmented using joinSegments=TRUE.");
  }
  if (regexpr(",", params$flavor, fixed=TRUE) != -1L) {
    throw(sprintf("Cannot bootstrap TCN and DH by segments if PSCNs are segmented using flavor=\"%s\".", params$flavor));
  }
  # Sanity check (same as above, but just in case)
  stopifnot(all(segs$tcnStart == segs$dhStart, na.rm=TRUE));
  stopifnot(all(segs$tcnEnd == segs$dhEnd, na.rm=TRUE));



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Find estimators
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get mean estimators used
  estList <- getMeanEstimators(fit, c("tcn", "dh"));
  avgTCN <- estList$tcn;
  avgDH <- estList$dh;
  estList <- NULL; # Not needed anymore


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get (x,TCN,BAF) data
  chromosome <- data$chromosome;
  x <- data$x;
  CT <- data$CT;
  betaT <- data[[by]];
  muN <- data$muN;
  rho <- data$rho;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Classify each locus as (i) heterozygous SNP, (ii) homozygous SNP,
  # or (iii) non-polymorphic loci
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying heterozygous & homozygous SNPs and non-polymorphic loci");
  nbrOfLoci <- length(CT);
  verbose && cat(verbose, "Number of loci: ", nbrOfLoci);

  # Identify SNPs
  hasDH <- !is.null(rho)
  if (hasDH) {
    isHet <- !is.na(rho)
    isSnp <- isHet
  } else {
    isSnp <- (!is.na(muN) & !is.na(betaT))
    isHet <- (isSnp & (muN == 1/2))
  }

  snps <- which(isSnp);
  nonSNPs <- which(!isSnp);
  nbrOfSNPs <- sum(isSnp);
  nbrOfNonSNPs <- sum(!isSnp);
  verbose && cat(verbose, "Number of SNPs: ", nbrOfSNPs);
  verbose && cat(verbose, "Number of non-SNPs: ", nbrOfNonSNPs);

  # Sanity checks
  stopifnot(length(intersect(snps, nonSNPs)) == 0L);

  # Heterozygous SNPs
  hets <- which(isSnp &  isHet)
  homs <- which(isSnp & !isHet);
  nbrOfHets <- length(hets);
  nbrOfHoms <- length(homs);

  if (!hasDH) {
    verbose && printf(verbose, "Number of heterozygous SNPs: %d (%.2f%%)\n",
                                        nbrOfHets, 100*nbrOfHets/nbrOfSNPs);
    verbose && printf(verbose, "Number of homozygous SNPs: %d (%.2f%%)\n",
                                        nbrOfHoms, 100*nbrOfHoms/nbrOfSNPs);
  }

  # Sanity checks
  stopifnot(length(intersect(hets, homs)) == 0L);
  stopifnot(nbrOfHets + nbrOfHoms == nbrOfSNPs);

  # Sanity checks
  stopifnot(length(isSnp) == nbrOfLoci);
  stopifnot(length(isHet) == nbrOfLoci);

  # Not needed anymore
  muN <- isSnp <- NULL
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Precalculate DH signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(rho)) {
    # Calculate DHs for heterozygous SNPs
    rho <- 2*abs(betaT - 1/2);

    # DH is by definition only defined for heterozygous SNPs.
    # For simplicity, we set it to be NA for non-heterozygous loci.
    rho[!isHet] <- NA_real_;

    data$rho <- rho;
  }

  # Not needed anymore
  betaT <- isHet <- NULL;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Resample (TCN,DH) within each segments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfSegments <- nrow(segs);

  # Allocate JxBx4 matrix M of bootstrap means
  dim <- c(nbrOfSegments, B, 4L);
  dimnames <- list(NULL, NULL, c("tcn", "dh", "c1", "c2"));
  M <- array(NA_real_, dim=dim, dimnames=dimnames);
  verbose && str(verbose, M);

  # Identify all loci with non-missing signals
  idxsCT <- which(!is.na(CT));
  idxsRho <- which(!is.na(rho));

  # Vectorized pre-adjustments
  for (field in c("tcnNbrOfLoci", "dhNbrOfLoci")) {
    counts <- segs[[field]];
    counts[is.na(counts)] <- 0L;
    segs[[field]] <- counts;
  }

  hasTcnLoci <- (is.finite(tcnSegRows[,1L]) & is.finite(tcnSegRows[,2L]));
  hasDhLoci <- (is.finite(dhSegRows[,1L]) & is.finite(dhSegRows[,2L]));

  # Identify "splitter" segments which have no data
  chrs <- segs[["chromosome"]];
  tcnIds <- segs[["tcnId"]];
  dhIds <- segs[["dhId"]];
  tcnMeans <- segs[["tcnMean"]];
  dhMeans <- segs[["dhMean"]];
  isSplitter <- (is.na(chrs) & is.na(tcnIds) & is.na(dhIds));

  # Get all segment indices except for "splitters"
  jjs <- seq(length=nbrOfSegments);
  jjs <- jjs[!isSplitter];

  for (jj in jjs) {
    chr <- chrs[jj];
    tcnId <- tcnIds[jj];
    dhId <- dhIds[jj];

    verbose && enter(verbose, sprintf("Segment #%d (chr %d, tcnId=%d, dhId=%d) of %d", jj, chr, tcnId, dhId, nbrOfSegments));

    # Sanity check
    if (.debug) stopifnot(!is.na(chr) && !is.na(tcnId) && !is.na(dhId));

    # Get the segment data
    segJJ <- segs[jj,,drop=FALSE];
    verbose && print(verbose, segJJ);

    nbrOfTCNs <- segJJ[,"tcnNbrOfLoci"];
    nbrOfDHs <- segJJ[,"dhNbrOfLoci"];
    verbose && cat(verbose, "Number of TCNs: ", nbrOfTCNs);
    verbose && cat(verbose, "Number of DHs: ", nbrOfDHs);
    if (.debug) {
      stopifnot(!is.na(nbrOfTCNs));
      stopifnot(!is.na(nbrOfDHs));
    }

    tcnSegRowJJ <- unlist(tcnSegRows[jj,], use.names=FALSE);
    dhSegRowJJ <- unlist(dhSegRows[jj,], use.names=FALSE);

    # Indices of all loci
    if (hasTcnLoci[jj]) {
      idxsAll <- tcnSegRowJJ[1L]:tcnSegRowJJ[2L];
    } else {
      idxsAll <- 0L;
    }

    verbose && str(verbose, idxsAll);
    verbose && print(verbose, hpaste(idxsAll), level=-120);
    verbose && str(verbose, idxsCT);
    verbose && print(verbose, hpaste(idxsCT), level=-120);

    # Keep only loci with finite TCNs
    idxsAll <- intersect(idxsAll, idxsCT);
    verbose && str(verbose, idxsAll);
    verbose && print(verbose, hpaste(idxsAll), level=-120);

    # Sanity check
    if (length(idxsAll) != nbrOfTCNs) {
      verbose && str(verbose, setdiff(idxsCT, idxsAll));
      throw("INTERNAL ERROR: length(idxsAll) != nbrOfTCNs: ", length(idxsAll), " != ", nbrOfTCNs);
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Identify loci used to calculate DH means
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Identify loci used to bootstrap DH means");

    if (hasDhLoci[jj]) {
      idxsDH <- dhSegRowJJ[1L]:dhSegRowJJ[2L];
      idxsDH <- intersect(idxsDH, hets);
      # Drop missing values
      idxsDH <- intersect(idxsDH, idxsRho);
    } else {
      idxsDH <- 0L;
    }

    verbose && cat(verbose, "Heterozygous SNPs to resample for DH:");
    verbose && str(verbose, idxsDH);

    # Sanity check
    if (.debug) stopifnot(length(idxsDH) == nbrOfDHs);

    verbose && exit(verbose);



    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Identify loci used to calculate TCN means
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Identify loci used to bootstrap TCN means");

    # Identify SNPs and non-SNPs
    idxsSNP <- intersect(snps, idxsAll);
    idxsNonSNP <- setdiff(idxsAll, idxsSNP);
    verbose && cat(verbose, "SNPs:");
    verbose && str(verbose, idxsSNP);
    verbose && cat(verbose, "Non-polymorphic loci:");
    verbose && str(verbose, idxsNonSNP);
    # Sanity check
    if (.debug) stopifnot(length(idxsSNP) + length(idxsNonSNP) == length(idxsAll));

    # Identify heterozygous and homozygous SNPs
    idxsHet <- intersect(idxsSNP, hets);
    if (nbrOfHoms > 0) {
      idxsHom <- intersect(idxsSNP, homs)
    } else {
      ## Happens when only DH is available
      idxsHom <- integer(0L)
    }

    # Drop missing values
    idxsNonSNP <- intersect(idxsNonSNP, idxsCT);
    idxsHet <- intersect(idxsHet, idxsCT);
    idxsHom <- intersect(idxsHom, idxsCT);
    idxsHetNonDH <- setdiff(idxsHet, idxsDH);

    verbose && cat(verbose, "Heterozygous SNPs to resample for TCN:");
    verbose && str(verbose, idxsHet);
    verbose && cat(verbose, "Homozygous SNPs to resample for TCN:");
    verbose && str(verbose, idxsHom);
        verbose && cat(verbose, "Non-polymorphic loci to resample for TCN:");
    verbose && str(verbose, idxsNonSNP);
    verbose && cat(verbose, "Heterozygous SNPs with non-DH to resample for TCN:");
    verbose && str(verbose, idxsHetNonDH);
    # Note that length(idxsHetNonDH) may differ from zero in case CT is non-missing
    # but rho is missing, e.g. CT = sum(c(thetaA,thetaB), na.rm=TRUE) and
    # thetaB is missing. /HB 2010-12-01

    idxsTCN <- sort(unique(c(idxsHet, idxsHom, idxsNonSNP)));
    verbose && cat(verbose, "Loci to resample for TCN:");
    verbose && str(verbose, idxsTCN);

    # Sanity check
    if (.debug) {
      stopifnot(length(idxsHet) + length(idxsHom) + length(idxsNonSNP) == nbrOfTCNs);
      stopifnot(length(intersect(idxsDH, idxsHetNonDH)) == 0L);
      stopifnot(length(idxsTCN) == nbrOfTCNs);
    }

    verbose && exit(verbose);


    # These numbers should be preserved when the resampling
    verbose && printf(verbose, "Number of (#hets, #homs, #nonSNPs): (%d,%d,%d)\n",
                      length(idxsHet), length(idxsHom), length(idxsNonSNP));


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Sanity checks
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (nbrOfTCNs > 0L) {
      # Sanity check
      ys <- CT[idxsTCN];
      mu <- avgTCN(ys, na.rm=TRUE);
      dMu <- (mu - tcnMeans[jj]);
      if (abs(dMu) > tol) {
        str(list(nbrOfTCNs=nbrOfTCNs, tcnNbrOfLoci=segJJ$tcnNbrOfLoci, mu=mu, tcnMean=tcnMeans[jj], dMu=dMu, "abs(dMu)"=abs(dMu), "range(x[units])"=range(x[idxsTCN])));
        throw(sprintf("INTERNAL ERROR: Incorrectly recalculated TCN mean for Segment #%d (chr %d, tcnId=%d, dhId=%d): %g != %g", jj, chr, tcnId, dhId, mu, tcnMeans[jj]));
      }
    }

    shouldHaveDHs <- (nbrOfDHs > 0L && !is.na(dhMeans[jj]));
    if (shouldHaveDHs) {
      # Sanity check
      ys <- rho[idxsDH];
      mu <- avgDH(ys, na.rm=TRUE);
      dMu <- (mu - dhMeans[jj]);
      if (abs(dMu) > tol) {
        str(list(nbrOfDHs=nbrOfDHs, dhNbrOfLoci=segJJ$dhNbrOfLoci, mu=mu, dhMean=dhMeans[jj], dMu=dMu, "abs(dMu)"=abs(dMu), "range(x[units])"=range(x[idxsDH])));
        throw(sprintf("INTERNAL ERROR: Incorrectly recalculated DH mean for Segment #%d (chr %d, tcnId=%d, dhId=%d): %g != %g", jj, chr, tcnId, dhId, mu, dhMeans[jj]));
      }
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Bootstrap while preserving (#hets, #homs, #nonSNPs)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Bootstrapping while preserving (#hets, #homs, #nonSNPs)");
    verbose && cat(verbose, "Number of bootstrap samples: ", B);

    if (!shouldHaveDHs) {
      idxsHetNonDH <- idxsDH;
    }

    nHoms <- length(idxsHom);
    nNonSNPs <- length(idxsNonSNP);
    nHetNonDHs <- length(idxsHetNonDH);

    # Defaults
    idxsDHBB <- NULL;

    # Bootstrap B times
    for (bb in seq(length=B)) {
      # (1) Bootstrap DHs
      if (shouldHaveDHs) {
        # (a) Resample heterozygous SNPs (=> resampled DH units)
        idxsDHBB <- resample(idxsDH, size=nbrOfDHs, replace=TRUE);

        # Extract signals
        rhoBB <- rho[idxsDHBB];

        # Calculate bootstrap mean
        M[jj,bb,"dh"] <- avgDH(rhoBB, na.rm=TRUE);
      } # if (shouldHaveDHs)

      # (2) Bootstrap TCNs
      if (nbrOfTCNs > 0L) {
        # (a) Resample non-DH hets SNPs
        idxsHetNonDHBB <- resample(idxsHetNonDH, size=nHetNonDHs, replace=TRUE);
        idxsHetBB <- c(idxsDHBB, idxsHetNonDHBB);

        # (a) Resample homozygous SNPs
        if (nbrOfHoms > 0) {
          idxsHomBB <- resample(idxsHom, size=nHoms, replace=TRUE)
        } else {
          idxsHomBB <- integer(0L)
        }

        # (b) Resample non-SNPs
        idxsNonSNPBB <- resample(idxsNonSNP, size=nNonSNPs, replace=TRUE);

        # (c) Resampled TCN units
        idxsTCNBB <- c(idxsHetBB, idxsHomBB, idxsNonSNPBB);

        # Sanity check
        if (.debug) {
          stopifnot(length(intersect(idxsDHBB, idxsHetNonDHBB)) == 0L);
          stopifnot(length(intersect(idxsHetBB, idxsHomBB)) == 0L);
          stopifnot(length(intersect(idxsHetBB, idxsNonSNPBB)) == 0L);
          stopifnot(length(intersect(idxsHomBB, idxsNonSNPBB)) == 0L);
        }

        # Extract signals
        CTBB <- CT[idxsTCNBB];

        # Calculate bootstrap mean
        M[jj,bb,"tcn"] <- avgTCN(CTBB, na.rm=TRUE);
      } # if (nbrOfTCNs > 0L)
    } # (for bb ...)
    verbose && exit(verbose);

    verbose && exit(verbose);
  } # for (jj ...)

  verbose && cat(verbose, "Bootstrapped segment mean levels");
  verbose && str(verbose, M);

  # Sanity check
  stopifnot(all(!is.nan(M)));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Add (C1,C2) bootstrap mean levels
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calculating (C1,C2) mean levels from (TCN,DH) mean levels");
  C1 <- (1-M[,,"dh", drop=FALSE]) * M[,,"tcn", drop=FALSE] / 2;
  C2 <- M[,,"tcn", drop=FALSE] - C1;
  M[,,"c1"] <- C1;
  M[,,"c2"] <- C2;
  verbose && str(verbose, M);
  # Sanity check
  stopifnot(dim(M)[1L] == nbrOfSegments);
  stopifnot(all(!is.nan(M)));
  # Not needed anymore
  C1 <- C2 <- NULL;
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Bootstrap polar (alpha,radius,manhattan) for change points
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calculating polar (alpha,radius,manhattan) for change points");
  C <- M[,,c("c1","c2"), drop=FALSE];

  # Calculate difference
  # (will be empty if nbrOfSegments == 1, but that's ok/intended)
  D <- C[-nbrOfSegments,,, drop=FALSE] - C[-1L,,, drop=FALSE];
  # Sanity check
  stopifnot(dim(D)[1L] == nbrOfSegments-1L);
  stopifnot(all(!is.nan(D)));
  C <- NULL; # Not needed anymore

  # Allocate array
  dimnames <- dimnames(D);
  dimnames[[3L]] <- c("alpha", "radius", "manhattan", "d1", "d2");
  dim <- dim(D);
  dim[3L] <- length(dimnames[[3L]]);
  P <- array(NA_real_, dim=dim, dimnames=dimnames);
  stopifnot(dim(P)[1L] == nbrOfSegments-1L);

  if (nbrOfSegments >= 2L) {
    verbose && str(verbose, D);
    P[,,"alpha"] <- atan2(D[,,2], D[,,1]); # Changepoint angles in (0,2*pi)
    P[,,"radius"] <- sqrt(D[,,2]^2 + D[,,1]^2);
    P[,,"manhattan"] <- abs(D[,,2]) + abs(D[,,1]);
    P[,,"d1"] <- D[,,1];
    P[,,"d2"] <- D[,,2];
  }
  alpha <- D <- NULL; # Not needed anymore
  verbose && cat(verbose, "Bootstrapped change points");
  verbose && str(verbose, P);

  # Sanity check
  stopifnot(dim(P)[1L] == nbrOfSegments-1L);
  stopifnot(all(!is.nan(P)));
  verbose && exit(verbose);

  boot <- list(segments=M, changepoints=P);

  # Cache?
  if (cache) {
    saveCache(boot, key=key, dirs=dirs);
    verbose && cat(verbose, "Saved results to cache.");
  }

  verbose && exit(verbose);

  boot;
}, private=TRUE) # bootstrapSegmentsAndChangepoints()



setMethodS3("findBootstrapSummaries", "PairedPSCBS", function(fit, what=c("segment", "changepoint"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'what':
  what <- match.arg(what);

  if (what == "segment") {
    data <- getSegments(fit);
  } else if (what == "changepoint") {
    data <- getChangePoints(fit);
  }

  grep("^[^_]+_[^_]+$", colnames(data), value=TRUE);
}, protected=TRUE) # findBootstrapSummaries()


setMethodS3("hasBootstrapSummaries", "PairedPSCBS", function(fit, ...) {
  fields <- findBootstrapSummaries(fit, ...);
  (length(fields) > 0L);
})

setMethodS3("clearBootstrapSummaries", "PairedPSCBS", function(fit, what=c("segment", "changepoint"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'what':
  what <- unique(match.arg(what, several.ok=TRUE));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Clearing bootstrap summaries");

  whats <- what;
  for (what in whats) {
    verbose && enter(verbose, sprintf("Clearing %ss", what));

    fields <- findBootstrapSummaries(fit, what=what, ...);
    if (what == "segment") {
      data <- getSegments(fit);
      data <- data[,!is.element(colnames(data), fields)];
      fit$output <- data;
    } else if (what == "changepoint") {
      data <- getChangePoints(fit);
      data <- data[,!is.element(colnames(data), fields)];
      fit$changepoints <- data;
    }
    # Sanity check
    fields <- findBootstrapSummaries(fit, what=what, ...);
    stopifnot(length(fields) == 0L);

    data <- fields <- NULL; # Not needed anymoew

    verbose && exit(verbose);
  } # for (what ...)

  verbose && exit(verbose);

  fit;
}, protected=TRUE) # clearBootstrapSummaries()


##############################################################################
# HISTORY
# 2013-10-22
# o SPEEDUP: Added argument 'cache' to bootstrapSegmentsAndChangepoints(),
#   which caches the results to file if cache=TRUE.
# 2013-10-21
# o Added find-, has- and clearBootstrapSummaries().
# o Added argument 'what' to bootstrapTCNandDHByRegion().
# o BUG FIX: The new bootstrapSegmentsAndChangepoints() would give
#   "Error in D[, , 2] : incorrect number of dimensions" if bootstrapping
#   was done on a single change point (==two segments).
# 2013-10-20
# o Now utilizing new getChangePoints().
# o Now calculating change-point angles in (0,2*pi) using atan2().
# o Added bootstrapSegmentsAndChangepoints(), which was extract from
#   internal code of bootstrapTCNandDHByRegion().  The latter now
#   utilizes the former.
# o BUG FIX: bootstrapTCNandDHByRegion() did not identify segment
#   "splitters" as intended.  This has had no impact on the results.
# 2013-04-23
# o SPEEDUP: Made bootstrapTCNandDHByRegion() much faster by adding
#   use.names=FALSE to two internal unlist() statements.
# 2013-04-22
# o Added argument 'probs' to bootstrapTCNandDHByRegion().
# 2013-04-09
# o Added Rdoc comments.
# 2013-02-09
# o BUG FIX: bootstrapTCNandDHByRegion() for PairedPSCBS did not bootstrap
#   for all available loci when calculating total CNs statistics, iff
#   the segment had been called run-of-homozygosity (ROH).
#   Thanks to Oscar Rueda at the Cancer Research UK Cambridge Institute
#   for reporting on this.
# 2013-02-07
# o Improved some verbose outputs of bootstrapTCNandDHByRegion().
# 2013-01-15
# o Now bootstrapTCNandDHByRegion() uses the params$avgTCN and params$avgDH
#   estimators, iff given.
# 2012-11-05
# o GENERALIZATION: Now bootstrapTCNandDHByRegion() works for more "flavors",
#   including the default ('tcn') used by segmentByNonPairedPSCBS().
# 2012-09-20
# o SPEEDUP: By default bootstrapTCNandDHByRegion() for PairedPSCBS no
#   longer do sanity checks within the bootstrap loop.  This significantly
#   speed up the method.  To run checks, use argument .debug=TRUE.
# 2012-02-26
# o BUG FIX: bootstrapTCNandDHByRegion() for PairedPSCBS would resample
#   from a subset of the intended TCNs, iff the DH mean was non-finite
#   while there were still heterozygous SNPs.  This introduced a bias in
#   the estimates, which was neglectable for large segments, but for very
#   small segments (a few loci) it could be relatively large.
# 2012-02-24
# o BUG FIX: bootstrapTCNandDHByRegion(..., force=TRUE) for PairedPSCBS
#   would give an error iff previous bootstrap estimates already existed.
# o Added argument 'force' to bootstrapTCNandDHByRegion().
# 2011-11-26
# o Now bootstrapTCNandDHByRegion() for PairedPSCBS preserves NAs for DH
#   and (C1,C2) quantiles, if the DH mean level is NA, which can happen
#   when a segment is called ROH.
# o An internal sanity check of bootstrapTCNandDHByRegion() for PairedPSCBS
#   would give an error if DH mean levels had been set to NA for segments
#   called ROH.
# 2011-11-24
# o BUG FIX: bootstrapTCNandDHByRegion() for PairedPSCBS would give
#   an error, if a segment did not have any TCN signals, which can
#   occur when known segments are specified for Paired PSCBS.
# 2011-08-08
# o Moved the sanity checks that tests the TCN and DH "segRows" from the
#   bootstrapTCNandDHByRegion() to segmentByPairedPSCBS().  This is the
#   first step to fix a failure in the sanity checks that could happend
#   iff one first run dropSegmentationOutliers().
# 2011-06-14
# o Updated code to recognize new column names.
# 2011-05-29
# o Renamed options to reflect new package name.
# 2010-12-03
# o BUG FIX: In rare cases the bootstrap sanity checks can indeed produce
#   an invalid 'range', more precisely where (range[,2] >= range[,1]) is
#   not true.  This can happen if there is no variation in the bootstrap
#   estimates.  Beause of this we allow for some tolerance.
# 2010-12-02
# o Now bootstrapTCNandDHByRegion() uses option
#   "psCBS/sanityChecks/tolerance".
# 2010-12-01
# o BUG FIX: bootstrapTCNandDHByRegion() did not always exclude the correct
#   loci.
# 2010-11-27
# o BUG FIX: bootstrapTCNandDHByRegion() would incorrectly include
#   non-polymorphic loci in the set of homozygous SNPs during resampling.
# 2010-11-26
# o BUG FIX: The statistical sanity checks of the bootstrap estimates would
#   give an error when only single-sided bootstrap confidence interval was
#   calculated.
# 2010-11-23
# o ROBUSTNESS: Added more sanity checks to bootstrapTCNandDHByRegion().
# o WORKAROUND: The precision of the mean levels of DNAcopy::segment()
#   is not great enough to always compare it to that of R's estimates.
# o BUG FIX: bootstrapTCNandDHByRegion() would give an error if there was
#   only one segment.
# 2010-11-22
# o BUG FIX: The DH segmentation and bootstrap incorrectly included
#   missing values, when subseting.
# o BUG FIX: Some sanity checks were incorrect.
# o BUG FIX: bootstrapTCNandDHByRegion() for PairedPSCBS would not correctly
#   detect if bootstrap results are already available.
# 2010-11-21
# o Added argument 'seed'.
# o Added bootstrapTCNandDHByRegion() for PairedPSCBS.
# o Created from PairedPSCBS.BOOT.R.
##############################################################################
