###########################################################################/**
# @set class=PairedPSCBS
# @RdocMethod getBootstrapLocusSets
# @alias getBootstrapLocusSets
#
# @title "Generates original and bootstrapped segment-specific index sets"
#
# \description{
#  @get "title", which can be used to calculate various bootstrap summaries,
#  e.g. segment mean levels.
# }
#
# @synopsis
#
# \arguments{
#   \item{B}{A non-negative @integer specifying the number of bootstrap samples.}
#   \item{by}{Should \code{betaTN} or \code{betaT} be used?}
#   \item{seed}{An (optional) @integer specifying the random seed to be
#     set before sampling indices.  The random seed is set to its original
#     state when exiting.  If @NULL, it is not set.}
#   \item{verbose}{See @see "R.utils::Verbose".}
#   \item{.validate}{If @TRUE, additional sanity checks are performed
#     to validate the correctness.  This is only needed for troubleshooting
#     if it is suspected there is a bug.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @list.
# }
#
# @author "HB"
#
# \seealso{
#   This is used internally by various bootstrap methods.
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("getBootstrapLocusSets", "PairedPSCBS", function(fit, B=1000L, by=c("betaTN", "betaT"), seed=NULL, verbose=FALSE, .validate=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'B':
  B <- Arguments$getInteger(B, range=c(0,Inf));

  # Argument 'by':
  by <- match.arg(by);

  # Argument 'seed':
  if (!is.null(seed)) {
    seed <- Arguments$getInteger(seed);
  }

  # Argument '.validate':
  .validate <- Arguments$getLogical(.validate);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  verbose && enter(verbose, "Bootstrapping (TCN,DH,C1,C2) segment mean levels");

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
  # Get signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get (x,TCN,BAF) data
  chromosome <- data$chromosome;
  x <- data$x;
  CT <- data$CT;
  betaT <- data[[by]];
  muN <- data$muN;
  rho <- data$rho
  hasDH <- !is.null(rho)

  # Not needed anymore
  data <- NULL;

  # Sanity checks
  stopifnot(all(!is.na(chromosome) & !is.na(x)));



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Classify each locus as (i) heterozygous SNP, (ii) homozygous SNP,
  # or (iii) non-polymorphic loci
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying heterozygous & homozygous SNPs and non-polymorphic loci");
  nbrOfLoci <- length(x);
  verbose && cat(verbose, "Number of loci: ", nbrOfLoci);

  # SNPs are identifies as those loci that have non-missing 'muN' (& betaTN')
  if (hasDH) {
    isHet <- !is.na(rho)
    isSnp <- isHet
  } else {
    isSnp <- (!is.na(muN) & !is.na(betaT))
    isHet <- isSnp & (muN == 1/2)
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
  hets <- which(isSnp &  isHet);
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
  muN <- isSnp <- NULL;
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Precalculate DH signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!hasDH) {
    # Calculate DHs for heterozygous SNPs
    rho <- 2*abs(betaT - 1/2);

    # DH is by definition only defined for heterozygous SNPs.
    # For simplicity, we set it to be NA for non-heterozygous loci.
    rho[!isHet] <- NA;
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

  # Identify all loci with non-missing signals
  idxsCT <- which(!is.na(CT));
  idxsRho <- which(!is.na(rho));

  # Not needed anymore
  CT <- rho <- NULL;


  # Vectorized pre-adjustments
  for (key in c("tcnNbrOfLoci", "dhNbrOfLoci")) {
    counts <- segs[[key]];
    counts[is.na(counts)] <- 0L;
    segs[[key]] <- counts;
    counts <- NULL; # Not needed anymore
  }

  hasTcnLoci <- (is.finite(tcnSegRows[,1L]) & is.finite(tcnSegRows[,2L]));
  hasDhLoci <- (is.finite(dhSegRows[,1L]) & is.finite(dhSegRows[,2L]));

  # Identify "splitter" segments which have no data
  chrs <- segs[["chromosome"]];
  tcnIds <- segs[["tcnId"]];
  dhIds <- segs[["dhId"]];
  dhMeans <- segs[["dhMean"]];
  isSplitter <- (is.na(chrs) & is.na(tcnIds) & is.na(dhIds));

  # Get all segment indices except for "splitters"
  jjs <- seq(length=nbrOfSegments);
  jjs <- jjs[!isSplitter];


  # Allocate list to hold the results
  res <- list(
    nbrOfSegments = nbrOfSegments,
    segments = jjs,
    nbrOfLoci = nbrOfLoci,
    loci = list(
      all = seq_len(nbrOfLoci),
      tcn = idxsCT,
      dh  = idxsRho
    ),
    by = by,
    seed = seed
  );
  locusSet <- vector("list", length=nbrOfSegments);


  # For each segment jj = 1, 2, ..., S
  for (jj in jjs) {
    chr <- chrs[jj];
    tcnId <- tcnIds[jj];
    dhId <- dhIds[jj];

    verbose && enter(verbose, sprintf("Segment #%d (chr %d, tcnId=%d, dhId=%d) of %d", jj, chr, tcnId, dhId, nbrOfSegments));

    # Sanity check
    if (.validate) {
      stopifnot(!is.na(chr) && !is.na(tcnId) && !is.na(dhId));
    }

    # Get the segment data
    segJJ <- segs[jj,,drop=FALSE];
    nbrOfTCNs <- segJJ[,"tcnNbrOfLoci"];
    nbrOfDHs <- segJJ[,"dhNbrOfLoci"];

    if (verbose) {
      print(verbose, segJJ);
      cat(verbose, "Number of TCNs: ", nbrOfTCNs);
      cat(verbose, "Number of DHs: ", nbrOfDHs);
    }
    if (.validate) {
      stopifnot(!is.na(nbrOfTCNs));
      stopifnot(!is.na(nbrOfDHs));
    }

    tcnSegRowJJ <- unlist(tcnSegRows[jj,], use.names=FALSE);
    dhSegRowJJ <- unlist(dhSegRows[jj,], use.names=FALSE);

    # Indices of all loci
    if (hasTcnLoci[jj]) {
      idxsAll <- tcnSegRowJJ[1L]:tcnSegRowJJ[2L];
    } else {
      idxsAll <- integer(0L);
    }

    if (verbose) {
      str(verbose, idxsAll);
      print(verbose, hpaste(idxsAll), level=-120);
      str(verbose, idxsCT);
      print(verbose, hpaste(idxsCT), level=-120);
    }

    # Keep only loci with finite TCNs
    idxsAll <- intersect(idxsAll, idxsCT);
    if (verbose) {
      str(verbose, idxsAll);
      print(verbose, hpaste(idxsAll), level=-120);
    }

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
      idxsDH <- integer(0L);
    }

    verbose && cat(verbose, "Heterozygous SNPs to resample for DH:");
    verbose && str(verbose, idxsDH);

    # Sanity check
    if (.validate) stopifnot(length(idxsDH) == nbrOfDHs);

    verbose && exit(verbose);



    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Identify loci used to calculate TCN means
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Identify loci used to bootstrap TCN means");

    # Identify SNPs and non-SNPs
    idxsSNP <- intersect(snps, idxsAll);
    idxsNonSNP <- setdiff(idxsAll, idxsSNP);
    if (verbose) {
      cat(verbose, "SNPs:");
      str(verbose, idxsSNP);
      cat(verbose, "Non-polymorphic loci:");
      str(verbose, idxsNonSNP);
    }
    # Sanity check
    if (.validate) stopifnot(length(idxsSNP) + length(idxsNonSNP) == length(idxsAll));

    # Identify heterozygous and homozygous SNPs
    idxsHet <- intersect(idxsSNP, hets);
    if (nbrOfHoms > 0) {
      idxsHom <- intersect(idxsSNP, homs);
    } else {
      idxsHom <- integer(0L)
    }

    # Drop missing values
    idxsNonSNP <- intersect(idxsNonSNP, idxsCT);
    idxsHet <- intersect(idxsHet, idxsCT);
    if (nbrOfHoms > 0) {
      idxsHom <- intersect(idxsHom, idxsCT)
    }
    idxsHetNonDH <- setdiff(idxsHet, idxsDH);

    if (verbose) {
      cat(verbose, "Heterozygous SNPs to resample for TCN:");
      str(verbose, idxsHet);
      cat(verbose, "Homozygous SNPs to resample for TCN:");
      str(verbose, idxsHom);
      cat(verbose, "Non-polymorphic loci to resample for TCN:");
      str(verbose, idxsNonSNP);
      cat(verbose, "Heterozygous SNPs with non-DH to resample for TCN:");
      str(verbose, idxsHetNonDH);
    }
    # Note that length(idxsHetNonDH) may differ from zero in case CT is non-missing
    # but rho is missing, e.g. CT = sum(c(thetaA,thetaB), na.rm=TRUE) and
    # thetaB is missing. /HB 2010-12-01

    idxsTCN <- sort(unique(c(idxsHet, idxsHom, idxsNonSNP)));
    if (verbose) {
      cat(verbose, "Loci to resample for TCN:");
      str(verbose, idxsTCN);
    }

    # Sanity check
    if (.validate) {
      stopifnot(length(idxsHet) + length(idxsHom) + length(idxsNonSNP) == nbrOfTCNs);
      stopifnot(length(intersect(idxsDH, idxsHetNonDH)) == 0L);
      stopifnot(length(idxsTCN) == nbrOfTCNs);
    }

    verbose && exit(verbose);


    # These numbers should be preserved when the resampling
    verbose && printf(verbose, "Number of (#hets, #homs, #nonSNPs): (%d,%d,%d)\n",
                      length(idxsHet), length(idxsHom), length(idxsNonSNP));

    # Workaround: ... Why? /HB 2013-10-22
    shouldHaveDHs <- (nbrOfDHs > 0L && !is.na(dhMeans[jj]));
    if (!shouldHaveDHs) {
      idxsHetNonDH <- idxsDH;
      stopifnot(all(idxsHetNonDH > 0L));
    }
    shouldHaveDHs <- NULL; # Not needed anymore

    nHoms <- length(idxsHom);
    nNonSNPs <- length(idxsNonSNP);
    nHetNonDHs <- length(idxsHetNonDH);

    locusSetJJ <- list(
      segment = segJJ,
      loci = list(
        all = idxsAll,
        snp = idxsSNP,
        tcn = idxsTCN,
        dh = idxsDH,
        nonSnp = idxsNonSNP,
        het = idxsHet,
        hom = idxsHom,
        hetNonDh = idxsHetNonDH
      )
    );

    # Sanity checks
    if (.validate) {
      loci <- locusSetJJ$loci;
      for (key in names(loci)) {
        idxs <- loci[[key]];
        # Assert positive indices
        stopifnot(all(idxs > 0L));
        # Assert a unique set of indices
        stopifnot(!any(duplicated(idxs)));
      }
      # Assert non-overlapping sets
      with(loci, {
        stopifnot(length(intersect(dh,  hetNonDh)) == 0L);
        stopifnot(length(intersect(het, hom)) == 0L);
        stopifnot(length(intersect(het, nonSnp)) == 0L);
        stopifnot(length(intersect(hom, nonSnp)) == 0L);
      });
      loci <- NULL; # Not needed anymore
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Identify bootstrap locus sets preserving (#hets, #homs, #nonSNPs)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Identifying bootstrap locus sets that preservs (#hets, #homs, #nonSNPs)");
    verbose && cat(verbose, "Number of bootstrap samples: ", B);

    nbrOfTCNs <- nbrOfDHs+nHoms+nNonSNPs;
    nbrOfHets <- nbrOfDHs+nHetNonDHs;

    # Allocate index matrices
    tcn      <- matrix(NA_integer_, nrow=nbrOfTCNs, ncol=B);
    dh       <- matrix(NA_integer_, nrow=nbrOfDHs, ncol=B);
    nonSnp   <- matrix(NA_integer_, nrow=nNonSNPs, ncol=B);
    het      <- matrix(NA_integer_, nrow=nbrOfHets, ncol=B);
    hom      <- matrix(NA_integer_, nrow=nHoms, ncol=B);
    hetNonDh <- matrix(NA_integer_, nrow=nHetNonDHs, ncol=B);

##    resample <- function(x, size, ...) {
##      stopifnot(size == length(x));
##      x;
##    } # resample()

    # (1) Bootstrap DH loci
    if (nbrOfDHs > 0L) {
      # (a) Resample heterozygous SNPs (=> resampled DH units)
      for (bb in seq_len(B)) {
        dh[,bb] <- resample(idxsDH, size=nbrOfDHs, replace=TRUE);
      }
    }

    # (2) Bootstrap other loci
    # (a) Resample non-DH hets SNPs
    if (nHetNonDHs > 0L) {
      for (bb in seq_len(B)) {
        idxsHetNonDHBB <- resample(idxsHetNonDH, size=nHetNonDHs, replace=TRUE);
        hetNonDh[,bb] <- idxsHetNonDHBB;
      }
    }

    # (b) Resample homozygous SNPs
    if (nHoms > 0L) {
      for (bb in seq_len(B)) {
        idxsHomBB <- resample(idxsHom, size=nHoms, replace=TRUE);
        hom[,bb] <- idxsHomBB;
      }
    }

    # (c) Resample non-SNPs
    if (nNonSNPs > 0L) {
      for (bb in seq_len(B)) {
        idxsNonSNPBB <- resample(idxsNonSNP, size=nNonSNPs, replace=TRUE);
        nonSnp[,bb] <- idxsNonSNPBB;
      }
    }

    # (d) Resampled hets
    if (nbrOfHets > 0L) {
      for (bb in seq_len(B)) {
        idxsDHBB <- dh[,bb];
        idxsHetNonDHBB <- hetNonDh[,bb];
        idxsHetBB <- c(idxsDHBB, idxsHetNonDHBB);
#        idxsHetBB <- sort(idxsHetBB);
        het[,bb] <- idxsHetBB;
      }
    }

    # (e) Update TCN loci
    if (nbrOfTCNs > 0L) {
      for (bb in seq_len(B)) {
        idxsHetBB <- het[,bb];
        idxsHomBB <- hom[,bb];
        idxsNonSNPBB <- nonSnp[,bb];
        idxsTCNBB <- c(idxsHetBB, idxsHomBB, idxsNonSNPBB);
#        idxsTCNBB <- sort(idxsTCNBB);
        tcn[,bb] <- idxsTCNBB;
      }
    }

    # Record
    locusSetJJ$bootstrap <- list(
      B    = B,
      loci = list(
        tcn      = tcn,
        dh       = dh,
        nonSnp   = nonSnp,
        het      = het,
        hom      = hom,
        hetNonDh = hetNonDh
      )
    );

    # Sanity check
    if (.validate) {
      loci <- locusSetJJ$loci;
      lociB <- locusSetJJ$bootstrap$loci;
      for (key in names(lociB)) {
        idxs <- loci[[key]];
        idxsB <- lociB[[key]];
        idxsB <- unique(sort(idxsB));
        stopifnot(all(is.element(idxsB, idxs)));
      }
      loci <- lociB <- NULL; # Not needed anymore
    }

    verbose && exit(verbose);

    # Record
    locusSet[[jj]] <- locusSetJJ;

    # Not needed anymore
    locusSetJJ <- NULL;
    verbose && exit(verbose);
  } # for (jj ...)

  # Sanity checks
  stopifnot(is.list(locusSet));
  stopifnot(length(locusSet) == nbrOfSegments);

  verbose && exit(verbose);

  res$locusSet <- locusSet;

  res;
}, protected=TRUE) # getBootstrapLocusSets()


##############################################################################
# HISTORY
# 2013-10-22
# o Added Rdoc comments.
# o Added getBootstrapLocusSets() for PairedPSCBS.
##############################################################################
