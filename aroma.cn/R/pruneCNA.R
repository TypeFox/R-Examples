# \section{Required method implementations}{
#  In order for this method to work, the following methods need to be
#  implemented for the class of argument \code{fit}:
#  \itemize{
#   \item \code{findAtomicAberrations()}
#   \item \code{mergeTwoSegments()}
#  }
# }
setMethodS3("pruneCNA", "AbstractCBS", function(fit, ..., maxGeneration=Inf, onAtomicIsland=NULL, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'maxGeneration':
  maxGeneration <- Arguments$getDouble(maxGeneration, range=c(0,Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  fitT <- fit;

  fitList <- list();
  gg <- 1;
  while (gg <= maxGeneration) {
    verbose && enter(verbose, sprintf("Generation #%d", gg));

    fitList[[gg]] <- fitT;

    nbrOfSegments <- nbrOfSegments(fitT);
    maxH <- nbrOfSegments-2L;

    # Done
    if (maxH < 0) {
      break;
    }

    hasChanged <- FALSE;
    for (hh in 0:maxH) {
#    for (hh in 0) {
      verbose && enter(verbose, sprintf("Block size H=%d of %d", hh, maxH));

      res <- findAtomicAberrations(fitT, H=hh, ..., verbose=verbose);
      verbose && str(verbose, res);

      # (i) Atomic islands?
      if (hh == 0) {
        atomicIslands <- res$atomicRegions;
      } else {
        atomicIslands <- res$atomicIslands;
      }

      if (length(atomicIslands) > 0) {
        verbose && printf(verbose, "Atomic islands found (H=%d):\n", hh);
        verbose && print(verbose, atomicIslands);

        # Overlapping atomic islands?
        if (length(atomicIslands) > 1) {
          regionsHH <- matrix(c(atomicIslands, atomicIslands+hh), ncol=2L, byrow=FALSE);
          colnames(regionsHH) <- c("from", "to");
          rownames(regionsHH) <- sprintf("Atomic island #%d", seq_along(atomicIslands));
          verbose && print(verbose, regionsHH);

          froms <- regionsHH[-1,"from"];
          tos <- regionsHH[-length(atomicIslands),"to"];
          isOverlapping <- (froms <= tos);
          verbose && printf(verbose, "Overlapping: %s\n", any(isOverlapping));

          if (any(isOverlapping)) {
            verbose && cat(verbose, "Overlapping atomic islands. Dropping only the first.");
            atomicIslands <- atomicIslands[1];
          }
        }

        # Drop atomic islands and merge flanking segments
        dropList <- list();
        atomicIslands <- sort(atomicIslands, decreasing=TRUE);
        for (kk in seq_along(atomicIslands)) {
          atomicIsland <- atomicIslands[kk];
          if (hh == 0) {
            atomicIslandTag <- sprintf("change point #%d", atomicIsland);
          } else if (hh == 1) {
            atomicIslandTag <- sprintf("segment #%d", atomicIsland);
          } else {
            atomicIslandTag <- sprintf("segments #%d-#%d", atomicIsland, atomicIsland+(hh-1L));
          }
          verbose && enter(verbose, sprintf("Atomic island #%d ('%s') of %d", kk, atomicIslandTag, length(atomicIslands)));
          n0 <- nbrOfSegments(fitT);
          verbose && cat(verbose, "Number of segments before: ", n0);
          fitTT <- fitT;

          if (hh > 0) {
            verbose && enter(verbose, sprintf("Dropping %d segments (%s)", hh, atomicIslandTag));
            fitTT <- dropRegions(fitTT, regions=atomicIsland, H=hh);
            nT <- nbrOfSegments(fitTT);
            verbose && exit(verbose);
            verbose && cat(verbose, "Number of segments left: ", nT);
            # Sanity check
            stopifnot(n0-nT == hh);
          }

          fitDrop <- fitTT$dropped;
          fitTT$dropped <- NULL;
          dropList[[kk]] <- fitDrop;

          verbose && enter(verbose, sprintf("Merging segments (#%d and #%d)", atomicIsland-1L, atomicIsland));
          fitTT <- mergeTwoSegments(fitTT, left=atomicIsland-1L);

          if (is.function(onAtomicIsland)) {
            onAtomicIsland(fit0=fitT, fit1=fitTT, fitD=fitDrop,
                           atomicIsland=atomicIsland, H=hh);
          }
          n1 <- nbrOfSegments(fitTT);
          verbose && exit(verbose);
          verbose && cat(verbose, "Number of segments left: ", n1);
          # Sanity check
          stopifnot(n0-n1 == hh+1);

          fitT <- fitTT;

          verbose && exit(verbose);
        } # for (kk ...)

        fitT$dropped <- dropList;
        fitT$atomicIslands <- rev(atomicIslands);
        fitT$H <- hh;

        hasChanged <- TRUE;

        verbose && exit(verbose);

        # Go to next generation
        break;
      } # if (length(atomicIslands) > 0)

      verbose && exit(verbose);
    } # for (hh ...)

    if (!hasChanged) {
      break;
    }

    # Next generation
    gg <- gg + 1L;

    verbose && exit(verbose);
  } # while(gg <= maxGeneration)

  class(fitList) <- c("PruneCNA", class(fitList));

  fitList;
}) # prunceCNA()



############################################################################
# HISTORY:
# 2012-06-05
# o Now pruneCNA() is for AbstractCBS, not just PairedPSCBS objects.
# 2011-01-18
# o Added class 'PruneCNA' to the return object of pruneCNA().
# o Now pruneCNA() returns pruned objects with the dropped segments
#   included in a separate list.
# o Added argument 'maxGeneration' to pruneCNA().
# 2010-09-08
# o Updated to use findAtomicAberrations().
# 2010-07-24
# o CLEAN UP: Now the notation of the code better reflect the algorithm.
# o Now findAtomicRegions() returns ambigous atomic regions too.
# o Added argument 'ylim'.
# 2010-07-20
# o Added argument 'debugPlot'.
# 2010-07-19
# o Added trial version of segmentByPruneCBS().
# o TO DO: Down-weight loci that were close to earlier
#   change points in the succeeding segmentations.
# o Added prototype version of findAtomicRegions().
# o Added prototype version of callByPruning().
# o Created.
############################################################################
