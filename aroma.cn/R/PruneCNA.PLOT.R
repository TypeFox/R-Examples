setMethodS3("plotTracks", "PruneCNA", function(this, nrow=length(this), ncol=1, byrow=FALSE, changepoints=TRUE, cex=1, col="#33cc66", xScale=1e-6, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'nrow' & 'ncol':
  if (!is.null(nrow)) {
    nrow <- Arguments$getIndex(nrow);
  }
  if (!is.null(ncol)) {
    ncol <- Arguments$getIndex(ncol);
  }

  # Argument 'byrow':
  byrow <- Arguments$getLogical(byrow);

  # Argument 'xScale':
  xScale <- Arguments$getDouble(xScale, range=c(0,Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # AD HOC
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fit <- this[[1]];
  if (inherits(fit, "PairedPSCBS")) {
    # Plot function
    plotFcn <- function(..., changepoints=TRUE, add=FALSE) {
      plotTracks2(..., changepoints=changepoints);
    }
    startKey <- "tcnStart";
    endKey <- "tcnEnd";
  } else if (inherits(fit, "CBS")) {
    plotFcn <- function(..., changepoints=TRUE, add=FALSE) {
      if (!add) {
        plotTracks(...);
      }
    }
    startKey <- "start";
    endKey <- "end";
  } else {
    throw("Non-supported segmentation result: ", class(fit)[1]);
  }


  # Setup panels  
  if (!is.null(nrow) && !is.null(nrow)) {
    subplots(nrow*ncol, nrow=nrow, ncol=ncol, byrow=byrow);
    par(mar=c(3,5,1,1));
  }

  nbrOfGenerations <- length(this);
  for (kk in seq_len(nbrOfGenerations)) { 
    verbose && enter(verbose, sprintf("Generation %d of %d", kk, nbrOfGenerations));

    fit <- this[[kk]];
    msg <- sprintf("Number of segments: %s", nbrOfSegments(fit));

    plotFcn(fit, ..., cex=cex, changepoints=changepoints, xScale=xScale);

    if (kk < nbrOfGenerations) {
      fitN <- this[[kk+1]];
      H <- fitN$H;
      verbose && cat(verbose, "H: ", H);

      nbrOfDrops <- length(fitN$atomicIslands);
      verbose && cat(verbose, "Number of dropped blocks: ", nbrOfDrops);

      # Anything dropped?
      if (!is.null(H)) {
        dropped <- fitN$dropped;
        if (H == 0) {
          verbose && enter(verbose, sprintf("Highlighting %d dropped change points (H=0)", nbrOfDrops));
          msg <- sprintf("%s. Dropping %d change points (H=0).", msg, nbrOfDrops);
          if (changepoints) {
            dropped <- extractRegions(fit, regions=fitN$atomicIslands);
            segs <- getSegments(dropped, splitters=TRUE);
            x <- segs[[startKey]];
            x <- x * xScale;
            abline(v=x, lwd=1, col=col);
            verbose && cat(verbose, "Change-point locations:");
            verbose && print(verbose, as.matrix(x));
          }

          verbose && exit(verbose);
        } else if (H == 1) {
          verbose && enter(verbose, sprintf("Highlighting %d dropped single-segment (H=1) blocks", nbrOfDrops));
          msg <- sprintf("%s. Dropping %s single-segment (H=1) blocks.", msg, nbrOfDrops);
          verbose && exit(verbose);
        } else {
          verbose && enter(verbose, sprintf("Highlighting %d dropped blocks.", nbrOfDrops));
          msg <- sprintf("%s. Dropping %s segments from %d blocks each containing H=%d segments", msg, H*nbrOfDrops, nbrOfDrops, H);
          verbose && exit(verbose);
        }
  
        if (H > 0) {
          verbose && cat(verbose, "Segments:");

          for (fd in dropped) { 
            plotFcn(fd, ..., col=col, cex=2*cex, xScale=xScale, 
                       changepoints=FALSE, add=TRUE);
          }
  
          # Highlight segments means again
##          plotFcn(fit, ..., changepoints=FALSE, xScale=xScale, add=TRUE);

          if (changepoints) {
            x <- sapply(dropped, FUN=function(fd) {
              segs <- getSegments(fd, splitters=TRUE);
              c(segs[[startKey]][1], segs[[endKey]][nrow(segs)]);
            });
            x <- t(x);
            colnames(x) <- c("blockStart", "blockEnd");
            x <- x[order(x[,1]),,drop=FALSE];
            x <- cbind(x, length=x[,2]-x[,1]);
            x <- x * xScale;
            verbose && print(verbose, x);
            abline(v=x[,1:2], lwd=1, col=col);
            yy <- par("usr")[3:4];
            for (rr in seq_len(nrow(x))) {
              xs <- x[rr,1:2];
#str(list(xs=xs, ys=ys));
              ys <- yy[1]*c(1,1);
              lines(x=xs, y=ys, lwd=4, col=col);
              ys <- yy[2]*c(1,1);
              lines(x=xs, y=ys, lwd=4, col=col);
            }
#            box();
          }
        } # if (H > 0)
      } # if (!is.null(dropped))
    } else {
      msg <- sprintf("%s.", msg);
    } # if (kk < nbrOfGenerations)

    stext(side=3, pos=0, msg);

    verbose && exit(verbose);
  } # for (kk ...)
}) # plotTracks()


############################################################################
# HISTORY:
# 2012-06-05 [HB]
# o Now plotTracks() for PruneCNA supports CBS segmentation results
#   in additional to PairedPSCBS ones.
# 2011-07-10 [HB]
# o Updated code to work with the new column names in PSCBS v0.11.0.
# 2011-01-19
# o Now plotTracks() of PruneCNA utilizes plotTracks2() for PairedPSCBS.
# o Dropped argument 'tracks'.
# 2011-01-18
# o Added plotTracks() for PruneCNA.
############################################################################
