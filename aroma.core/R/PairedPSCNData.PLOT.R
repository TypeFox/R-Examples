###########################################################################/**
# @set "class=PairedPSCNData"
# @RdocMethod plotTracks
# @alias plotTracks
#
# @title "Plots parental specific copy numbers along the genome"
#
# \description{
#  @get "title" for one or more chromosomes.
#  It is possible to specify what type of tracks to plot.
#  Each type of track is plotted in its own panel.
# }
#
# @synopsis
#
# \arguments{
#   \item{tracks}{A @character @vector specifying what types of tracks to plot.}
#   \item{pch}{The type of the scatter points, if any.}
#   \item{col}{The color of the scatter points, if any.}
#   \item{cex}{The size of the scatter points, if any.}
#   \item{grid}{If @TRUE, horizontal lines are displayed.}
#   \item{xlim}{(Optional) The genomic range to plot.}
#   \item{Clim}{The range of copy numbers.}
#   \item{Blim}{The range of allele B fractions (BAFs) and 
#     decrease of heterozygosity (DHs).}
#   \item{xScale}{The scale factor used for genomic positions.}
#   \item{...}{Not used.}
#   \item{add}{If @TRUE, the panels plotted are added to the existing plot,
#     otherwise a new plot is created.}
#   \item{subplots}{If @TRUE, then subplots are automatically setup.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns nothing.
# }
# 
# @author
#
# @keyword IO
# @keyword internal
#*/########################################################################### 
setMethodS3("plotTracks", "PairedPSCNData", function(x, tracks=c("tcn", "dh", "tcn,c1,c2", "tcn,c1", "tcn,c2", "c1,c2", "betaN", "betaT", "betaTN")[1:3], pch=".", col=NULL, cex=1, grid=FALSE, xlim=NULL, Clim=c(0,6), Blim=c(0,1), xScale=1e-6, ..., add=FALSE, subplots=!add && (length(tracks) > 1), verbose=FALSE) {

  # To please R CMD check
  this <- x;
 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'this':
  if (nbrOfChromosomes(this) > 1) {
    return(plotTracksManyChromosomes(this, tracks=tracks, pch=pch, Clim=Clim, Blim=Blim, xScale=xScale, ..., add=add, verbose=verbose));
  }

  # Argument 'tracks':
  knownTracks <- c("tcn", "dh", "tcn,c1,c2", "tcn,c1", "tcn,c2", "c1,c2", "betaN", "betaT", "betaTN");
  tracks <- match.arg(tracks, choices=knownTracks, several.ok=TRUE);
  tracks <- unique(tracks);

  # Argument 'grid':
  grid <- Arguments$getLogical(grid);

  # Argument 'xlim':
  if (!is.null(xlim)) {
    xlim <- Arguments$getNumerics(xlim, length=c(2,2));
  }

  # Argument 'xScale':
  xScale <- Arguments$getNumeric(xScale, range=c(0,Inf));

  # Argument 'add':
  add <- Arguments$getLogical(add);

  # Argument 'subplots':
  subplots <- Arguments$getLogical(subplots);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Plotting ", class(this)[1], " tracks");
  stopifnot(nbrOfChromosomes(this) == 1);

  # Extract the locus-level signals (including virtual ones)
  data <- as.data.frame(this);

  chromosome <- data$chromosome[1];
  verbose && cat(verbose, "Chromosome: ", chromosome);

  nbrOfLoci <- nrow(data);
  verbose && cat(verbose, "Number of loci: ", nbrOfLoci);

  x <- data$x;
  CT <- data$CT;  # <== data$C?!? /HB 2012-04-16
  betaT <- data$betaT;
  betaN <- data$betaN;
  betaTN <- data$betaTN;
  muN <- data$muN;

  if (chromosome != 0) {
    chrTag <- sprintf("Chr%02d", chromosome);
  } else {
    chrTag <- "";
  }

  if (xScale != 1) {
    x <- xScale * x;
    if (!is.null(xlim)) {
      xlim <- xScale * xlim;
    }
  }

  if (subplots) {
    subplots(length(tracks), ncol=1);
    par(mar=c(1,4,1,2)+1);
  }

  # Color loci by genotype
  colMu <- c("gray", "black")[(muN == 1/2) + 1];

  for (tt in seq_along(tracks)) {
    track <- tracks[tt];
    verbose && enter(verbose, sprintf("Track #%d ('%s') of %d", 
                                             tt, track, length(tracks)));

    pchT <- pch;
    colT <- col;

    if (track == "tcn") {
      colT <- ifelse(is.null(colT), "black", colT);
      if (add) {
        points(x, CT, pch=pchT, col=colT, cex=cex);
      } else {
        plot(x, CT, pch=pchT, col=colT, cex=cex, xlim=xlim, ylim=Clim, ylab="TCN");
        stext(side=3, pos=1, chrTag);
        if (grid) {
          abline(h=seq(from=0, to=Clim[2], by=2), lty=3, col="gray");
          abline(h=0, lty=1, col="black");
        }
      }
    }
  
    if (is.element(track, c("tcn,c1,c2", "tcn,c1", "tcn,c2", "c1,c2"))) {
      colT <- ifelse(is.null(colT), "black", colT);
      subtracks <- strsplit(track, split=",", fixed=TRUE)[[1]];
      ylab <- paste(toupper(subtracks), collapse=", ");
      if (add) {
        points(x, CT, pch=pchT, cex=cex, col=colT);
      } else {
        plot(x, CT, pch=pchT, cex=cex, col=colT, xlim=xlim, ylim=Clim, ylab=ylab);
        stext(side=3, pos=1, chrTag);
        if (grid) {
          abline(h=seq(from=0, to=Clim[2], by=2), lty=3, col="gray");
          abline(h=0, lty=1, col="black");
        }
      }
    }
  
    if (track == "betaN") {
      colT <- ifelse(is.null(colT), colMu, colT);
      if (add) {
        points(x, betaN, pch=pchT, cex=cex, col="black");
      } else {
        plot(x, betaN, pch=pchT, cex=cex, col=colT, xlim=xlim, ylim=Blim, ylab=expression(BAF[N]));
        stext(side=3, pos=1, chrTag);
      }
    }
  
    if (track == "betaT") {
      colT <- ifelse(is.null(colT), colMu, colT);
      if (add) {
        points(x, betaT, pch=pchT, cex=cex, col="black");
      } else {
        plot(x, betaT, pch=pchT, cex=cex, col=colT, xlim=xlim, ylim=Blim, ylab=expression(BAF[T]));
        stext(side=3, pos=1, chrTag);
      }
    }
  
    if (track == "betaTN") {
      colT <- ifelse(is.null(colT), colMu, colT);
      if (add) {
        points(x, betaTN, pch=pchT, cex=cex, col="black");
      } else {
        plot(x, betaTN, pch=pchT, cex=cex, col=colT, xlim=xlim, ylim=Blim, ylab=expression(BAF[T]^"*"));
        stext(side=3, pos=1, chrTag);
      }
    }
  
    if (track == "dh") {
      isSnp <- (!is.na(betaTN) & !is.na(muN));
      isHet <- isSnp & (muN == 1/2);
      naValue <- as.double(NA);
      rho <- rep(naValue, length=nbrOfLoci);
      rho[isHet] <- 2*abs(betaTN[isHet]-1/2);
      colT <- ifelse(is.null(colT), colMu[isHet], colT);
      if (add) {
        points(x, rho, pch=pchT, cex=cex, col="black");
      } else {
        plot(x, rho, pch=pchT, cex=cex, col=colT, xlim=xlim, ylim=Blim, ylab="DH");
        stext(side=3, pos=1, chrTag);
      }
    }

    verbose && exit(verbose);
  } # for (tt ...)

  verbose && exit(verbose);

  invisible();  
}) # plotTracks()
 

############################################################################
# HISTORY:
# 2012-04-16
# o Added plotTracks() for PairedPSCNData; adopted from ditto for 
#   PairedPSCBS of the PSCBS package.
# o Created.
############################################################################
