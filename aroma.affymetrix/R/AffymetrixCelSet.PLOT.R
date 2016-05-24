###########################################################################/**
# @set "class=AffymetrixCelSet"
# @RdocMethod plotDensity
#
# @title "Plots the densities of all samples"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{subset}{The subset of probes to considered \emph{before} any
#     filtering by probe type is applied.
#     If a @vector of @doubles, the cell indices.
#     If a scalar @double in [0,1], the fraction of cells, which can
#     be used to speed up the plotting if approximate densities are
#     acceptable.
#     if @NULL, all cells are considered.
#   }
#   \item{types}{The type of probes to include, e.g. \code{"all"},
#     \code{"pmmm"}, \code{"pm"}, and \code{"mm"}.}
#   \item{...}{Additional arguments passed to
#      @see "plotDensity.AffymetrixCelFile".}
#   \item{col}{A @vector of colors for each of the arrays.}
#   \item{lty}{A @vector of line types for each of the arrays.}
#   \item{lwd}{A @vector of line widths for each of the arrays.}
#   \item{add}{If @FALSE, a new plot is created, otherwise the generated
#     graphics is added to the current plot.}
# }
#
# \value{
#  Returns nothing.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("plotDensity", "AffymetrixCelSet", function(this, subset=NULL, types=NULL, ..., col=seq_along(this), lty=NULL, lwd=NULL, annotate=TRUE, add=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'subset':
  if (is.null(subset)) {
  } else if (is.numeric(subset)) {
  }

  nbrOfArrays <- length(this);

  # Argument 'col':
  if (is.null(col)) {
    col <- seq_len(nbrOfArrays);
  } else {
    col <- rep(col, length.out=nbrOfArrays);
  }

  # Argument 'lty':
  if (!is.null(lty))
    lty <- rep(lty, length.out=nbrOfArrays);

  # Argument 'lwd':
  if (!is.null(lwd))
    lwd <- rep(lwd, length.out=nbrOfArrays);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify the subset of probes to be updated
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);
  verbose && enter(verbose, "Identifying subset of probes");
  suppressWarnings({
    subset <- identifyCells(cdf, indices=subset, types=types,
                                                    verbose=less(verbose));
  })
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Plot densities
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (kk in seq_len(nbrOfArrays)) {
    df <- this[[kk]];
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", kk,
                                                getName(df), nbrOfArrays));

    verbose && cat(verbose, "Using cell indices:");
    verbose && str(verbose, subset);

    add <- add || (kk > 1);
    plotDensity(df, subset=subset, ..., col=col[kk], lty=lty[kk],
              lwd=lwd[kk], annotate=FALSE, add=add, verbose=less(verbose));

    if (annotate) {
      stextChipType(getChipType(this));
      stextSize(df, size=length(subset));
      annotate <- FALSE;
    }

    verbose && exit(verbose);
  }
})




############################################################################
# HISTORY:
# 2009-09-17
# o Now argument 'subset' of plotDensity() of AffymetrixCelSet defaults
#   to NULL (all probes).  Before it was 1/2 (a fraction).
# o Added verbose output to plotDensity() to AffymetrixCelSet show what
#   cells are used.
# 2007-04-16
# o Added more verbose output.
# 2006-05-16
# o Created.
############################################################################
