###########################################################################/**
# @set "class=PSCBS"
# @RdocMethod append
#
# @title "Appends one segmentation result to another"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{x, other}{The two @see "PSCBS" objects to be combined.}
#  \item{other}{A @see "PSCBS" object.}
#  \item{addSplit}{If @TRUE, a "divider" is added between chromosomes.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @see "PSCBS" object of the same class as argument \code{x}.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("append", "PSCBS", function(x, other, addSplit=TRUE, ...) {
  # To please R CMD check
  this <- x;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'other':
  other <- Arguments$getInstanceOf(other, "PSCBS");
  for (field in c("data", "output")) {
    dataA <- this[[field]]
    dataB <- other[[field]]
    namesA <- colnames(dataA)
    namesB <- colnames(dataB)
    if (!all(namesA == namesB)) {
      throw(sprintf("Cannot merge %s objects. Arguments 'other' and 'this' has different sets of columns in field '%s': {%s} [n=%d] != {%s} [n=%d]", class(this)[1], field, paste(namesA, collapse=", "), length(namesA), paste(namesB, collapse=", "), length(namesB)))
    }
  }

  # Argument 'addSplit':
  addSplit <- Arguments$getLogical(addSplit);


  # Allocate results
  res <- this;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locus data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res$data <- rbind(this$data, other$data);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Segmentation data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  indexOffset <- nrow(this$data);
  fields <- c("output", "tcnSegRows", "dhSegRows");
  for (field in fields[-1]) {
    other[[field]] <- other[[field]] + indexOffset;
  }

  splitter <- if (addSplit) NA else NULL;
  for (field in fields) {
    res[[field]] <- rbind(this[[field]], splitter, other[[field]]);
    rownames(res[[field]]) <- NULL;
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ksT <- this$params$knownSegments;
  ksT$length <- NULL;  # In case it's been added
  ksO <- other$params$knownSegments;
  ksO$length <- NULL;  # In case it's been added
  res$params$knownSegments <- rbind(ksT, ksO);

  # Sanity check
  ns <- sapply(res[fields], FUN=nrow);
  stopifnot(all(ns == ns[1]));

  res;
}) # append()


setMethodS3("extractChromosomes", "PSCBS", function(x, chromosomes, ...) {
  # To please R CMD check
  this <- x;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'chromosomes':
  disallow <- c("NaN", "Inf");
  chromosomes <- Arguments$getIntegers(chromosomes, range=c(0,Inf), disallow=disallow);
  stopifnot(all(is.element(chromosomes, getChromosomes(this))));

  # Always extract in order
  chromosomes <- unique(chromosomes);
  chromosomes <- sort(chromosomes);

  # Allocate results
  res <- this;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locus data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chromosome <- NULL; rm(list="chromosome"); # To please R CMD check
  res$data <- subset(res$data, chromosome %in% chromosomes);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Segmentation data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify rows to subset
  rows <- which(is.element(res$output$chromosome, chromosomes));
  for (field in c("output", "tcnSegRows", "dhSegRows")) {
    res[[field]] <- res[[field]][rows,,drop=FALSE];
  }

  # Identify chromosome offsets
  chrStarts <- match(getChromosomes(this), this$data$chromosome);
  chrEnds <- c(chrStarts[-1]-1L, nrow(this$data));
  chrLengths <- chrEnds - chrStarts + 1L;

  chrLengthsExcl <- chrLengths;

  keep <- match(chromosomes, getChromosomes(this));
  chrLengthsExcl[keep] <- 0L;
  cumChrLengthsExcl <- cumsum(chrLengthsExcl);

  shifts <- cumChrLengthsExcl[keep];
  stopifnot(all(is.finite(shifts)));

  # Adjust indices
  for (cc in seq(along=chromosomes)) {
    chromosome <- chromosomes[cc];
    shift <- shifts[cc];
    # Nothing to do?
    if (shift == 0) next;
    for (field in c("tcnSegRows", "dhSegRows")) {
      segRows <- res[[field]];
      rows <- which(res$output$chromosome == chromosome);
      segRows[rows,] <- segRows[rows,] - shift;
      res[[field]] <- segRows;
    }
  }

  res;
}, protected=TRUE)



############################################################################
# HISTORY:
# 2012-09-21
# o ROBUSTNESS: Now append() for CBS and PSCBS drops column 'length'
#   from 'knownSegments', iff it exists.
# 2011-10-20
# o Now append() for PSCBS also appends '...$params$knownSegments'.
# 2011-10-02
# o Now the CBS class extends the AbstractCBS class.
# o Added print() and as.data.frame() to PSCBS.
# o Added getSegments() to PSCBS.
# o DOCUMENTATION: Added Rdoc for several PSCBS methods.
# o Added a PSCBS constructor with documentation.
# 2010-12-01
# o Now also extractByChromosomes() and append() for PSCBS recognizes
#   fields 'tcnLociToExclude' and 'dhLociToExclude'.
# o BUG FIX: extractByChromosome() for PSCBS would call it self instead
#   of extractByChromosomes().
# 2010-11-26
# o Added extractByChromosomes() for PSCBS.
# 2010-09-26
# o getChromosomes() no longer returns NA divers.
# 2010-09-24
# o Added append() and more for PSCBS objects.
############################################################################
