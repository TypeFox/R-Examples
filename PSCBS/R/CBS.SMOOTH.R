###########################################################################/**
# @set class=CBS
# @RdocMethod getSmoothLocusData
#
# @title "Gets smoothed locus-level data"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{fit}{An @see "CBS" object.}
#   \item{by}{A @numeric scalar specifying the bin size.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @data.frame where the
#   first three columns are 'chromosome', 'x' (position),
#   and 'count' (number of loci average over for the given bin),
#   and the remaining ones are the smoothed locus-level data.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("getSmoothLocusData", "CBS", function(fit, by, ...) {
  # Argument 'by':
  by <- Arguments$getNumeric(by, range=c(0,Inf));

  chromosomes <- getChromosomes(fit);
  data <- getLocusData(fit);

  chromosome <- NULL; rm(list="chromosome"); # To please R CMD check

  dataS <- NULL;
  for (kk in seq_along(chromosomes)) {
    chr <- chromosomes[kk];
    dataT <- subset(data, chromosome == chr);
    x <- dataT$x;
    y <- dataT$y;
    rx <- range(x, na.rm=TRUE);
    bx <- seq(from=rx[1], to=rx[2], by=by);
    xS <- bx[-1] - by/2;
    yS <- binMeans(y=y, x=x, bx=bx);
    count <- attr(yS, "count");
    yS[count == 0L] <- NA_real_;
    attr(yS, "count") <- NULL;
    dataTS <- data.frame(chromosome=chr, x=xS, count=count, y=yS);
    dataS <- rbind(dataS, dataTS);
  } # for (kk ...)

  dataS;
}, protected=TRUE) # getSmoothLocusData()


############################################################################
# HISTORY:
# 2013-10-09
# o Now getSmoothLocusData() for CBS also returns column 'count'.
# 2013-04-18
# o Added getSmoothLocusData() for CBS.
# o Created.
############################################################################
