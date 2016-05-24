setMethodS3("drawCnRegions", "HaarSeg", function(object, ...) {
  cnr <- extractCopyNumberRegions(object, ...);
  drawLevels(cnr, ...);
})


setMethodS3("extractCopyNumberRegions", "HaarSeg", function(object, ...) {
  output <- object$output;
  regions <- output$SegmentsTable;
  data <- object$data;
  x <- data$x;

  # Identify the locus indices where the regions starts and ends
  starts <- regions[,1];
  counts <- regions[,2];
  ends <- starts+counts-1;

  # Translate to genomic positions
  starts <- x[starts];
  ends <- x[ends];

  # Get the mean levels of each region
  means <- regions[,3];

  CopyNumberRegions(
    chromosome=data$chromosome,
    start=starts, 
    stop=ends, 
    mean=means,
    count=counts
  );
})


setMethodS3("extractRawCopyNumbers", "HaarSeg", function(object, ...) {
  data <- object$data;
  RawCopyNumbers(cn=data$M, x=data$x, chromosome=data$chromosome);
})


##############################################################################
# HISTORY:
# 2009-05-14
# o Moved extractRawCopyNumbers() for profileCGH from aroma.affymetrix.
# 2009-05-10
# o Moved to aroma.core v1.0.6. Source file: HaarSegModel.R.
# 2008-01-26
# o Reordered constructor arguments.
# 2008-12-31
# o Removing non-finite data points before passing to haarSeg().
# 2008-12-17
# o Now using the HaarSeg package (put together by HB).
# 2008-12-16
# o Created.
############################################################################## 
