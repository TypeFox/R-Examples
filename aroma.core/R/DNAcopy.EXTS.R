setMethodS3("extractCopyNumberRegions", "DNAcopy", function(object, ...) {
  output <- object$output;

  CopyNumberRegions(
    chromosome=output[["chrom"]], 
    start=output[["loc.start"]], 
    stop=output[["loc.end"]], 
    mean=output[["seg.mean"]],
    count=output[["num.mark"]]
  );
})


setMethodS3("extractRawCopyNumbers", "DNAcopy", function(object, ...) {
  data <- object$data;
  chromosome <- unique(data$chrom);
  chromosome <- Arguments$getIndex(chromosome);
  RawCopyNumbers(cn=data[[3]], x=data$maploc, chromosome=chromosome);
})
  

setMethodS3("drawCnRegions", "DNAcopy", function(this, ...) {
  cnr <- extractCopyNumberRegions(this, ...);
  drawLevels(cnr, ...);
})


##############################################################################
# HISTORY:
# 2009-05-14
# o Moved extractRawCopyNumbers() for DNAcopy from aroma.affymetrix.
# o Moved extractCopyNumberRegions() for DNAcopy from aroma.affymetrix.
# 2009-05-10
# o Moved to aroma.core v1.0.6.  Source file: DNAcopy.drawCnRegions.R.
# 2008-05-21
# o Now extractRawCopyNumbers() adds 'chromosome' to the returned object.
# 2007-08-20
# o Created.
##############################################################################
