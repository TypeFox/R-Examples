setMethodS3("segmentByCBS", "CNA", function(y, ...) {
  # To please R CMD check
  cnData <- y;

  # Extract signals of interest
  chromosome <- cnData$chrom;
  x <- cnData$maploc;
  y <- cnData[,3];
  signalType <- attr(cnData, "data.type");
  sampleName <- colnames(cnData)[3];

# str(list(y=y, chromosome=chromosome, x=x));

  fit <- segmentByCBS(y=y, chromosome=chromosome, x=x, ...);
  sampleName(fit) <- sampleName;

  fit;
}) # segmentByCBS()


#############################################################################
# HISTORY:
# 2011-09-04
# o Added segmentByCBS for CNA objects.
# o Created.
#############################################################################
