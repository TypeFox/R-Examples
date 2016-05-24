setMethodS3("readCfhUnits", "default", function(pathname, snps=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pathname':
  pathname <- Arguments$getReadablePathname(pathname, mustExist=TRUE);

  # Argument 'units':
  if (!is.null(snps))
    snps <- Arguments$getIndices(snps);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "Reading data from CFH file");
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieve file header
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  hdr <- readCfhHeader(pathname);
  verbose && cat(verbose, "Header:");
  verbose && str(verbose, hdr);
  nbrOfBytes <- hdr$nbrOfBytes;
  nbrOfSnps <- hdr$nbrOfSnps;
  bytesPerSnp <- hdr$bytesPerSnp;

  # Validating units
  if (!is.null(snps))
    snps <- Arguments$getIndices(snps, max=nbrOfSnps);

  map <- matrix(1:(bytesPerSnp*nbrOfSnps), nrow=bytesPerSnp);
  map <- map + hdr$dataOffset;

  # Read subset of SNPs?
  if (!is.null(snps))
    map <- map[,snps,drop=FALSE];
  

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieve data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read all data
  raw <- readBin(pathname, what=raw(), n=nbrOfBytes);

  # Bytes 1:8 contains (thetaA,thetaB) as floats
  rr <- c(1:8,10:13);
  ncol <- length(rr) / 4;
  verbose && cat(verbose, "Number of integers: ", ncol);

  map <- map[rr,,drop=FALSE];
  theta <- readBin(raw[map], what=double(), size=4, endian="little", 
                                                           n=ncol*ncol(map));
  theta <- matrix(theta, ncol=ncol, byrow=TRUE);

#  colnames(theta) <- c("A", "B");

  verbose && exit(verbose);
  
  theta;
})


############################################################################
# HISTORY:
# 2007-04-06
# o Created.
############################################################################
