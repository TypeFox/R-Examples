setMethodS3("readCfnUnits", "default", function(pathname, cnagId=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pathname':
  pathname <- Arguments$getReadablePathname(pathname, mustExist=TRUE);

  # Argument 'units':
  if (!is.null(cnagId))
    cnagId <- Arguments$getIndices(cnagId);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "Reading data from CFN file");
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
  if (!is.null(cnagId)) {
    cnagId <- Arguments$getIndices(cnagId, range=hdr$cnagIdRange);
  }

  map <- matrix(1:(bytesPerSnp*nbrOfSnps), nrow=bytesPerSnp);
  map <- map + hdr$dataOffset - bytesPerSnp;

  # Read subset of SNPs?
  if (!is.null(cnagId)) {
    rr <- cnagId - hdr$cnagIdRange[1] + 1;
    map <- map[,rr,drop=FALSE];
  }
  

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieve data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read all data
  raw <- readBin(pathname, what=raw(), n=nbrOfBytes);

  # Bytes 13:16 contains (M) as floats
  rrM <- 13:16;
  # Bytes 17:20 contains (C) as integers
  rrC <- 17:20;
  rr <- 1:12;
rr <- rrM;
  ncol <- length(rr) / 4;
  map <- map[rr,,drop=FALSE];
  theta <- readBin(raw[map], what=double(), size=4, endian="little", 
#  theta <- readBin(raw[map], what=integer(), size=4, endian="little", 
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
