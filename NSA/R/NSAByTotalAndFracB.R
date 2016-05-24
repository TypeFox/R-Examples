###########################################################################/**
# @set "class=matrix"
# @RdocMethod NSAByTotalAndFracB
# @alias NSAByTotalAndFracB
# 
# @title "Finds normal regions within tumoral samples (total,fracB)"
#
# \description{
#  @get "title", where total is the total (non-polymorphic) signal and
#  fracB is the allele B fraction.
#  It is only loci with a non-missing (@NA) fracB value that are
#  considered to calculate the normal regions.
# }
#
# @synopsis
#
# \arguments{
#  \item{data}{An Jx2 @matrix, where J is the number of loci and
#                      2 is total and fracB.}
#  \item{...}{Additional arguments passed to fitNSA().}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns an Jx2 @numeric @array.
# }
#
#
#*/###########################################################################
setMethodS3("NSAByTotalAndFracB", "matrix", function(data, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'data':
  if (!is.matrix(data)) {
    throw("Argument 'data' is not an array: ", class(data)[1]);
  }
  dims <- dim(data);
  dimnames <- dimnames(data);
  if (length(dims) != 2) {
    throw("Argument 'data' is not a 2-dimensional array: ", 
                                                paste(dims, collapse="x"));
  }
  if (dims[2] != 2) {
    throw("Argument 'data' is not a Jx2-dimensional array: ", 
                                                paste(dims, collapse="x"));
  }
  if (!is.null(dimnames[[2]])) {
    if (dimnames[[2]][2]=="freqB"){
      dimnames(data)[[2]][2] = "fracB";  
      dimnames[[2]][2]="fracB";
    }    
    if (!identical(dimnames[[2]], c("total", "fracB"))) {
      throw("If given, the names of the allele (2nd) dimension of the Jx2-dimensional array (argument 'data') have to be 'total' & 'fracB': ", paste(dimnames[[2]], collapse=", "));
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  verbose && enter(verbose, "NSAByTotalAndFracB()");
  verbose && cat(verbose, "(total,fracB) signals:");
  verbose && str(verbose, data);

  # Number of elements in data
  nbrLoci <- nrow(data);

  verbose && enter(verbose, "Identifying SNPs (non-missing fracB)");
  nok <- is.na(data[,"fracB"]);
  snps <- which(!nok);
  verbose && printf(verbose, "Number of SNPs: %d (%.2f%%)\n",
                            length(snps), 100*length(snps)/dim(data)[1]);
  verbose && exit(verbose);
  dataC <- data;
  dataSNPs <- data[snps,];
  dataC[snps,] <- fitNSA(data[snps,],  ..., verbose=verbose);

  verbose && enter(verbose, "Identifying CN probes in normal regions");
  nok <- which(nok);

  if(length(nok)>0){
    dataC[nok,] <- NA;
    dataC <- fitNSAcnPs(dataC,..., verbose=verbose);
  }

  verbose && exit(verbose);

  rm(snps,data,nok); # Not needed anymore
  
  verbose && cat(verbose, "Calibrated (total,fracB) signals:");
  verbose && str(verbose, dataC);

  verbose && exit(verbose);

  dataC;
}) # NSAByTotalAndFracB()


###########################################################################
# HISTORY:
# 2010-06-29 [MO]
# o Created from calmateByTotalAndFracB().
###########################################################################
