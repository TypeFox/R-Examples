###########################################################################/**
# @set "class=matrix"
# @RdocMethod snpsNByTotalAndFracB
# @alias snpsNByTotalAndFracB
# 
# @title "Normalize allele-specific copy numbers (total,fracB)"
#
# \description{
#  @get "title", where total is the total (non-polymorphic) signal and
#  fracB is the allele B fraction. It normalizes the data snp by snp.
# }
#
# @synopsis
#
# \arguments{
#  \item{data}{An JxI @matrix containing the copy number values, where J is the number of loci,
#                      and I is the number of samples.}
#  \item{references}{A @logical or @numeric @matrix specifying which
#     samples should be used as the reference set.  
#     By default, all samples are considered.}
#  \item{...}{Optional argument passed to fitSample(). This argument "Threshold" is used to determine if
#     one region is considered normal or not. Initially set to .135.}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns an JxI @numeric @matrix.
# }
#*/###########################################################################
setMethodS3("snpsNByTotalAndFracB", "matrix", function(data, references = NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'data':
  if (!is.array(data)) {
    throw("Argument 'data' is not an matrix: ", class(data)[1]);
  }
  dim <- dim(data);
  dimnames <- dimnames(data);
  if (length(dim) != 2) {
    throw("Argument 'data' is not a matrix: ", 
                                                paste(dim, collapse="x"));
  }
  
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  verbose && enter(verbose, "snpsNByTotalAndFracB()");
  verbose && cat(verbose, "(total,fracB) signals:");
  verbose && str(verbose, data);

  verbose && enter(verbose, "Identifying SNPs (non-finite values)");

  nok <- is.na(data);
  nok <- rowAlls(nok);
  loci <- which(!nok);
  dataS <- data[loci,,drop=FALSE];

  refS <- references[loci,];
  verbose && printf(verbose, "Number of loci: %d (%.2f%%)\n",
                            length(loci), 100*length(loci)/dim(data)[1]);
  verbose && exit(verbose);

  verbose && enter(verbose, "Normalizing the data by loci");  

  nbrOfLoci <- dim(dataS)[1];
  verbose && cat(verbose, "Number of loci: ", nbrOfLoci);
  verbose && printf(verbose, "Number of loci left:");
  
  if(sum(nok) < nrow(data)){
    fit <- fitSNPsN(dataS, references = refS, ...);
    stopifnot(identical(dim(fit), dim(dataS)));  
    verbose && cat(verbose, "Normalized by loci (total) signals:");
    verbose && str(verbose, fit);
    data[loci,] <- fit;
    rm(fit); #not needed
  }
  rm(loci, dataS);
  
#  verbose && enter(verbose, "Calibrating non-polymorphic probes");
#  dataC[nok,"total",] <- fitCalMaTeCNprobes(data[nok,"total",], references=references);
#  verbose && str(verbose, dataC);
#  verbose && exit(verbose);

  verbose && exit(verbose);

  data;
}) # snpsNByTotalAndFracB()


###########################################################################
# HISTORY:
# 2010-06-28 [MO]
# o Created from NSAByTotalAndFracB.R.
###########################################################################
