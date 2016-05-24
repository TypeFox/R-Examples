###########################################################################/**
# @set "class=numeric"
# @RdocMethod sampleNByTotalAndFracB
# @alias sampleNByTotalAndFracB
# 
# @title "Normalize total copy numbers by samples (total,fracB)"
#
# \description{
#  @get "title", where total is the total (non-polymorphic) signal and
#  fracB is the allele B fraction.
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
#  \item{...}{Optional argument passed to fitSNPsN(). This argument "Threshold" is used to determine if
#     one region is considered normal or not. Initially set to .135.}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns an Jx2 @numeric @array.
# }
#
#*/###########################################################################
setMethodS3("sampleNByTotalAndFracB", "numeric", function(data, references=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'data':
  if (!is.numeric(data)) {
    throw("Argument 'data' is not numeric: ", class(data)[1]);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "sampleNByTotalAndFracB()");
  verbose && cat(verbose, "('total' signals:");
  verbose && str(verbose, data);
  verbose && cat(verbose, "'normal regions' signals:");
  verbose && str(verbose, references);
  
  verbose && enter(verbose, "Identifying loci (non-missing/finite values)");
  nok <- is.na(data);
  loci <- which(!nok);
  verbose && printf(verbose, "Number of loci: %d (%.2f%%)\n",
                            length(loci), 100*length(loci)/dim(data)[1]);
  theta <- data[loci,drop=FALSE];
  refs <- references[loci];
  verbose && str(verbose, theta);
  verbose && str(verbose, refs);  
  verbose && exit(verbose);

  dataC <- data;
  verbose && enter(verbose, "Normalizing the signals.");
  dataC[loci] <- fitSample(theta, references = refs, ..., verbose=verbose);
  rm(data, theta, refs, loci); # Not needed anymore

  verbose && cat(verbose, "Calibrated (total,fracB) signals:");
  verbose && str(verbose, dataC);
  verbose && exit(verbose);
  
  verbose && exit(verbose);
  
  dataC;
}) # sampleNByTotalAndFracB()


###########################################################################
# HISTORY:
# 2010-06-29 [MO]
# o Created.
###########################################################################
