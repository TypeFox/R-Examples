###########################################################################/**
# @set "class=array"
# @RdocMethod calmateByTotalAndFracB
# @alias calmateByTotalAndFracB
# 
# @title "Normalize allele-specific copy numbers (total,fracB)"
#
# \description{
#  @get "title", where total is the total (non-polymorphic) signal and
#  fracB is the allele B fraction.
#  It is only loci with a non-missing (@NA) fracB value that are
#  considered to be SNPs and normalized by CalMaTe.  The other loci
#  are left untouched.
# }
#
# @synopsis
#
# \arguments{
#  \item{data}{An Jx2xI @numeric @array, where J is the number of loci,
#              2 is total and fracB (in that order, if unnamed), and 
#              I is the number of samples.}
#  \item{references}{A @logical or @numeric @vector specifying which
#     samples should be used as the reference set.  
#     By default, all samples are considered. If not NULL at least 3 samples.}
#  \item{...}{Additional arguments passed to @seemethod "calmateByThetaAB".}
#  \item{refAvgFcn}{(optional) A @function that takes a JxI @numeric @matrix
#     an argument \code{na.rm} and returns a @numeric @vector of length J.
#     It should calculate some type of average for each of the J rows, e.g.
#     @see "matrixStats::rowMedians".  
#     If specified, then the total copy numbers of the calibrated ASCNs
#     are standardized toward (twice) the average of the total copy numbers
#     of the calibrated reference ASCNs.}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns an Jx2xI @numeric @array
#   with the same dimension names as argument \code{data}.
# }
#
# @examples "../incl/calmateByTotalAndFracB.Rex"
#
# \references{
#  [1] @include "../incl/OrtizM_etal_2012.Rd" \cr 
# }
#
# \seealso{
#  To calibrate (thetaA,thetaB) or (CA,CB) signals, 
#  see @seemethod "calmateByThetaAB".
# }
#*/###########################################################################
setMethodS3("calmateByTotalAndFracB", "array", function(data, references=NULL, ..., refAvgFcn=NULL, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'data':
  if (!is.array(data)) {
    throw("Argument 'data' is not an array: ", class(data)[1]);
  }
  dim <- dim(data);
  dimnames <- dimnames(data);
  if (length(dim) != 3) {
    throw("Argument 'data' is not a 3-dimensional array: ", 
                                                paste(dim, collapse="x"));
  }
  if (dim[2] != 2) {
    throw("Argument 'data' is not a Jx2xI-dimensional array: ", 
                                                paste(dim, collapse="x"));
  }
  if (!is.null(dimnames[[2]])) {
    # Backward compatibility
    if (dimnames(data)[[2]][2] == "freqB"){
      dimnames(data)[[2]][2] <- "fracB";
    }
    if (!identical(dimnames(data)[[2]], c("total", "fracB"))) {
      throw("If given, the names of the allele (2nd) dimension of the Jx2xI-dimensional array (argument 'data') have to be 'total' & 'fracB': ", paste(dimnames(data)[[2]], collapse=", "));
    }
  }

  # From here on we force dimension names on the 2nd dimension
  dimnames(data)[[2]] <- c("total", "fracB");

  nbrOfSamples <- dim[3];
  if (nbrOfSamples < 3) {
    throw("Argument 'data' contains less than three samples: ", nbrOfSamples);
  }

  if (nbrOfSamples < 3) {
    throw("Argument 'data' contains less than three samples: ", nbrOfSamples);
  }

  # Argument 'references':
  if (is.null(references)) {
    # The default is that all samples are used to calculate the reference.
    references <- seq(length=nbrOfSamples);
  } else if (is.logical(references)) {
    if (length(references) != nbrOfSamples) {
      throw("Length of argument 'references' does not match the number of samples in argument 'data': ", length(references), " != ", nbrOfSamples);
    }
    references <- which(references);
  } else if (is.numeric(references)) {
    references <- as.integer(references);
    if (any(references < 1 | references > nbrOfSamples)) {
      throw(sprintf("Argument 'references' is out of range [1,%d]: %d", nbrOfSamples), length(references));
    }
  }
  if (length(references) < 3) {
    throw("Argument 'reference' specify less than three reference samples: ", length(references));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "calmateByTotalAndFracB()");
  verbose && cat(verbose, "(total,fracB) signals:");
  verbose && str(verbose, data);

  verbose && enter(verbose, "Identifying SNPs (non-missing fracB)");
  nok <- is.na(data[,"fracB",,drop=FALSE]);
  dim(nok) <- dim(nok)[-2]; # Drop the 2nd dimension
  nok <- rowAlls(nok);
##   save(nok, file="nok.Rdata");
  snps <- which(!nok);
  verbose && printf(verbose, "Number of SNPs: %d (%.2f%%)\n",
                            length(snps), 100*length(snps)/dim(data)[1]);
  verbose && exit(verbose);

  verbose && enter(verbose, "Transforming SNPs to (thetaA, thetaB)");
  theta <- data[snps,,,drop=FALSE];
  theta <- totalAndFracB2ThetaAB(theta, verbose=less(verbose, 5));
  verbose && str(verbose, theta);
  verbose && exit(verbose);

  thetaC <- calmateByThetaAB(theta, references=references, ..., verbose=verbose);  
  rm(theta); # Not needed anymore

  verbose && enter(verbose, "Backtransforming SNPs to (total, fracB)");
  dataC <- data;
  dataC[snps,,] <- thetaAB2TotalAndFracB(thetaC, verbose=less(verbose, 5));
  verbose && str(verbose, dataC);
 
  rm(snps, thetaC); # Not needed anymore
  verbose && exit(verbose);

  verbose && enter(verbose, "Calibrating non-polymorphic probes");
  # Extract total CNs
  units <- which(nok);
  rm(nok);
  theta <- data[units,"total",,drop=FALSE];
  dim(theta) <- dim(theta)[-2]; # Drop the 2nd dimension
  thetaC <- fitCalMaTeCNprobes(theta, references=references);
  rm(theta); # Not needed anymore

  dataC[units,"total",] <- thetaC;

##  aux <- dataC[units,,,drop=FALSE];
##  save(aux,file="dataC.Rdata")
  verbose && str(verbose, dataC[units,,,drop=FALSE]);

  rm(units, thetaC); # Not needed anymore
  verbose && exit(verbose);
  

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Standardize toward a custom average of the references?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(refAvgFcn)) {
    verbose && enter(verbose, "Standardize total copy numbers toward the average reference signals");
    # Extract reference total copy number signals
    yCR <- dataC[,"total",references,drop=FALSE];
    dim(yCR) <- dim(yCR)[-2]; # Drop 2nd dimension
    # Calculate the average
    yCR <- refAvgFcn(yCR, na.rm=TRUE);
    # Standardize total copy numbers to this average
    dataC[,"total",] <- (2/yCR) * dataC[,"total",,drop=FALSE];
    rm(yCR);
    verbose && exit(verbose);
  }

  # Enforce the same dimension names as the input data
  dimnames(dataC) <- dimnames;

  verbose && cat(verbose, "Calibrated (total,fracB) signals:");
  verbose && str(verbose, dataC);

  verbose && exit(verbose);

  dataC;
}) # calmateByTotalAndFracB()


###########################################################################
# HISTORY:
# 2011-12-15 [HB]
# o CLEANUP: Tidied up the validation of argument 'references' and
#   improved the corresponding error messages.
# 2011-12-07 [MO]
# o Number of references has to be at least 3.
# 2011-03-18 [HB]
# o BUG FIX: calmateByTotalAndFracB() required that the 2nd dimension
#   of argument 'data' had names "total" and "fracB".
# 2010-08-05 [HB]
# o ROBUSTNESS: Now calmateByTotalAndFracB() asserts that there is at 
#   least two samples.
# o BUG FIX: calmateByTotalAndFracB() assumed that there where enough
#   units and samples so that subsetting would not drop singleton 
#   dimension.  Now we use drop=FALSE everywhere.
# 2010-08-02 [HB]
# o Added argument 'refAvgFcn' to calmateByTotalAndFracB().
# o CLEANUP: Removed save() calls used for debugging.
# 2010-06-22 [MO]
# o Now calmateByTotalAndFracB() calibrates also non-polymorphic loci.
# 2010-06-18 [HB]
# o Now calmateByTotalAndFracB() handles also non-polymorphic loci.
# 2010-06-04 [MO]
# o Created.
###########################################################################
