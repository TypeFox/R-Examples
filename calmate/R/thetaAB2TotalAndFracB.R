###########################################################################/**
# @set "class=array"
# @RdocMethod thetaAB2TotalAndFracB
# @alias thetaAB2TotalAndFracB
# @alias totalAndFracB2ThetaAB
# @alias totalAndFracB2ThetaAB.array
# 
# @title "Converts an Jx2xI array between (thetaA,thetaB) and (total,fracB) formats"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{data}{An Jx2xI @numeric array, where J is the number of SNPs,
#          2 is the number of alleles, and I is the number of samples.}
#  \item{...}{Not used.}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# @examples "../incl/thetaAB2TotalAndFracB.Rex"
#
# \value{
#   Returns an Jx2xI @numeric array.
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("thetaAB2TotalAndFracB", "array", function(data, ..., verbose=FALSE) {
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
    if (!identical(dimnames[[2]], c("A", "B"))) {
      throw("If given, the names of the allele (2nd) dimension of the Jx2xI-dimensional array (argument 'data') have to be 'A' & 'B': ", paste(dimnames[[2]], collapse=", "));
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "thetaAB2TotalAndFracB()");
  verbose && cat(verbose, "ASCN signals:");
  verbose && str(verbose, data);

  dataT <- data;
  dimnames(dataT)[[2]] <- c("total", "fracB");
  dataT[,"total",] <- data[,"A",,drop=FALSE] + data[,"B",,drop=FALSE];
  dataT[,"fracB",] <- data[,"B",,drop=FALSE] / dataT[,"total",,drop=FALSE];
  verbose && str(verbose, dataT);

  verbose && exit(verbose);

  dataT;
}) # thetaAB2TotalAndFracB()



setMethodS3("totalAndFracB2ThetaAB", "array", function(data, ..., verbose=FALSE) {
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
    if (!identical(dimnames[[2]], c("total", "fracB"))) {
      throw("If given, the names of the allele (2nd) dimension of the Jx2xI-dimensional array (argument 'data') have to be 'total' & 'fracB': ", paste(dimnames[[2]], collapse=", "));
    }
  }
  

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "totalAndFracB2ThetaAB()");
  verbose && cat(verbose, "(total,fracB) signals:");
  verbose && str(verbose, data);

  dataT <- data;
  dimnames(dataT)[[2]] <- c("A", "B");
  dataT[,"B",] <- data[,"total",,drop=FALSE] * data[,"fracB",,drop=FALSE];
  dataT[,"A",] <- data[,"total",,drop=FALSE] - dataT[,"B",,drop=FALSE];
  verbose && str(verbose, dataT);

  verbose && exit(verbose);

  dataT;
}, protected=TRUE) # totalAndFracB2ThetaAB()


###########################################################################
# HISTORY:
# 2010-05-19
# o Created.
###########################################################################
