###########################################################################/**
# @RdocClass AbstractPSCNData
#
# @title "The AbstractPSCNData class"
#
# \description{
#  @classhierarchy
#
#  A AbstractPSCNData object holds parent-specific copy number data.
# }
# 
# @synopsis
#
# \arguments{
#   \item{chromosome}{(Optional) An @integer scalar (or a @vector of length J),
#        which can be used to specify which chromosome each locus belongs to
#        in case multiple chromosomes are segments.
#        This argument is also used for annotation purposes.}
#   \item{x}{Optional @numeric @vector of J genomic locations.
#            If @NULL, index locations \code{1:J} are used.}
#   \item{isSNP}{An optional @logical @vector of length J specifying
#        whether each locus is a SNP or not (non-polymorphic loci).}
#   \item{mu}{An optional @numeric @vector of J genotype calls in 
#        \{0,1/2,1\} for AA, AB, and BB, respectively, 
#        and @NA for non-polymorphic loci.}
#   \item{...}{Optional named locus-specific signal @vectors of length J.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("AbstractPSCNData", function(chromosome=NULL, x=NULL, isSNP=NULL, mu=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(chromosome)) {
    # Is first argument special?
    if (is.data.frame(chromosome)) {
      data <- chromosome;

      chromosome <- data$chromosome;
      x <- data$x;
      isSNP <- data$isSNP;
      mu <- data$mu;
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Required arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Argument 'chromosome':
    disallow <- c("Inf");
    chromosome <- Arguments$getIntegers(chromosome, range=c(0,Inf), disallow=disallow);
    nbrOfLoci <- length(chromosome);
    length2 <- rep(nbrOfLoci, times=2);
  

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Optional arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Argument 'isSNP':
    if (!is.null(isSNP)) {
      isSNP <- Arguments$getLogicals(isSNP, length=length2, disallow="NA");
    }
  
    # Argument 'mu':
    if (!is.null(mu)) {
      mu <- Arguments$getDoubles(mu, range=c(0,1), length=length2, disallow="Inf");
    }
  }

##  cat("AbstractPSCNData...\n");

  this <- extend(AbstractCNData(chromosome=chromosome, x=x, isSNP=isSNP, mu=mu, ...), "AbstractPSCNData");

##print(head(as.data.frame(this)));
##print(class(this));
##  cat("AbstractPSCNData...done\n");

  this;
})



setMethodS3("getSNPFields", "AbstractPSCNData", function(this, ...) {
  data <- as.data.frame(this);
  fields <- colnames(data);
  fields <- grep("^(beta|mu)", fields, value=TRUE);
  fields;
}, protected=TRUE)



setMethodS3("callSNPs", "AbstractPSCNData", function(this, method=c("some", "all"), force=FALSE, ..., verbose=FALSE) {
  data <- this;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'method':
  method <- match.arg(method);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

 
  isSNP <- data$isSNP;
  if (!force && !is.null(isSNP)) {
    verbose && cat(verbose, "Already done. Skipping.");
    verbose && exit(verbose);
    return(data);
  }

  verbose && enter(verbose, "Calling SNPs based on BAFs and genotype calls");
  snpFields <- getSNPFields(this);
  verbose && cat(verbose, "SNP fields: ", hpaste(snpFields));

  if (length(snpFields) == 0) {
    throw("Cannot identify SNPs, because there are no SNP fields.");
  }

  # First, assume none of the loci are SNPs
  nbrOfLoci <- nrow(data);
  if (method == "some") {
    isSNP <- rep(FALSE, times=nbrOfLoci);
    fcn <- function(x,y) { x | y };
  } else if (method == "all") {
    isSNP <- rep(TRUE, times=nbrOfLoci);
    fcn <- function(x,y) { x & y };
  }

  # Then remove those for which we would expect to have SNP signals
  for (field in snpFields) {
    signals <- data[[field]];
    isSNP <- fcn(isSNP, !is.na(signals));
  } # for (ff ...)

  # Sanity check
  stopifnot(all(is.finite(isSNP)));

  nbrOfSNPs <- sum(isSNP);
  verbose && printf(verbose, "Number of SNPs: %d (%.4g%% of %d loci)\n", nbrOfSNPs, 100*nbrOfSNPs/nbrOfLoci, nbrOfLoci); 

  data$isSNP <- isSNP;
  verbose && exit(verbose);

  data;
}, protected=TRUE)



############################################################################
# HISTORY:
# 2012-03-10
# o Renames isSNP() to callSNPs().
# 2012-02-29
# o Added isSNP() for AbstractPSCNData.
# o Created.
############################################################################
