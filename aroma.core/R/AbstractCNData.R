###########################################################################/**
# @RdocClass AbstractCNData
#
# @title "The AbstractCNData class"
#
# \description{
#  @classhierarchy
#
#  An AbstractCNData object holds copy number data.
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
#   \item{y}{Optional @numeric @vector of J genomic locations.}
#   \item{...}{Optional named locus-specific signal @vectors of length J.}
#   \item{name}{Optional @character string.}
#   \item{.virtuals}{(Internal) a @list with virtual column name @functions.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("AbstractCNData", function(chromosome=NULL, x=NULL, y=NULL, ..., name=NULL, .virtuals=NULL) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(chromosome)) {
    # Argument 'chromosome':
    if (is.data.frame(chromosome)) {
      data <- chromosome;
      chromosome <- data$chromosome;
      x <- data$x;
      y <- data$y;
    }

    # Any extra locus signals?
    args <- list(...);
    keep <- !sapply(args, FUN=is.null);
    args <- args[keep];

    if (is.null(y)) {
      y <- args$y;
      args$y <- NULL;
    }
  } else {
    args <- NULL;
  }

# str(list(chromosome=chromosome, x=x, cn=cn));


##  cat("AbstractCNData...\n");

  this <- extend(RawGenomicSignals(chromosome=chromosome, x=x, y=y, name=name), "AbstractCNData");

  this <- setVirtualColumnFunctions(this, .virtuals);

## print(class(this));

  if (length(args) > 0) {
    class <- class(this);
    for (key in names(args)) {
      values <- args[[key]];
      this[[key]] <- values;
    }
    class(this) <- class;
  }
##print(class(this));


##print(head(as.data.frame(this)));
##  cat("AbstractCNData...done\n");

  this;
})


setMethodS3("getPlatform", "AbstractCNData", function(this, ...) {
  attr(this, "platform");
})



setMethodS3("setPlatform", "AbstractCNData", function(this, platform, ...) {
  platform <- Arguments$getCharacter(platform);
  attr(this, "platform") <- platform;
  invisible(this);
})


setMethodS3("getChipType", "AbstractCNData", function(this, fullname=TRUE, ...) {
  chipType <- attr(this, "chipType");
  if (!fullname) {
    chipType <- gsub(",.*", "", chipType);
  }
  chipType;
})

setMethodS3("setChipType", "AbstractCNData", function(this, chipType, ...) {
  chipType <- Arguments$getCharacter(chipType);
  attr(this, "chipType") <- chipType;
  invisible(this);
})


setMethodS3("getLocusData", "AbstractCNData", function(this, ...) {
  as.data.frame(this);
})



setMethodS3("hasKnownPositions", "AbstractCNData", function(this, ...) {
  data <- as.data.frame(this, ...);
  ok <- (!is.na(data$chromosome) & !is.na(data$x));
  ok;
}, protected=TRUE)



setMethodS3("orderAlongGenome", "AbstractCNData", function(this, ...) {
  sort(this, ...);
}, protected=TRUE)



setMethodS3("findLargeGaps", "AbstractCNData", function(chromosome, ...) {
  pkg <- "PSCBS";
  require(pkg, character.only=TRUE) || throw("Package not loaded: ", pkg);

  # To please R CMD check
  this <- chromosome;

  data <- as.data.frame(this, ...);
  findLargeGaps(chromosome=data$chromosome, x=data$x, ...);
}) # findLargeGaps()



############################################################################
# HISTORY:
# 2012-03-09
# o Started to inherit from RawGenomicSignals of aroma.core v2.4.10,
#   which is now a data.frame.
# 2012-02-29
# o Added getLocusData() for AbstractCNData.
# o Added findLargeGaps() for AbstractCNData.
# o Created.
############################################################################
