setMethodS3("extractChromosomalDataFrame", "ChipEffectFile", function(this, units=NULL, ..., chromosomes=NULL, orderBy=c("chromosome", "physicalPosition"), decreasing=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'chromosomes':
  if (!is.null(chromosomes)) {
    chromosomes <- Arguments$getIndices(chromosomes, max=999);
    chromosomes <- sort(unique(chromosomes));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read annotation data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading (chromosome, position):");
  cdf <- getCdf(this);
  gi <- getGenomeInformation(cdf, verbose=less(verbose, 10));
  gp <- readDataFrame(gi, units=units, verbose=less(verbose, 20));
  verbose && str(verbose, gp);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Subset by chromosomes?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(chromosomes)) {
    verbose && enter(verbose, "Subsetting by chromosome(s)");
    verbose && cat(verbose, "Chromosomes:");
    verbose && print(verbose, chromosomes);
    keep <- which(gp[,1] %in% chromosomes);
    gp <- gp[keep,,drop=FALSE];
    if (is.null(units)) {
      units <- keep;
    } else {
      units <- units[keep];
    }
    verbose && cat(verbose, "Selected units:");
    verbose && str(verbose, units);
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Order by?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(orderBy)) {
    verbose && enter(verbose, "Ordering data");
    verbose && cat(verbose, "Order by:");
    verbose && print(verbose, orderBy);
    args <- lapply(orderBy, FUN=function(cc) gp[,cc,drop=FALSE]);
    args <- c(args, decreasing=decreasing);
    o <- do.call(order, args)
    gp <- gp[o,,drop=FALSE];
    if (is.null(units)) {
      units <- o;
    } else {
      units <- units[o];
    }
    # Not needed anymore
    o <- NULL;
    verbose && cat(verbose, "Ordered units:");
    verbose && str(verbose, units);
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading data");
  verbose && cat(verbose, "Units:");
  verbose && str(verbose, units);
  data <- extractDataFrame(this, units=units, ..., verbose=less(verbose, 1));
  verbose && str(verbose, data);
  verbose && exit(verbose);

  # Expand (chromosome, position) map?
  if (!identical(units, data$unit)) {
    verbose && enter(verbose, "Expanding (chromosome, position) map");
    idxs <- match(data$unit, units);
    verbose && cat(verbose, "Expansion map:");
    verbose && str(verbose, idxs);
    gp <- gp[idxs,,drop=FALSE];
    verbose && str(verbose, gp);
    verbose && exit(verbose);
  }
  # Not needed anymore
  units <- NULL;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Insert (chromosome, position)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify where to insert
  afterColumns <- c("unitName", "groupName", "unit", "group", "cell");
  at <- match(afterColumns, colnames(data));
  at <- max(at, na.rm=TRUE);

  dataHead <- data[,1:at,drop=FALSE];
  dataTail <- data[,(at+1):ncol(data),drop=FALSE];
  # Not needed anymore
  data <- at <- NULL;

  data <- cbind(dataHead, gp, dataTail);

  # Not needed anymore
  dataHead <- dataTail <- gp <- NULL;

  data;
}, protected=TRUE)


# The exact same code can be used for both ChipEffectFile and ChipEffectSet.
setMethodS3("extractChromosomalDataFrame", "ChipEffectSet", extractChromosomalDataFrame.ChipEffectFile, protected=TRUE)



############################################################################
# HISTORY:
# 2008-12-29
# o Created.
############################################################################
