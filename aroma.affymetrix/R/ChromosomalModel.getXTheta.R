setMethodS3("getPcuTheta", "ChromosomalModel", function(this, chromosome, reorder=TRUE, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Extracting (position, theta+)");

  # Get (position, chipType, unit) map
  pcu <- getPositionChipTypeUnit(this, chromosome=chromosome, 
                                                  verbose=less(verbose, 20));

  # Get list of chip-effect sets
  cesList <- getSets(this);

  # Allocate return structure
  naValue <- as.double(NA);
  theta <- matrix(naValue, nrow=nrow(pcu), ncol=length(this));

  for (kk in seq_along(cesList)) {
    verbose && enter(verbose, "Chip type #", kk, " of ", length(cesList));
    ces <- cesList[[kk]];
    idxs <- which(as.integer(pcu[,"chipType"]) == kk);  
    units <- pcu[idxs,"unit"];
    verbose && cat(verbose, "Units: ");
    verbose && str(verbose, units);

    verbose && enter(verbose, "Reading data across chip-effect files");
    theta[idxs,] <- extractMatrix(ces, units=units, verbose=less(verbose, 20));
    verbose && exit(verbose);

    verbose && exit(verbose);
  }

  verbose && exit(verbose);

  list(pcu=pcu, theta=theta);
}, protected=TRUE);  # getPcuTheta()


##############################################################################
# HISTORY:
# 2008-07-20
# o Updated the following methods to preallocate matrixes with the correct
#   data type to avoid coercing later: getPcuTheta().
# 2007-09-25
# o Moved getPcuTheta() to ChromosomalModel.
# 2007-09-24
# o Added getXTheta().
##############################################################################
