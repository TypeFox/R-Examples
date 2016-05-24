setMethodS3("getImage", "AromaUnitTabularBinaryFile", function(this, transforms=NULL, xrange=c(0,Inf), yrange=xrange, zrange=c(0,sqrt(2^16)), field, levels=NULL, zoom=1, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'field':
#  field <- match.arg(field);
  
  # Argument 'zoom':
  zoom <- Arguments$getDouble(zoom, range=c(0,Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Read data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Getting ", class(this)[1], " image");

  verbose && enter(verbose, "Reading (unit,group,cell) map");
  cdf <- getCdf(this);
  ugcMap <- getUnitGroupCellMap(cdf, verbose=less(verbose, 5));
  verbose && exit(verbose);

  verbose && enter(verbose, "Reading data");
  values <- readDataFrame(this, fields=field, ..., verbose=less(verbose,1));
  verbose && str(verbose, values);
  values <- values[[field]];
#  verbose && str(verbose, values);
  z <- vector(mode=mode(values), 1);
  z <- matrix(z, nrow=nbrOfRows(cdf), ncol=nbrOfColumns(cdf));
  z[indexByRow(z, ugcMap[,"cell"])] <- values[ugcMap[,"unit"]];
  # Not needed anymore
  values <- ugcMap <- NULL;
  verbose && summary(verbose, as.vector(z));
  verbose && printf(verbose, "RAM: %.1fMB\n", object.size(z)/1024^2);
  verbose && exit(verbose);

  verbose && enter(verbose, "Transforming data");
  dim <- dim(z);
  mode <- mode(z);
  if (mode == "character") {
    if (is.null(levels)) {
      z <- factor(z);
    } else {
      z <- factor(z, levels=levels);
    }
    z <- as.integer(z);
  }
  dim(z) <- dim;
  verbose && str(verbose, z);
  verbose && exit(verbose);

  verbose && enter(verbose, "Creating Image");
  img <- getImage(z, scale=zoom, lim=zrange, ..., verbose=less(verbose, 1));
  verbose && print(verbose, img);
  verbose && exit(verbose);

  verbose && exit(verbose);

  # Return the 'field'
  attr(img, "field") <- field;

  img;
})


############################################################################
# HISTORY:
# 2011-01-30
# o CLEAN UP: getImage() for AromaUnitTabularBinaryFile no longer explicitly
#   require the EBImage but instead calls aroma.core's getImage(). 
#   The latter method depends on EBimage for now.
# 2008-04-14
# o Created from AffymetrixCdfFile.PLOT.R.
############################################################################
