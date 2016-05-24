setConstructorS3("AromaAffymetrix", function(...) {
  extend(AromaPackage("aroma.affymetrix"), "AromaAffymetrix");
})

setMethodS3("fixSearchPath", "AromaAffymetrix", function(this, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # RULES
  # 2008-02-27:
  # o affy must be after R.huge, otherwise the former overrides the
  #   generic colnames() function of the latter.
  # o affyPLM must be after aroma.affymetrix.
  # o EBImage must be after aroma.affymetrix.
  # 2008-08-27:
  # o affy must be after aroma.light, otherwise the former overrides
  #   the generic plotDensity() function of the latter.
  # 2009-01-10:
  # o oligo must be after aroma.affymetrix, otherwise the former
  #   overrides generic justSNPRMA().
  # 2011-06-07:
  # o ggplot2 must be after aroma.affymetrix, otherwise the former
  #   overrides generic rescale().
  # 2013-08-03:
  # o affxparser should be after aroma.affymetrix so that the
  #   generic writeCdf() of aroma.affymetrix is found.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Figure out which of our packages (aroma.core, aroma.light etc.) is
  # last on the search path.
  aheadPkgs <- c("aroma.affymetrix", "aroma.light", "R.huge", "R.oo");

  # Problematic package that must be after the packages on the search path
  behindPkgs <- c("affxparser", "affy", "affyPLM", "EBImage", "oligo", "ggplot2");

  res <- fixSearchPathInternal(this, aheadPkgs=aheadPkgs,
                                           behindPkgs=behindPkgs, ...);

  # Make rescale() for numeric:s compatible for ggplot2 users.
  if (isPackageLoaded("ggplot2") && !exists("rescale.numeric", mode="function")) {
    # Get ggplot2::rescale() without triggering a R CMD check warning.
    ggplot2 <- Package("ggplot2");
    fcn <- get("rescale", mode="function", envir=getEnvironment(ggplot2));
    # Add the rescale() for numeric
    assign("rescale.numeric", fcn, envir=getEnvironment(this));
  }

  # Return the package actually moved
  invisible(res);
})


setMethodS3("update", "AromaAffymetrix", function(object, ..., verbose=FALSE) {
  # To please R CMD check
  this <- object;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Checking for and install updates");
  verbose && cat(verbose, "Package: ", getName(this));
  verbose && printf(verbose, "Current version: v%s (%s)\n",
                                        getVersion(this), getDate(this));


  state <- 0;

  url <- "http://www.braju.com/R/hbLite.R";
  verbose && enter(verbose, "Trying to download update script");
  verbose && cat(verbose, "URL: ", url);
  hbInstall <- NULL;
  tryCatch({
    suppressWarnings({
      source(url);
    })
    state <- 1;
    verbose && exit(verbose);
  }, error = function(ex) {
    verbose && exit(verbose, suffix="failed");
    throw(ex);
  })

  verbose && enter(verbose, "Launching update command");
  verbose && printf(verbose, "Call: hbInstall(\"%s\")", getName(this));
  tryCatch({
    hbInstall(getName(this));
    state <- 2;
  }, error = function(ex) {
    verbose && exit(verbose, suffix="failed");
    throw(ex);
  })

  verbose && cat(verbose, "Package has been updated.");

  verbose && exit(verbose);

  invisible(state);
})



setMethodS3("getDefaultSettings", "AromaAffymetrix", function(this, ...) {
  defaults <- list(
    memory = list(
      ram = 1,
      gcArrayFrequency = 50
    ),

    rules = list(
      allowAsciiCdfs = FALSE,
      allowDChipAnnotationFiles = FALSE
    ),

    output = list(
      # Should the checksum be reported when print():ing files?
      checksum = FALSE,

      # Max number of arrays for which to report timestamps
      timestampsThreshold = 500
    ),

    models = list(
      RmaPlm = list(
       # Number of cells *and* arrays for using median polish
        medianPolishThreshold  = c( 500, 6),
       # Number of cells *and* arrays for skipping unit group
        skipThreshold          = c(5000, 1)
      )
    )
  );

  defaults;
}, protected=TRUE)



############################################################################
# HISTORY:
# 2013-08-03
# o Now affxparser is moved behind aroma.affymetrix on the search path.
# 2011-06-07
# o Now the search path is adjusted such that 'ggplot2' comes after
#   'aroma.affymetrix', because the former overrides the generic
#   rescale() function of the latter with a non-generic function.
# 2010-06-07
# o Added setting 'rules/allowDChipAnnotationFiles'. Still to be asserted
#   by the code.
# 2009-05-16
# o Updated fitSearchPath() to utilize new fitSearchPathInternal().
# o Now AromaAffymetrix inherits from AromaPackage.
# 2009-01-10
# o Now the oligo package is forced to be after aroma.affymetrix.
# 2008-08-27
# o Now the affy, affyPLM, and EBImage packages are forced to be after all
#   of aroma.affymetrix, aroma.light, and R.huge on the search() path.
# 2007-12-13
# o Added update() and patch() to the AromaAffymetrix Package class.
# 2007-03-06
# o Created.
############################################################################
