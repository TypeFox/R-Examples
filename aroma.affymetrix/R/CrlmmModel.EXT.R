setMethodS3("getPlatformDesignDB", "CrlmmModel", function(this, ..., verbose=FALSE) {
  requireNamespace("oligoClasses") || throw("Package not loaded: oligoClasses")
  db <- oligoClasses::db


  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting Platform Design Database");
  ds <- getDataSet(this);
  cdf <- getCdf(ds);
  chipType <- getChipType(cdf, fullname=FALSE);
  verbose && cat(verbose, "Chip type: ", chipType);
  pdPkgName <- .cleanPlatformName(chipType);
  verbose && cat(verbose, "Plaform Design package: ", pdPkgName);

  require(pdPkgName, character.only=TRUE) || throw("Package not loaded: ", pdPkgName);

  pdDB <- db(get(pdPkgName, mode="S4"));
  verbose && print(verbose, pdDB);
  verbose && exit(verbose);
  pdDB;
}, private=TRUE)



setMethodS3("getCrlmmPriors", "CrlmmModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting CRLMM priors");
  ds <- getDataSet(this);
  cdf <- getCdf(ds);
  chipType <- getChipType(cdf, fullname=FALSE);
  verbose && cat(verbose, "Chip type: ", chipType);

  res <- NULL;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Trying oligoParams package
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pkgName <- "oligoParams";
  if (isPackageInstalled(pkgName) && require(pkgName, character.only=TRUE)) {
    # To trick R CMD check
    getCrlmmSnpNames <- NULL; rm(list="getCrlmmSnpNames");
    verbose && enter(verbose, "Querying oligoParams");
    tryCatch({
      res <- getCrlmmSnpNames(chipType, tags="SNPs",
                                           verbose=less(verbose, -20));
    }, error=function(ex) {})
    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Trying PD package
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(res)) {
    verbose && enter(verbose, "Querying PD package");

    pdPkgName <- .cleanPlatformName(chipType);
    verbose && cat(verbose, "Platform Design (PD) package: ", pdPkgName);

    # Load target from PD package
    path <- system.file(package=pdPkgName);
    if (path == "") {
      throw("Cannot load HapMap reference target quantiles. Package not installed: ", pdPkgName);
    }

    verbose && enter(verbose, "Loading CRLMM priors etc");
    path <- file.path(path, "extdata");
    path <- Arguments$getReadablePath(path);
    filename <- sprintf("%sCrlmmInfo.rda", pdPkgName);
    pathname <- Arguments$getReadablePathname(filename, path=path);
    verbose && cat(verbose, "Pathname: ", pathname);
    key <- sprintf("%sCrlmm", pdPkgName);
    res <- loadToEnv(pathname)[[key]];
    verbose && exit(verbose);

    verbose && exit(verbose);
  }

  verbose && cat(verbose, "Loaded data:");
  verbose && capture(verbose, ll(envir=res));

  verbose && exit(verbose);

  res;
}, protected=TRUE);



setMethodS3("getCrlmmSNPs", "CrlmmModel", function(this, flavor=c("oligoPD", "oligoCDF"), ..., verbose=FALSE) {
  requireNamespace("DBI") || throw("Package not loaded: DBI")
  dbGetQuery <- DBI::dbGetQuery


  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Identifying SNP according to oligo::CRLMM");
  verbose && cat(verbose, "Flavor: ", flavor);


  if (flavor == "oligoCDF") {
    verbose && enter(verbose, "Querying the CDF");
    # Immitates oligo
    ds <- getDataSet(this);
    cdf <- getCdf(ds);
    # Identify all SNP* units.
    units <- indexOf(cdf, pattern="^SNP");
    verbose && exit(verbose);
  } else if (flavor == "oligoPD") {
    verbose && enter(verbose, "Querying");
    res <- NULL;

    pkgName <- "oligoParams";
    if (isPackageInstalled(pkgName) && require(pkgName, character.only=TRUE)) {
      # To trick R CMD check
      getCrlmmSnpNames <- NULL; rm(list="getCrlmmSnpNames");
      ds <- getDataSet(this);
      cdf <- getCdf(ds);
      chipType <- getChipType(cdf, fullname=FALSE);

      verbose && enter(verbose, "Querying oligoParams");
      tryCatch({
        res <- getCrlmmSnpNames(chipType, tags="SNPs",
                                             verbose=less(verbose, -20));
      }, error=function(ex) {})
      verbose && exit(verbose);
    }

    if (is.null(res)) {
      verbose && enter(verbose, "Querying the PD package");

      pdDB <- getPlatformDesignDB(this, verbose=less(verbose,1));
      verbose && print(verbose, pdDB);

      res <- dbGetQuery(pdDB, "SELECT man_fsetid FROM featureSet WHERE man_fsetid LIKE 'SNP%' ORDER BY man_fsetid")[[1]];
      verbose && str(verbose, res);
      verbose && exit(verbose);
    }

    verbose && enter(verbose, "Mapping to CDF unit indices");
    ds <- getDataSet(this);
    cdf <- getCdf(this);
    units <- indexOf(cdf, names=res);
    names(units) <- res;
    verbose && exit(verbose);

    verbose && exit(verbose);
  }

  verbose && str(verbose, units);

  verbose && exit(verbose);

  units;
}, private=TRUE)  # getCrlmmSNPs()


setMethodS3("getCrlmmSNPsOnChrX", "CrlmmModel", function(this, flavor=c("oligoPD", "oligoCDF"), ..., verbose=FALSE) {
  requireNamespace("DBI") || throw("Package not loaded: DBI")
  dbGetQuery <- DBI::dbGetQuery


  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Identifying all SNPs on ChrX");

  if (flavor == "oligoCDF") {
    verbose && enter(verbose, "Querying the CDF and UGP");

    # Immitates oligo
    ds <- getDataSet(this);
    cdf <- getCdf(ds);

    # Identify units on ChrX.
    gi <- getGenomeInformation(cdf);
    units <- getUnitsOnChromosome(gi, 23);

    # Identify which are SNP* units.
    unitNames <- getUnitNames(cdf, units=units);
    units <- units[grep("^SNP", unitNames)];
    # Not needed anymore
    unitNames <- NULL;

    verbose && exit(verbose);
  } else if (flavor == "oligoPD") {
    verbose && enter(verbose, "Querying");
    res <- NULL;

    pkgName <- "oligoParams";
    if (isPackageInstalled(pkgName) && require(pkgName, character.only=TRUE)) {
      # To trick R CMD check
      getCrlmmSnpNames <- NULL; rm(list="getCrlmmSnpNames");
      ds <- getDataSet(this);
      cdf <- getCdf(ds);
      chipType <- getChipType(cdf, fullname=FALSE);

      verbose && enter(verbose, "Querying oligoParams");
      tryCatch({
        res <- getCrlmmSnpNames(chipType, tags="SNPs,ChrX",
                                             verbose=less(verbose, -20));
      }, error=function(ex) {})
      verbose && exit(verbose);
    }

    if (is.null(res)) {
      verbose && enter(verbose, "Querying the PD package");
      pdDB <- getPlatformDesignDB(this, verbose=less(verbose,1));
      verbose && print(verbose, pdDB);

      res <- dbGetQuery(pdDB, "SELECT man_fsetid FROM featureSet WHERE man_fsetid LIKE 'SNP%' AND chrom = 'X'")[[1]];
      verbose && str(verbose, res);
    }

    verbose && enter(verbose, "Mapping to CDF unit indices");
    ds <- getDataSet(this);
    cdf <- getCdf(this);
    units <- indexOf(cdf, names=res);
    names(units) <- res;
    verbose && exit(verbose);

    verbose && exit(verbose);
  }

  verbose && str(verbose, units);

  verbose && exit(verbose);

  units;
}, private=TRUE) # getCrlmmSNPsOnChrX()



setMethodS3("getCrlmmSplineParameters", "CrlmmModel", function(this, flavor=c("oligoPD"), ..., verbose=FALSE) {
  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Retrieving spline parameters");
  verbose && cat(verbose, "Flavor: ", flavor);

  ds <- getDataSet(this);
  cdf <- getCdf(ds);
  chipType <- getChipType(cdf, fullname=FALSE);
  verbose && cat(verbose, "Chip type: ", chipType);

  if (flavor == "oligoPD") {
    verbose && enter(verbose, "Querying the PD package");

    pdPkgName <- .cleanPlatformName(chipType);
    verbose && cat(verbose, "Platform Design (PD) package: ", pdPkgName);

    # Load target from PD package
    path <- system.file(package=pdPkgName);
    if (path == "") {
      throw("Cannot load spline parameters. Package not installed: ", pdPkgName);
    }

    path <- file.path(path, "extdata");
    path <- Arguments$getReadablePath(path);
    filename <- sprintf("%s.spline.params.rda", pdPkgName);
    pathname <- Arguments$getReadablePathname(filename, path=path);
    verbose && cat(verbose, "Pathname: ", pathname);
    res <- loadToEnv(pathname);
    verbose && exit(verbose);
  }

  verbose && cat(verbose, "Loaded data:");
  verbose && print(verbose, ll(envir=res));

  verbose && str(verbose, res);

  verbose && exit(verbose);

  res;
}, private=TRUE)


############################################################################
# HISTORY:
# 2009-01-12
# o Added getCrlmmSplineParameters().
# o Warnings about missing (optional) oligoParams is no longer reported.
# 2009-01-10
# o Now data in oligoParams is used, if installed.
# o Updated to work with latest aroma.core and aroma.affymetrix.
# 2009-01-07
# o Overriding getAsteriskTags() and getRootPath().
# 2008-12-08
# o Now setup is much more like fit() for ProbeLevelModel.
# 2008-12-07
# o Starting to make use of AromaCrlmmBinarySet.
# o Created CrlmmModel from justCRLMMv2().
# 2008-12-05
# o Created from justCRLMMv2() of oligo v1.7.3.
############################################################################
