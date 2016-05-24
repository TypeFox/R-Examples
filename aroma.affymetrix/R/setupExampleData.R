###########################################################################/**
# @set class=AromaAffymetrix
# @RdocMethod setupExampleData
#
# @title "Setups example data in the current directory"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{dataset, chipType}{@character strings specifying the data set
#    and the chip type to install.}
#   \item{dirs}{A @character @vector specifying which directories to setup.}
#   \item{mustWork}{If @TRUE, an error is thrown if the requested data set
#    could not be installed, otherwise not.}
#   \item{validate}{If @TRUE, the installed files are also validated,
#    otherwise not.}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns (invisibly) @TRUE if all requested data was installed,
#   otherwise @FALSE.
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("setupExampleData", "AromaAffymetrix", function(pkg, dataset=NULL, chipType=NULL, dirs=c("annotationData", "rawData"), mustWork=TRUE, validate=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataset':
  # BACKWARD compatibility
  if (is.null(dataset)) dataset <- "FusionSDK_HG-Focus"
  dataset <- Arguments$getCharacter(dataset)

  # Argument 'chipType':
  if (!is.null(chipType)) {
    chipType <- Arguments$getCharacter(chipType)
  }

  # Argument 'dirs':
  dirs <- Arguments$getCharacters(dirs)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify data set
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # BACKWARD compatibility
  if (dataset == "HG-Focus") dataset <- "FusionSDK_HG-Focus"

  ## Check 'data set'
  knownDataSets <- c("FusionSDK_HG-Focus", "FusionSDK_Test3")
  if (!is.element(dataset, knownDataSets)) {
    throw(sprintf("Data set '%s' is not among the set of known data sets: ", dataset, paste(sQuote(knownDataSets), collapse=", ")))
  }

  ## Default chip type?
  if (is.null(chipType)) {
    if (dataset == "FusionSDK_HG-Focus") {
      chipType <- "HG-Focus"
    } else if (dataset == "FusionSDK_Test3") {
      chipType <- "Test3"
    } else {
      throw("There is no known default chip type for this data set: ", dataset)
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check needed packages
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (pkg in c("AffymetrixDataTestFiles", "affxparser")) {
    if (!isPackageInstalled(pkg)) {
      if (!mustWork) return(invisible(FALSE))
      throw("Package not installed: ", pkg)
    }
  }

  path0 <- system.file(package="AffymetrixDataTestFiles", mustWork=TRUE)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # annotationData/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pathR <- "annotationData"
  if (is.element(pathR, dirs)) {
    pathD <- file.path(pathR, "chipTypes", chipType)
    filename <- sprintf("%s.CDF", chipType)
    pathnameD <- Arguments$getReadablePathname(filename, path=pathD, mustExist=FALSE)

    # Missing?
    if (!isFile(pathnameD)) {
      pathnameD <- Arguments$getWritablePathname(pathnameD)
      pathS <- file.path(path0, pathD)
      pathS <- Arguments$getReadablePath(pathS)

      ## Locate the/a CDF
      pathnameS <- findFiles(pattern=sprintf("^%s$", filename), paths=pathS, recursive=TRUE, firstOnly=TRUE)
      if (is.null(pathnameS)) {
        throw(sprintf("No such file under '%s': '%s'", pathS, filename))
      }

      ## Setup CDF
      cdfS <- AffymetrixCdfFile(pathnameS)

      # Convert to binary XDA/v4?
      formatS <- getFileFormat(cdfS)
      if (regexpr("v4", formatS) != -1L) {
        # Link
        createLink(link=pathnameD, target=pathnameS)
      } else {
        # Convert source ASCII CDF to binary CDF
        .convertCdf(pathnameS, pathnameD, .validate=FALSE, verbose=TRUE)
      }

      # Sanity check
      if (validate) {
        cdfD <- AffymetrixCdfFile(pathnameD)
        formatD <- getFileFormat(cdfD)
        stopifnot(regexpr("v4", formatD) != -1L)
      }
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # rawData/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pathR <- "rawData"
  if (is.element(pathR, dirs)) {
    pathD <- file.path(pathR, dataset, chipType)
    pathD <- Arguments$getWritablePath(pathD)
    pathS <- file.path(path0, pathD)
    pathS <- Arguments$getReadablePath(pathS)

    ## Locate the CEL files (if more than one subdirectory/file format type
    ## then use the first one found)
    pathnameS <- findFiles(pattern="[.]CEL$", paths=pathS, recursive=TRUE, firstOnly=TRUE)
    pathS <- dirname(pathnameS)

    ## Link to all source files
    pathnamesS <- list.files(path=pathS, pattern="[.]CEL$", full.names=TRUE)
    for (ii in seq_along(pathnamesS)) {
      pathnameS <- pathnamesS[ii]
      pathnameD <- file.path(pathD, basename(pathnameS))
      createLink(link=pathnameD, target=pathnameS, skip=TRUE)
    }

    # Sanity check
    if (validate && isDirectory("annotationData")) {
      cels <- AffymetrixCelSet$byName(dataset, chipType=chipType)
    }
  }

  invisible(TRUE)
}, protected=TRUE) # setupExampleData()


############################################################################
# HISTORY:
# 2015-01-23
# o Added support for dataset='FusionSDK_Test3' in addition to the already
#   supported dataset='FusionSDK_HG-Focus'.
# o Added arguments 'dataset' and 'chipType'.
# 2014-06-24
# o Created from aroma.seq ditto.
############################################################################
