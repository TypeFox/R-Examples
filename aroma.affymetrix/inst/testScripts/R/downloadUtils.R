library("aroma.affymetrix");

# Local functions
getDataSet <- function(path, ...) {
  chipType <- basename(path);
  cdf <- tryCatch({
    AffymetrixCdfFile$byChipType(chipType);
  }, error=function(ex) {
    AffymetrixCdfFile$byChipType(chipType, tags=".*");
  });
  ds <- AffymetrixCelSet$byPath(path, cdf=cdf);
  ds;
}

getRawDataSetPath <- function(dataSet=NULL, tags=NULL, chipType=NULL, ...) {
  if (is.null(dataSet) && is.null(chipType)) return(NULL);
  dataSetF <- paste(c(dataSet, tags), collapse=",");
  path <- file.path("rawData", dataSetF, chipType);
  path <- Arguments$getWritablePath(path);
  path;
} # getRawDataSetPath()



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# HapMap
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
getHapMapUrlPath <- function(chipType=c("Mapping50K_Hind240", "Mapping50K_Xba240", "GenomeWideSNP_6"), ...) {
  # Argument 'chipType':
  chipType <- match.arg(chipType);

  dirs <- c("Mapping50K_Hind240"="affy100k", "Mapping50K_Xba240"="affy100k",
            "GenomeWideSNP_6"="hapmap3_affy6.0");
  dir <- dirs[chipType];

  sprintf("http://hapmap.ncbi.nlm.nih.gov/downloads/raw_data/%s", dir);
} # getHapMapUrlPath()
 

# url <- getHapMapUrl("CEU_NA06985", chipType="Mapping50K_Hind240");
getHapMapUrl <- function(sampleName, chipType=c("Mapping50K_Hind240", "Mapping50K_Xba240"), ...) {
  # Argument 'sampleName':
  sampleName <- Arguments$getCharacters(sampleName);

  # Argument 'chipType':
  chipType <- match.arg(chipType);

  tags <- c("Mapping50K_Hind240"="_HIND", "Mapping50K_Xba240"="_XBA");
  tag <- tags[chipType];

  # Create URL
  urlPath <- getHapMapUrlPath(chipType);
  sprintf("%s/%s%s.CEL.gz", urlPath, sampleName, tag);
} # getHapMapUrl()


# downloadHapMapSample(dataSet, chipType=chipType, sampleName=sampleName);
downloadHapMapSample <- function(dataSet, tags=NULL, chipType, sampleName, ..., chipTypeAliases=NULL, gunzip=TRUE, skip=TRUE) {
  path <- getRawDataSetPath(dataSet=dataSet, tags=tags, chipType=chipType, ...);
  url <- getHapMapUrl(sampleName, chipType=chipType);
  filenameGZ <- basename(url);
  filename <- gsub("[.]gz$", "", filenameGZ);
  pathname <- file.path(path, filename);
  if (skip && isFile(pathname)) {
    return(pathname);
  }

  pathD <- dirname(path);
  pathnameGZ <- downloadFile(url, filename=filenameGZ, path=pathD);
  if (regexpr("[.]gz$", pathnameGZ) != -1) {
    gunzip(pathnameGZ);
  } 

  opwd <- getwd();
  on.exit(if (!is.null(opwd)) setwd(opwd));
  setwd(pathD);
  pathname <- arrangeCelFilesByChipType(filename, path=".", aliases=chipTypeAliases);
  setwd(opwd); opwd <- NULL;
  pathname <- file.path(path, filename);
  pathname;
} # downloadHapMapSample()


downloadHapMapSamples <- function(..., sampleNames) {
  pathnames <- sapply(sampleNames, FUN=function(sampleName) {
    downloadHapMapSample(..., sampleName=sampleName);
  });

  path <- dirname(pathnames[1]);
  getDataSet(path);
} # downloadHapMapSamples()


downloadHapMapDataSet <- function(dataSet, tags=NULL, chipType=c("GenomeWideSNP_6"), ..., chipTypeAliases=NULL, skip=TRUE) {
  # Argument 'chipType':
  chipType <- match.arg(chipType);

  population <- dataSet;
  dataSet <- paste(c("HapMap", population), collapse=",");
  path <- getRawDataSetPath(dataSet=dataSet, tags=tags, chipType=chipType, ...);
  # Already downloaded?
  ds <- getDataSet(path);
  if (skip && nbrOfFiles(ds) > 0) {
    return(ds);
  }

  urlPath <- getHapMapUrlPath(chipType);
  filenameD <- sprintf("%s.tgz", population);
  pathD <- dirname(path);
  url <- file.path(urlPath, filenameD);
  pathnameD <- downloadFile(url, filename=filenameD, path=pathD);
  untar(pathnameD, exdir=pathD);

  # Arrange by chip type
  opwd <- getwd();
  on.exit(if (!is.null(opwd)) setwd(opwd));
  setwd(pathD);
  pathnames <- list.files(pattern="[.](cel|CEL)$");
  pathnamesD <- arrangeCelFilesByChipType(pathnames, path=".", aliases=chipTypeAliases);
  setwd(opwd); opwd <- NULL;

  getDataSet(path);
} # downloadHapMapDataSet()




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# GEO
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
getGeoDataSetURL <- function(dataSet, ...) {
  sprintf("ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/series/%s/%s_RAW.tar", dataSet, dataSet);
} # getGeoDataSetURL()

getGeoDataFileURL <- function(sampleName, ext, ...) {
  # Argument 'sampleName':
  sampleName <- Arguments$getCharacter(sampleName);
  pattern <- "^GSM[0-9]+$";
  if (regexpr(pattern, sampleName) == -1) {
    throw("Invalid GEO sample name: ", sampleName);
  }

  # Argument 'ext':
  ext <- Arguments$getCharacter(ext);

  # Create GEO download URL
  urlRoot <- "ftp://ftp.ncbi.nlm.nih.gov/pub/geo/DATA/supplementary/samples";
  filename <- sprintf("%s.%s.gz", sampleName, ext);
  sampleSet <- gsub("[0-9]{3}$", "nnn", sampleName);
  url <- file.path(urlRoot, sampleSet, sampleName, filename);

  url;
} # getGeoDataSetURL()


downloadGeoRawDataSet <- function(dataSet, tags=NULL, chipType, ..., chipTypeAliases=NULL, url=getGeoDataSetURL(dataSet), gunzip=TRUE, skip=TRUE) {
  path <- getRawDataSetPath(dataSet=dataSet, tags=tags, chipType=chipType, ...);

  # Already downloaded?
  ds <- getDataSet(path);
  if (skip && nbrOfFiles(ds) > 0) {
    return(ds);
  }

  pathname <- downloadFile(url, path="downloads/", ...);
  pathD <- dirname(path);
  untar(pathname, exdir=pathD);
  if (gunzip) {
    pathnames <- list.files(path=pathD, pattern="[.]gz$", full.names=TRUE);
    sapply(pathnames, FUN=gunzip);
  }

  # Arrange by chip type
  opwd <- getwd();
  on.exit(if (!is.null(opwd)) setwd(opwd));
  setwd(pathD);
  pathnames <- list.files(pattern="[.](cel|CEL)$");
  pathnamesD <- arrangeCelFilesByChipType(pathnames, path=".", aliases=chipTypeAliases);
  setwd(opwd); opwd <- NULL;

  getDataSet(path);
} # downloadGeoRawDataSet()


downloadGeoRawDataFile <- function(dataSet, tags=NULL, chipType, sampleName, ext="CEL", url=getGeoDataFileURL(sampleName, toupper(ext)), ..., gunzip=TRUE, skip=TRUE) {
  # Argument 'url':
  url <- Arguments$getCharacter(url);

  path <- getRawDataSetPath(dataSet=dataSet, tags=tags, chipType=chipType, ...);

  filenameGZ <- basename(url);
  pathnameGZ <- file.path(path, filenameGZ);
  pathname <- gsub("[.]gz$", "", pathnameGZ);

  # Already downloaded?
  if (skip && isFile(pathname)) {
    return(pathname);
  }

  pathnameGZ <- tryCatch({
    suppressWarnings({
      downloadFile(url, path=path, ...);
    });
  }, error = function(ex) {
    urlL <- getGeoDataFileURL(sampleName, tolower(ext));
    downloadFile(urlL, path=path, ...);
  });

  if (gunzip) {
    gunzip(pathnameGZ);
    pathname <- gsub("[.]gz$", "", pathnameGZ);
  }

  pathname;
} # downloadGeoRawDataFile()


downloadGeoRawDataFiles <- function(..., sampleNames, skip=TRUE) {
  # Argument 'sampleNames':
  sampleNames <- Arguments$getCharacters(sampleNames);

  # Already downloaded?
  path <- getRawDataSetPath(...);
  if (skip) {
    ds <- getDataSet(path);
    hasAll <- all(is.element(sampleNames, getNames(ds)));
    if (hasAll) {
      return(ds);
    }
  }

  pathnames <- sapply(sampleNames, FUN=function(sampleName) {
    downloadGeoRawDataFile(..., sampleName=sampleName);
  });

  ds <- getDataSet(path);

  ds;
} # downloadGeoRawDataFiles()



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Affymetrix
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
downloadAffymetrixFile <- function(pathname, ..., urlRoot=Sys.getenv("AFFY_URLROOT", "http://www.affymetrix.com/Auth"), skip=TRUE) {
  # Argument 'pathname':
  filename <- basename(pathname);
  filename <- Arguments$getFilename(filename);

  # Download URL
  url <- file.path(urlRoot, pathname);
  downloadFile(url, filename=filename, path="downloads/", skip=skip);
} # downloadAffymetrixFile()


downloadAffymetrixNetAffxCsvFile <- function(pathname, ..., skip=TRUE) {
  # Argument 'pathname':
  filename <- basename(pathname);
  filename <- Arguments$getFilename(filename);
  pattern <- "^([^.]*)[.]((cn[.])*na[0-9]{2})[.].*.csv$";
  stopifnot(regexpr(pattern, filename) != -1L);

  # Infer chip type and NA tag
  chipType <- gsub(pattern, "\\1", filename);
  naTag <- gsub(pattern, "\\2", filename);

  tags <- sprintf(".%s", naTag);

  # Already available?
  pathnameD <- AffymetrixNetAffxCsvFile$findByChipType(chipType, tags=tags);
  if (skip && isFile(pathnameD)) {
    return(pathnameD);
  }

  # Destination path
  path <- file.path("annotationData", "chipTypes", chipType);
  path <- Arguments$getWritablePath(path);

  # Download file
  pathnameZ <- sprintf("%s.zip", pathname);
  pathnameZ <- downloadAffymetrixFile(pathnameZ, ...);

  # Unzip to destination
  unzip(pathnameZ, exdir=path);

  # Assert
  pathnameD <- AffymetrixNetAffxCsvFile$findByChipType(chipType, tags=tags);
  stopifnot(!is.null(pathnameD));

  pathnameD;
} # downloadAffymetrixNetAffxCsvFile()


downloadAffymetrixDataSet <- function(dataSet, tags=NULL, chipType=chipType, ..., skip=TRUE) {
  path <- getRawDataSetPath(dataSet, chipType=chipType);
  ds <- getDataSet(path);
  if (skip && nbrOfArrays(ds) > 0) {
    return(ds);
  }

  pathname <- downloadAffymetrixFile(..., skip=skip);
  pathnames <- unzip(pathname, exdir=path);

  # Did we extract a standalone root directory?  If so, "drop" the root.
  paths <- list.files(path=path, full.names=TRUE);
  if (length(paths) == 1 && isDirectory(paths)) {
    pathS <- paths[1];
    # Move all files in up one level
    pathnamesT <- gsub(pathS, path, pathnames, fixed=TRUE);
    for (kk in seq_along(pathnames)) {
      file.rename(from=pathnames[kk], to=pathnamesT[kk]);
    }
  }

  ds <- getDataSet(path);
  ds;
} # downloadAffymetrixDataSet();


############################################################################
# HISTORY:
# 2012-11-29
# o Added argument 'path' to downloadAffymetrixFile().
# 2012-11-28
# o Added downloadAffymetrixNetAffxCsvFile().
# 2012-09-12
# o BUG FIX: downloadGeoRawDataFile(..., gunzip=TRUE) would give an
#   error that it could not gunzip the downloaded file, iff the filename
#   extension was lower case, e.g. *.cel.gz instead of *.CEL.gz.
# 2012-09-02
# o Added downloadAffymetrixDataSet().
# o Added downloadGeoRawDataFile().
# 2012-09-01
# o Added argument 'chipTypeAliases' to several download functions.
# 2012-08-31
# o Added as a utility test script of the aroma.affymetrix package.
# 2012-08-23
# o Fixed up HapMap downloader.
# 2012-08-22
# o Added GEO and HapMap downloader.
# o Created.
############################################################################
