###########################################################################/**
# @RdocClass SNPsNormalization
#
# @title "The SNPsNormalization class"
#
# \description{                                                                                                                                                                             
#  @classhierarchy
#
#  This class represents the SNPs normalization method [1], which 
#  scales copy number values SNP by SNP.
# }
# 
# @synopsis 
#
# \arguments{
#   \item{data}{A named @list with data set named \code{"total"} and
#     \code{"fracB"} where the former should be of class
#     @see "aroma.core::AromaUnitTotalCnBinarySet" and the latter of
#     class @see "aroma.core::AromaUnitFracBCnBinarySet".  The
#     two data sets must be for the same chip type, have the same
#     number of samples and the same sample names.}
#   \item{tags}{Tags added to the output data sets.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
#
# \references{
#   [1] ...
# }
#
# \seealso{
#   Low-level versions of the SNPs normalization method is available
#   via @see "snpsNByTotalAndFracB.matrix" methods.
# }
#
#*/###########################################################################
setConstructorS3("SNPsNormalization", function(data=NULL, tags="*", ...) {
  # Validate arguments
  if (!is.null(data)) {
    if (!is.list(data)) {
      throw("Argument 'data' is not a list: ", class(data)[1]);
    }
    reqNames <- c("total", "fracB");
    ok <- is.element(reqNames, names(data));
    if (!all(ok)) {
      throw(sprintf("Argument 'data' does not have all required elements (%s): %s", paste(reqNames, collapse=", "), paste(reqNames[!ok], collapse=", ")));
    }
    data <- data[reqNames];

    # Assert correct classes
    className <- "AromaUnitTotalCnBinarySet";
    ds <- data$total;
    if (!inherits(ds, className)) {
      throw(sprintf("The 'total' data set is not of class %s: %s", className, class(ds)[1]));
    }

    className <- "AromaUnitFracBCnBinarySet";
    ds <- data$fracB;
    if (!inherits(ds, className)) {
      throw(sprintf("The 'fracB' data set is not of class %s: %s", className, class(ds)[1]));
    }

    # Assert that the chip types are compatile
    if (getChipType(data$total) != getChipType(data$fracB)) {
      throw("The 'total' and 'fracB' data sets have different chip types: ", 
            getChipType(data$total), " != ", getChipType(data$fracB));
    }

    # Assert that the data sets have the same number data files
    nbrOfFiles <- nbrOfFiles(data$total)
    if (nbrOfFiles != nbrOfFiles(data$fracB)) {
      throw("The number of samples in 'total' and 'fracB' differ: ", 
            nbrOfFiles, " != ", nbrOfFiles(data$fracB));
    }

    # Assert that the data sets have the same samples
    if (!identical(getNames(data$total), getNames(data$fracB))) {
      throw("The samples in 'total' and 'fracB' have different names.");
    }
  }

  # Arguments '...':
  args <- list(...);
  if (length(args) > 0) {
    argsStr <- paste(names(args), collapse=", ");
    throw("Unknown arguments: ", argsStr);
  }

  this <- extend(Object(...), "SNPsNormalization",
    .data = data
  );
  setTags(this, tags);

  this; 
})


setMethodS3("as.character", "SNPsNormalization", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);

  dsList <- getDataSets(this);
  s <- c(s, sprintf("Data sets (%d):", length(dsList)));
  for (kk in seq(along=dsList)) {
    ds <- dsList[[kk]];
    s <- c(s, sprintf("<%s>:", capitalize(names(dsList)[kk])));
    s <- c(s, as.character(ds));
  } 
 
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)



setMethodS3("getAsteriskTags", "SNPsNormalization", function(this, collapse=NULL, ...) {
  tags <- "NSAN,SNPN";

  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  }
  
  tags;
}, private=TRUE) 


setMethodS3("getName", "SNPsNormalization", function(this, ...) {
  dsList <- getDataSets(this);
  ds <- dsList$total;
  getName(ds);
}) 



setMethodS3("getTags", "SNPsNormalization", function(this, collapse=NULL, ...) {
  # "Pass down" tags from input data set
  dsList <- getDataSets(this);
  ds <- dsList$total;
  tags <- getTags(ds, collapse=collapse);

  # Get class-specific tags
  tags <- c(tags, this$.tags);

  # Update default tags
  tags[tags == "*"] <- getAsteriskTags(this, collapse=",");

  # Collapsed or split?
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  } else {
    tags <- unlist(strsplit(tags, split=","));
  }

  if (length(tags) == 0)
    tags <- NULL;

  tags;
})


setMethodS3("setTags", "SNPsNormalization", function(this, tags="*", ...) {
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));
    tags <- tags[nchar(tags) > 0];
  }
  
  this$.tags <- tags;
})


setMethodS3("getFullName", "SNPsNormalization", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})


setMethodS3("getDataSets", "SNPsNormalization", function(this, ...) {
  this$.data;
})
 

setMethodS3("getRootPath", "SNPsNormalization", function(this, ...) {
  "totalAndFracBData";
})


setMethodS3("getPath", "SNPsNormalization", function(this, create=TRUE, ...) {
  # Create the (sub-)directory tree for the data set

  # Root path
  rootPath <- getRootPath(this);

  # Full name
  fullname <- getFullName(this);

  # Chip type    
  dsList <- getDataSets(this);
  ds <- dsList$total;
  chipType <- getChipType(ds, fullname=FALSE);

  # The full path
  path <- filePath(rootPath, fullname, chipType, expandLinks="any");

  # Verify that it is not the same as the input path
  inPath <- getPath(ds);
  if (getAbsolutePath(path) == getAbsolutePath(inPath)) {
    throw("The generated output data path equals the input data path: ", path, " == ", inPath);
  }

  # Create path?
  if (create) {
    if (!isDirectory(path)) {
      mkdirs(path);
      if (!isDirectory(path))
        throw("Failed to create output directory: ", path);
    }
  }

  path;
})


setMethodS3("nbrOfFiles", "SNPsNormalization", function(this, ...) {
  dsList <- getDataSets(this);
  ds <- dsList$total;
  nbrOfFiles(ds);
})


setMethodS3("getOutputDataSets", "SNPsNormalization", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 

  res <- this$.outputDataSets;
  if (is.null(res)) {
    res <- allocateOutputDataSets(this, ..., verbose=less(verbose, 10));
    this$.outputDataSets <- res;
  }
  res;
}) 


setMethodS3("allocateOutputDataSets", "SNPsNormalization", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 

  verbose && enter(verbose, "Retrieve/allocation output data sets");

  dsList <- getDataSets(this);
  path <- getPath(this);

  res <- list();
  for (kk in seq(along=dsList)) {
    ds <- dsList[[kk]];
    verbose && enter(verbose, sprintf("Data set #%d ('%s') of %d", 
                                 kk, getName(ds), length(dsList)));

    for (ii in seq(ds)) {
      df <- getFile(ds, ii);
      verbose && enter(verbose, sprintf("Data file #%d ('%s') of %d", 
                                        ii, getName(df), nbrOfFiles(ds)));

      filename <- getFilename(df);
      pathname <- Arguments$getWritablePathname(filename, path=path, 
                                                       mustNotExist=FALSE);
      # Skip?
      if (isFile(pathname)) {
        verbose && cat(verbose, "Already exists. Skipping.");
        verbose && exit(verbose);
        next;
      }

      # Create temporary file
      pathnameT <- sprintf("%s.tmp", pathname);
      pathnameT <- Arguments$getWritablePathname(pathnameT, mustNotExist=TRUE);

      # Copy source file
      copyFile(getPathname(df), pathnameT);

      # Make it empty by filling it will missing values
      # AD HOC: We should really allocate from scratch here. /HB 2010-06-21
      # if they are "fracB" they can be just copied and leave like that.
      
      dfT <- newInstance(df, pathnameT);
      if (kk == 1)
        dfT[,1] <- NA;

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Renaming temporary file
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Renaming temporary output file");
      file.rename(pathnameT, pathname);
      if (!isFile(pathname)) {
        throw("Failed to rename temporary file ('", pathnameT, "') to final file ('", pathname, "')");
      }
      verbose && exit(verbose);

      verbose && cat(verbose, "Copied: ", pathname);
      verbose && exit(verbose);
    } # for (ii ...)

    dsOut <- byPath(ds, path=path, ...);

    # AD HOC: The above byPath() grabs all *.asb files. /HB 2010-06-20
##    keep <- is.element(sapply(dsOut, getFilename), sapply(ds, getFilename));
    filenamesOut <- basename(getPathnames(dsOut));
    filenames <- basename(getPathnames(ds));
    keep <- is.element(filenamesOut, filenames);
    dsOut <- extract(dsOut, keep);

    res[[kk]] <- dsOut;
    rm(ds, dsOut);
    
    verbose && exit(verbose);
  } # for (kk ...)

  names(res) <- names(dsList);

  this$.outputDataSets <- res;

  verbose && exit(verbose);

  res;
}, protected=TRUE)

setMethodS3("findUnitsTodo", "SNPsNormalization", function(this, units, ..., verbose=FALSE) {
  # Argument 'verbose':

  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 

#  verbose && enter(verbose, "Finding units to do");
  dsList <- getOutputDataSets(this);

  # The last data set is always updated last
  ds <- dsList[[1]];
  verbose && print(verbose, ds);

  # The last file (in lexicographic ordering) is always updated last
  fullnames <- getFullNames(ds);
  verbose && str(verbose, fullnames);

  o <- order(fullnames, decreasing=TRUE);
  idx <- o[1];
  df <- getFile(ds, idx);
  verbose && print(verbose, df);

  # Read all values
  values <- df[units,1,drop=TRUE];

  verbose && cat(verbose, "Number of units: ", length(values));

  # Identify all missing values
  nok <- is.na(values);
  units <- which(nok);
  verbose && printf(verbose, "Number of units to do: %d (%.2f%%)\n", 
                       length(units), 100*length(units)/length(values));

  verbose && cat(verbose, "Units to do (with missing values):");
  verbose && str(verbose, units);

 # verbose && exit(verbose);

  units;
})

###########################################################################/**
# @RdocMethod process
#
# @title "Scaling of the copy number values SNP by SNP"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{this}(Object gerenated by SNPsNormalization.)
#   \item{references}{@logical or @numeric @matrix saying which samples should be considered as normal, previously
#                     calculated by @see "NSANormalization.".}
#   \item{units}{@numeric @vector indicating the specific units to scale. The default value is "remaining".}
#   \item{force}{@logical flag indicating if the already scaled units have to be scaled again. Initially set to FALSE.}
#   \item{...}{Additional arguments passed to @see "SNPsNormalization".}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @double @vector.
# }
#   
# \seealso{
#   @seeclass
# }
#*/###########################################################################

setMethodS3("process", "SNPsNormalization", function(this, references = NULL, units="remaining", ..., force=FALSE, ram=NULL, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  dsList <- getDataSets(this);
  dsTCN <- dsList$total;
  
  # Argument 'references':
  if (!is.list(references)) {
    throw("Argument 'references' is not a list: ", class(references)[1]);
  }                           
  reqName <- c("normalReg");
  ok <- is.element(reqName, names(references));
  if (!all(ok)) {
    throw(sprintf("Argument 'references' does not have all required elements (%s): %s", paste(reqName, collapse=", "), paste(reqName[!ok], collapse=", ")));
  }
  references <- references[reqName];

  className <- "AromaUnitFracBCnBinarySet";
  ds <- references$normalReg;
  if (!inherits(ds, className)){
    throw(sprintf("The 'normalReg' data set is not of class %s: %s", className, class(ds)[1]));
  }
  if(nbrOfFiles(dsTCN) != nbrOfFiles(ds)){
    throw("The number of samples in 'data' and 'references' differ: ", 
          nbrOfFiles, " != ", nbrOfFiles(ds));  
  }
  if (!identical(getNames(dsTCN), getNames(references$normalReg))) {
    throw("The samples in 'data' and 'references' have different names.");
  }

  rsNR <- references$normalReg;
  
  # Argument 'units':
  df <- getFile(dsTCN, 1);
  nbrOfUnits <- nbrOfUnits(df);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Finding units to do
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Finding Units To Do");
  if (is.null(units) || units == "remaining") {
    units <- seq(length=nbrOfUnits);
  }else{
    units <- Arguments$getIndices(units, max=nbrOfUnits);
  }
  if(!force){
    units <- findUnitsTodo(this, units, verbose=verbose);
  }
  nbrOfUnits <- length(units);
  
  verbose && exit(verbose);

  # Argument 'ram':
  ram <- getRam(aromaSettings, ram);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "SNPs normalization of ASCNs");
  nbrOfFiles <- nbrOfFiles(this);
  verbose && cat(verbose, "Number of arrays: ", nbrOfFiles);
  verbose && printf(verbose, "Number of units to do: %d (%.2f%%)\n", 
                      length(units), 100*length(units)/nbrOfUnits(df));
  verbose && cat(verbose, "Units:");
  verbose && str(verbose, units);

  chipType <- getChipType(dsTCN, fullname=FALSE);
  verbose && cat(verbose, "Chip type: ", chipType);
  rm(dsList);

  sampleNames <- getNames(dsTCN);
  dimnames <- list(NULL, sampleNames, c("total", "fracB"));

  outPath <- getPath(this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocate output data sets
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- getOutputDataSets(this, verbose=less(verbose, 5));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Process in chunks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calculating number of units to fit per chunk");
  verbose && cat(verbose, "RAM scale factor: ", ram);

  bytesPerChunk <- 100e6;       # 100Mb
  verbose && cat(verbose, "Bytes per chunk: ", bytesPerChunk);

  bytesPerUnitAndArray <- 2*8;  # Just a rough number; good enough?
  verbose && cat(verbose, "Bytes per unit and array: ", bytesPerUnitAndArray);

  bytesPerUnit <- nbrOfFiles * bytesPerUnitAndArray;
  verbose && cat(verbose, "Bytes per unit: ", bytesPerUnit);

  unitsPerChunk <- ram * bytesPerChunk / bytesPerUnit;
  unitsPerChunk <- as.integer(max(unitsPerChunk, 1));
  unitsPerChunk <- min(unitsPerChunk, nbrOfUnits);
  verbose && cat(verbose, "Number of units per chunk: ", unitsPerChunk);

  nbrOfChunks <- ceiling(nbrOfUnits / unitsPerChunk);
  verbose && cat(verbose, "Number of chunks: ", nbrOfChunks);
  verbose && exit(verbose);

  idxs <- 1:nbrOfUnits;
  head <- 1:unitsPerChunk;

  count <- 1;
  while (length(idxs) > 0) {

    verbose && enter(verbose, sprintf("Chunk #%d of %d", count, nbrOfChunks));
    if (length(idxs) <= unitsPerChunk) {
      head <- 1:length(idxs);
    }
    uu <- idxs[head];
    verbose && cat(verbose, "Units: ");
    verbose && str(verbose, units[uu]);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Reading (total,normalReg) data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Reading (total) data");
    total <- extractMatrix(dsTCN, units=units[uu], verbose=less(verbose,25));
    verbose && str(verbose, total);
    verbose && exit(verbose);
        
    verbose && enter(verbose, "Reading references");
    normalReg <- extractMatrix(rsNR, units=units[uu], verbose=less(verbose,25));
    verbose && str(verbose, normalReg);
    verbose && exit(verbose);
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Normalizing by SNPs
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    verbose && enter(verbose, "Normalizing");
    dataN <- snpsNByTotalAndFracB(total, references=normalReg, ..., 
             verbose=less(verbose,5));

    fit <- attr(dataN, "modelFit");
    verbose && str(verbose, fit);
    verbose && str(verbose, dataN);
    verbose && exit(verbose);

    rm(normalReg);  # Not needed anymore
    gc <- gc();
    verbose && print(verbose, gc);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Storing normalized data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Storing normalized data");

    ds <- res[[1]];
    verbose && enter(verbose, sprintf("Data set ('%s')", 
                                         getName(ds)));

    # Store in lexicograph ordering
    fullnames <- getFullNames(ds);
    idxss <- order(fullnames, decreasing=FALSE);
    for (ii in idxss) {
      df <- getFile(ds, ii);
      verbose && enter(verbose, sprintf("Data file #%d ('%s') of %d", 
                                      ii, getName(df), nbrOfFiles(ds)));

      signals <- dataN[,ii];
      verbose && cat(verbose, "Signals:");
      verbose && str(verbose, signals);
      df[units[uu],1] <- signals;
      rm(signals);

      verbose && exit(verbose);
    } # for (ii ...)
    verbose && exit(verbose);

    verbose && exit(verbose);
    # Next chunk
    idxs <- idxs[-head];
    count <- count + 1;

    verbose && exit(verbose);
  } # while(length(idxs) > 0)

  verbose && exit(verbose);
  invisible(res);
})


############################################################################
# HISTORY:
# 2010-06-28 [MO]
# o Created from CalMaTeNormalization.R.
############################################################################
