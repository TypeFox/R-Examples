###########################################################################
# Title:
# Author: Henrik Bengtsson
###########################################################################


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Loading support files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Find the pathname and directory of this script
library("R.utils");
pathname <- names(findSourceTraceback())[1];
path <- dirname(pathname);

# Loading include files
sourceTo("001.include.R", path=path);
#sourceTo("002.datasets.R", path=path);


figPath <- file.path("figures");
figPath <- Arguments$getWritablePath(figPath);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Choose data set & chip type
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rootPath <- "totalAndFracBData";
rootPath <- Arguments$getReadablePath(rootPath);
dataSets <- list.files(path=rootPath, full.names=FALSE, all.files=FALSE);

# Skip already process data sets
##keep <- grep(",CMTN", dataSets, invert=TRUE);
## dataSets <- dataSets[keep];

if (length(chipTypes) == 0) {
  throw("No data sets available: ", path);
}
dataSet <- textMenu(dataSets, title="Choose data set:", value=TRUE);


path <- file.path(rootPath, dataSet);
path <- Arguments$getReadablePath(path, mustExist=TRUE);
chipTypes <- list.files(path=path, full.names=FALSE, all.files=FALSE);
if (length(chipTypes) == 0) {
  throw("No chip types available for this data set: ", path);
}
if (length(chipTypes) == 1) {
  chipType <- chipTypes[1];
} else {
  chipType <- textMenu(chipTypes, title="Choose chip type:", value=TRUE);
}


dsList <- loadSets(dataSet, chipType=chipType, verbose=verbose);
verbose && print(verbose, dsList);

dsTags <- getTags(dsList$total, collapse=",");
dsTags <- gsub("ACC(|,ra)-XY,BPN,-XY,RMA,FLN,-XY", "ASCRMAv2", dsTags);
dsTags <- gsub("ACC(|,ra),-XY,BPN,-XY,AVG,FLN,-XY", "ASCRMAv2", dsTags);
dsTags <- gsub("CMTN", "CalMaTe", dsTags);
  

ugp <- getAromaUgpFile(dsList$total);
print(ugp);

chromosomes <- getChromosomes(ugp)[19:20];
print(chromosomes);

sampleNames <- getNames(dsList$total)[1];
print(sampleNames);

col <- "#999999";

dataSetT <- paste(c(getName(dsList[[1]]), dsTags), collapse=",");
figPathT <- file.path(figPath, dataSetT, chipType);
figPathT <- Arguments$getWritablePath(figPathT);

for (ii in seq(along=sampleNames)) {
  sampleName <- sampleNames[ii];
  verbose && enter(verbose, sprintf("Sample #%d ('%s') of %d", ii, sampleName, length(sampleNames)));

  for (kk in seq(along=chromosomes)) {
    chr <- chromosomes[kk];
    chrTag <- sprintf("chr%02d", chr);
  
    verbose && enter(verbose, sprintf("Chromosome #%d ('%s') of %d", kk, chrTag, length(chromosomes)));
  
    data <- extractSignals(dsList, sampleName=sampleName, chromosome=chr, verbose=verbose); 
    verbose && print(verbose, data);
  
    # Extract (gamma, beta)
    gamma <- getSignals(data$tcn);
    beta <- getSignals(data$baf);
    x <- data$tcn$x / 1e6;

    fullname <- fullname(sampleName, chrTag, "TCN");
    devEval("png", name=fullname, width=840, height=300, {
      plotTrackTCN(gamma, x=x, col=col, sampleName=sampleName, chrTag=chrTag, dataSet=dataSet, chipType=chipType);
    }, path=figPathT);

    fullname <- fullname(sampleName, chrTag, "BAF");
    devEval("png", name=fullname, width=840, height=300, {
      plotTrackBAF(beta, x=x, col=col, sampleName=sampleName, chrTag=chrTag, dataSet=dataSet, chipType=chipType);
    }, path=figPathT);
 
    verbose && exit(verbose);
  } # for (kk ...)

  verbose && exit(verbose);
} # for (ii ...)

 

###########################################################################
# HISTORY:
# 2011-08-01 [HB]
# o Now the figures are written to a data-set specific directory.
# o BUG FIX: Forgot to define 'col'.
# 2011-03-14 [HB]
# o Created.
###########################################################################
