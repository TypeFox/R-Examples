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
pathT <- file.path(path, "R/");
sourceDirectory(path=pathT);

anTags <- c("dens", "hets"); 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calibrated or not?
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
calTag <- textMenu(calTags, title="Choose calibration method:", value=TRUE);
if (calTag == "<none>") calTag <- NULL;



dsList <- loadSets(dataSet, tags=c(tags, calTag), chipType=chipType, verbose=verbose);
verbose && print(verbose, dsList);

# Sanity check
# stopifnot(all(sapply(dsList, FUN=length) == 40));


# Translated data set tags
dsTags <- tags(fullname(dataSet, getTags(dsList[[1]])));
dsTags <- gsub("ACC,-XY,BPN,-XY,RMA,FLN,-XY", "ASCRMAv2", dsTags);
dsTags <- gsub("ACC,-XY,BPN,-XY,AVG,FLN,-XY", "ASCRMAv2", dsTags);
dsTags <- gsub("CMTN", "CalMaTe", dsTags);
dsTags <- dropTags(dsTags, drop="refs=N");


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extracting one tumor-normal pair
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
verbose && enter(verbose, "Extracting tumor-normal pair of interest");
pairs <- getPairs(dsList);
pair <- pairs[sampleName,,drop=TRUE];
verbose && cat(verbose, "Tumor: ", pair["tumor"]);
verbose && cat(verbose, "Normal: ", pair["normal"]);
pairName <- paste(pair, collapse="vs"); 
verbose && cat(verbose, "Pair: ", pairName);
verbose && exit(verbose);
  

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segments to be studied
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for (kk in seq(along=segments)) {
  segName <- names(segments)[kk];
  segment <- segments[[segName]];
  chr <- segment[1];
  region <- segment[2:3];

  chrTag <- sprintf("chr%02d", chr);
  segTag <- sprintf("%s:%g-%gMb", chrTag, region[1], region[2]);
  
  verbose && enter(verbose, sprintf("Segment #%d ('%s') of %d", kk, segTag, length(segments)));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extracting signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dataList <- lapply(pair, FUN=function(name) {
    extractSignals(dsList, sampleName=name, chromosome=chr, region=region*1e6, verbose=verbose);
  });
  names(dataList) <- names(pair);

  # Extract (gammaT, gammaN)
  gammaT <- getSignals(dataList$tumor$tcn);
  gammaN <- getSignals(dataList$normal$tcn);
  
  # Extract (betaT, betaN)
  betaT <- getSignals(dataList$tumor$baf);
  betaN <- getSignals(dataList$normal$baf);
  muN <- callNaiveGenotypes(betaN);
  x <- dataList$normal$tcn$x / 1e6;

  col <- "black";
  if (is.element("hets", anTags)) {
    col <- rep(col, times=length(muN));
    col[muN != 1/2] <- hetCol;
  }

  for (tbn in c(FALSE, TRUE)) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Plot (betaN,betaT)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    beta <- cbind(N=betaN, T=betaT);

    tagsT <- dsTags;
    if (tbn) {
      verbose && enter(verbose, "TumorBoosting");
      betaTN <- normalizeTumorBoost(betaT=betaT, betaN=betaN);
      verbose && str(verbose, betaTN);
      verbose && exit(verbose);
      beta[,"T"] <- betaTN;
      tagsT <- fullname(tagsT, "TBN");
    }

    ascn <- cbind(A=(1-beta[,"T"])*gammaT, B=beta[,"T"]*gammaT);

    tags <- c(gsub(":","_",segTag), tagsT, "BvsB", anTags);
    toPNG(name=sampleName, tags=tags, width=840, {
      plotBvsB(beta, col=col, sampleName=sampleName, segName=segName, segTag=segTag, dataSet=dataSet, tagsT=tagsT, chipType=chipType);
    })

    tags <- c(gsub(":","_",segTag), tagsT, "ASCN", anTags);
    toPNG(name=sampleName, tags=tags, width=840, {
      plotASCN(ascn, col=col, sampleName=sampleName, segName=segName, segTag=segTag, dataSet=dataSet, tagsT=tagsT, chipType=chipType);
    })

    tags <- c(gsub(":","_",segTag), tagsT, "BvsC", anTags);
    toPNG(name=sampleName, tags=tags, width=840, {
      plotCvsB(ascn, col=col, sampleName=sampleName, segName=segName, segTag=segTag, dataSet=dataSet, tagsT=tagsT, chipType=chipType);
    })


    tags <- c(chrTag, tagsT, "TCN", anTags);
    toPNG(name=sampleName, tags=tags, width=840, height=300, {
      plotTrackTCN(gammaT, x=x, cex=2, col=col, sampleName=sampleName, chrTag=chrTag, dataSet=dataSet, tagsT=tagsT, chipType=chipType);
    });

    tags <- c(chrTag, tagsT, "BAFT", anTags);
    toPNG(name=sampleName, tags=tags, width=840, height=300, {
      plotTrackBAF(beta[,"T"], x=x, cex=2, col=col, sampleName=sampleName, chrTag=chrTag, dataSet=dataSet, tagsT=tagsT, chipType=chipType);
    });

    tags <- c(chrTag, tagsT, "BAFN", anTags);
    toPNG(name=sampleName, tags=tags, width=840, height=300, {
      plotTrackBAF(beta[,"N"], x=x, cex=2, col=col, sampleName=sampleName, chrTag=chrTag, dataSet=dataSet, tagsT=tagsT, chipType=chipType);
    });

    tags <- c(chrTag, tagsT, "DH", anTags);
    toPNG(name=sampleName, tags=tags, width=840, height=300, {
      plotTrackDH(beta[,"T"], x=x, muN=muN, cex=2, col=col, sampleName=sampleName, chrTag=chrTag, dataSet=dataSet, tagsT=tagsT, chipType=chipType);
    });

    tags <- c(gsub(":","_",segTag), tagsT, "C1C2", anTags);
    toPNG(name=sampleName, tags=tags, width=840, height=300, {
      plotTrackC1C2(ascn, x=x, muN=muN, cex=2, sampleName=sampleName, segName=segName, segTag=segTag, dataSet=dataSet, tagsT=tagsT, chipType=chipType);
    })

    tags <- c(gsub(":","_",segTag), tagsT, "C1C2", "smooth", anTags);
    toPNG(name=sampleName, tags=tags, width=840, height=300, {
      plotTrackC1C2(ascn, x=x, muN=muN, cex=2, smooth="2Mb", sampleName=sampleName, segName=segName, segTag=segTag, dataSet=dataSet, tagsT=tagsT, chipType=chipType);
    })
  } # for (tbn ...)

  verbose && exit(verbose);
} # for (kk ...)



###########################################################################
# HISTORY:
# 2011-03-09 [HB]
# o Created from PN's code in the online vignette.
###########################################################################
