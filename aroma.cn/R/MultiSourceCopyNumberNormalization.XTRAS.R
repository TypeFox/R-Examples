# \examples{\dontrun{
#  par(mar=c(3.2,4.2,4.0,1.1), mgp=c(3,0.6,0), cex.axis=1.4, tcl=-0.2);
#  plotSmoothedPairsOne(mscn, sampleName="TCGA-09-1670");
# }}

setMethodS3("plotSmoothedPairsOne", "MultiSourceCopyNumberNormalization", function(this, sampleName, units=NULL, ..., Mlim=c(-3,3), gap=0.1, cex.labels=1.5, oma=c(1,1,1,1)*2, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  panel <- function(x, y, cex=1, ...) { 
     xx <- M[,x,drop=TRUE];
     yy <- M[,y,drop=TRUE];

     abline(a=0, b=1, lty=2); 
     points(xx,yy, cex=cex, ...);

     xx <- fit$s[,x,drop=TRUE];
     yy <- fit$s[,y,drop=TRUE];

     points(xx,yy, col="red", cex=1.5*cex, ...);
   } # panel()

    
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'sampleName':
  sampleName <- Arguments$getCharacter(sampleName);
  allSampleNames <- getAllNames(this);
  if (!is.element(sampleName, allSampleNames)) {
    throw("No such sample: ", sampleName);
  }

  # Argument 'fit':


  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 

  # Argument 'force':
  force <- Arguments$getLogical(force);



  verbose && enter(verbose, "Plot smoothed pairs for sample"); 

  verbose && enter(verbose, "Get smoothed data");
  dsList <- getSmoothedDataSets(this, verbose=less(verbose, 5));
  verbose && print(verbose, dsList);
  verbose && exit(verbose); 

  verbose && enter(verbose, "Generating data source labels");
  tags <- Reduce(intersect, lapply(dsList, FUN=getTags));
  sites <- sapply(dsList, FUN=function(ds) setdiff(getTags(ds), tags));
  if (is.matrix(sites)) sites <- sites[1,];
  platforms <- sapply(dsList, FUN=function(ds) {
    footer <- readFooter(getFile(ds, 1));
    footer$sourceDataFile$platform;
  });
  chipTypes <- sapply(dsList, FUN=function(ds) {
    footer <- readFooter(getFile(ds, 1));
    footer$sourceDataFile$chipType;
  });
  verbose && exit(verbose); 
 
  verbose && enter(verbose, "Extracting data files");
  dfList <- extractTupleOfDataFiles(this, dsList=dsList, name=sampleName, 
                                                 verbose=less(verbose, 1));
  verbose && print(verbose, dfList);

  fullnames <- sapply(dfList, FUN=getFullName);
  fullnames <- unname(fullnames);
  fullnames <- gsub(",ratios", "", fullnames);
  fullnames <- gsub(",log2ratio", "", fullnames);
  fullnames <- gsub(",total", "", fullnames);
  verbose && cat(verbose, "Full names:");
  verbose && print(verbose, fullnames);

  verbose && exit(verbose);

  verbose && enter(verbose, "Get UGP file for smoothed data");
  # All smoothed data sets have the same UGP file
  ugp <- getAromaUgpFile(dfList[[1]]);
  verbose && print(verbose, ugp);
  verbose && exit(verbose); 
  
  if (is.null(units)) {
    verbose && enter(verbose, "Identifying units on autosomal chromosomes");
    units <- getUnitsOnChromosomes(ugp, chromosomes=1:22, drop=TRUE);
    verbose && cat(verbose, "Units on autosomal chromosomes:");
    verbose && str(verbose, units);
    verbose && exit(verbose); 
  }

  verbose && enter(verbose, "Fitting MSCN model");
  fit <- fitOne(this, dfList=dfList, verbose=less(verbose, 1));
  verbose && str(verbose, fit);
  verbose && exit(verbose); 

  verbose && enter(verbose, "Extract log2 CN ratios");
  M <- sapply(dfList, FUN=function(df) {
    extractMatrix(df, units=units, drop=TRUE);
  });
  verbose && str(verbose, M);
  stopifnot(is.matrix(M));
  verbose && exit(verbose);

  verbose && enter(verbose, "Generating data source labels");
  colnames(M) <- paste(platforms, chipTypes, sep="\n"); 
  verbose && print(verbose, colnames(M));
  verbose && exit(verbose);

  legend <- c("Data sources:", fullnames);
  legend <- c(legend, sprintf("Number of data points per panel: %d", length(units)));
  n <- length(legend);
  legend <- paste(legend, collapse="\n");

  verbose && enter(verbose, "Plotting");
  I <- matrix(seq_len(ncol(M)), nrow=1);
  colnames(I) <- colnames(M);
  pairs(I, pch=".", cex=2, lower.panel=NULL, upper.panel=panel, 
       xlim=Mlim, ylim=Mlim, gap=gap, cex.labels=cex.labels, oma=oma, ...);
  stext(side=1, pos=0.05, line=-3, margin=c(0,0), cex=1.8, sampleName);
  stext(side=1, pos=0.05, line=-0.01, margin=c(0,-3), cex=0.7, legend);
  verbose && exit(verbose);

  verbose && exit(verbose);
}, protected=TRUE)



############################################################################
# HISTORY:
# 2010-01-14
# o Added getBacktransforms() and plotBacktransforms().
# o Added plotSmoothedPairsOne() to MultiSourceCopyNumberNormalization.
############################################################################
