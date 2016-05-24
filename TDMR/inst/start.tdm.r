###############################################################################################
# start.tdm.r: 
#     Script to load TDMR software. The list tdm should be def'd by the user beforehand.
#     If it does not exist, a minimal list tdm will be def'd here. 
#     (TDMR is either sourced from R files in tdm$tdmPath or its library is loaded, if is.null(tdm$tdmPath))
#
#     start.tdm.r is not needed for the normal TDMR user (who loads the library TDMR). It is only
#     needed for the TDMR developer who wants to start TDMR from R sources. 
#
###############################################################################################

#
# script part to initialize tdm
# (this part is now normally done in tdmBigLoop, it is only needed here in case of a 'developer start' (from R sources)
#
if (is.na(match("tdm",ls()))) tdm <- list();
if (is.null(tdm$theSpotPath)) tdm$theSpotPath <- NA;
if (is.null(tdm$theRsfaPath)) tdm$theRsfaPath <- NA;
if (is.null(tdm$parallelCPUs)) tdm$parallelCPUs=1;


#
# function part to load TDMR functions or library and SPOT functions or library
#
if (is.null(tdm$tdmPath)) {
    cat("Loading installed library TDMR \n");
    require(TDMR);
    #if (tdm$parallelCPUs>1) sfLibrary(TDMR); # --- obsolete now, we load library(TDMR) in bigLoopStep (tdmBigLoop.r)
} else {
    cat("Sourcing TDMR from R files in",tdm$tdmPath,"\n");
    source(paste(tdm$tdmPath,"inst/source.tdm.r",sep="/"),local=TRUE);
    source.tdm(tdm$tdmPath,tdmParallelCPUs=tdm$parallelCPUs,theSpotPath=tdm$theSpotPath,theRsfaPath=tdm$theRsfaPath);
    if (tdm$parallelCPUs>1) {
     stop("For tdm$parallelCPUs>1 it is required to use the *library* version of TDMR. Consider to set tdm$tdmPath=NULL.")
     #sfExport(list=c("tdmRandomSeed"));      # --- obsolete now
    }
}
