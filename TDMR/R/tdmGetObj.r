######################################################################################
#----- deprecated, only for downward compatibility. For actual runs, a simple
#-----   bst <- envT$bst
#-----   res <- envT$res
#----- should be sufficient.
#
# tdmGetObj: 
#
#- Helper fct to get \code{bst} or \code{res} object from \code{envT}.
#-   
#-   If \code{envObj} is not NULL, return envObj.
#-   If \code{envObj} is NULL, try to read the corresponding file from \code{theTuner/objFileName} ( 
#-   if dir \code{theTuner} does not exist, from \code{objFileName} in current dir) and return the data frame read. 
#-   If \code{envObj} is NULL, this function should only be called when tdm$fileMode==T (or was =T in the prior tuning run), 
#-   otherwise the files might be missing or contain old information.
#-
#-   @param envObj      object, either \code{envT$bst} or \code{envT$res}, or NULL
#-   @param objFileName alternative file to get the result from (.bst or .res file)
#-   @param theTuner    alternative dir where to look for objFileName
#-   @param tdm         here only needed for tdm$optsVerbosity
#-
#-   @return envObj
#-
#- @seealso \code{\link{unbiasedRun}}, \code{\link{tdmBigLoop}}
#- @export
######################################################################################
tdmGetObj <- function(envObj,objFileName, theTuner,tdm) {
  if (is.null(envObj)) {
    tunedir = ifelse(file.exists(theTuner),paste(theTuner,"/",sep=""), "");
    # downward compatibility: if subdir envT$theTuner (e.g. "spot") is not present, take the current dir
    tBstFile = paste(tunedir,objFileName,sep="");
    if (!file.exists(tBstFile))
      stop(paste("Could not find file",tBstFile,"in directory",getwd()));
    suffix = unlist(strsplit(objFileName,".",fixed=T))[2]
  	if (tdm$optsVerbosity>0) writeLines(paste("Loading", suffix, "file data from:", tBstFile,"in directory",getwd()), con=stderr());
    obj <- read.table(tBstFile, sep=" ", header = TRUE);
  } else {
    obj <- envObj;
  }
  obj;
}

