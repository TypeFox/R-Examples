#'Stream Subspace Instances from Disc
#'
#'A DSD object to stream Data Points from a .arff file.
#'@param file a string that contains the path to the file
#'@export
DSD_SubspaceARFFStream <- function(file) {
  jref <- rJava::.jcall("DataStreamBuilder",
                        returnSig="Lmoa/streams/clustering/SubspaceARFFStream;",
                        method="buildSubspaceARFFStream",
                        file)
  
  numAtts <- rJava::.jcall("DataStreamAccessor","I","getNumAtts",jref)
  res <- structure(list(description="A data stream read from a file.",
                        javaObj=jref,d=numAtts
  ),
  class=c("DSD_SubspaceARFFStream","DSD_SubspaceMOA","DSD"))
  
  return(res)
}