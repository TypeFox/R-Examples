#'@export
get_points.DSD_SubspaceMOA <- function(x, n=1,
                                       outofpoints=c("stop", "warn", "ignore"),
                                       cluster=FALSE, 
                                       class=FALSE, ...) {
  if(rJava::is.jnull(x$javaObj)) stop("The java object in the DSD is null. This is probably due to a restart that has reset the JVM. Try recreating the DSD_RandomRBFSubspaceGeneratorEvents")
  res <- rJava::.jcall("DataStreamAccessor",returnSig="[[D",method="getPoints",x$javaObj,as.integer(n),simplify=T)
  res <- data.frame(res)
  #Ground truth is normally put into the last column
  ground_truth_classes <- res[,ncol(res)]
  if(cluster) {
    #Add the ground truth as an attribute
    attr(res,"cluster") <- ground_truth_classes
  }
  if(!class) {
    #Drop the column containing the class of the data point
    return(res[,1:(ncol(res)-1)])
  } else if (class) {
    #rename the last column to class
    names(res)[ncol(res)] <- "class"
    return(res)
  }
}