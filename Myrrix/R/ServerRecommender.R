#' @title Object of class ServerRecommender
#' @description The core implementation of Recommender and furthermore MyrrixRecommender 
#' that lies inside the Serving Layer.\cr
#' This is the recommendation engine class.\cr
#' 
#' The ServerRecommender has a local model, allowing it to build recommendation models based on data which are locally
#' stored on your disk.\cr
#' Next to the local serving, it also allows to build recommendation models based on data which is distributed on 
#' Hadoop. Special out-of-the-box classes exists to let it run on CDH, AWS and on Hadoop clusters. If you run
#' the ServerRecommender in a distributed mode, we assume that you have set up the Computation layer already. 
#' This R package allows to ingest new data, update the model, get recommendations and similarities 
#' based on the recommendation engine which is running.
#' 
#'
#' @section Slots: 
#'  \describe{
#'    \item{\code{recommender}:}{A java object of class net.myrrix.online.ServerRecommender}
#'  }
#' @param bucket character string with the bucket that Serving Layer is using for instances
#' @param instanceID character string with the instance ID that the Serving Layer is serving. May be 0 for local mode.
#' @param localInputDir character string with the local input and model file directory
#' @param partition integer with the partition number in a partitioned distributed mode. 0 if not partitioned.
#' @param allPartitions reference to an object that can describe all partitions; only used to get their count 
#' (seee http://myrrix.com/docs/serving/javadoc/index.html)
#' @name ServerRecommender-class 
#' @rdname ServerRecommender-class
#' @aliases ServerRecommender-class
#' @exportClass ServerRecommender
#' @importFrom rJava .jnew
#' @usage Local setup: new("ServerRecommender", localInputDir)
#' @usage Distributed setup: new("ServerRecommender", bucket, instanceID, localInputDir, partition, allPartitions)
#' @examples
#' recommendationengine <- new("ServerRecommender", localInputDir=tempdir())
#' recommendationengine
setClass(Class="ServerRecommender",
         representation=representation(recommender="jobjRef"))
setMethod(f="initialize", signature="ServerRecommender",
          definition = function(.Object, localInputDir=NULL, bucket=NULL, instanceID = NULL, partition = NULL, allPartitions = NULL) {
            print(nargs())
            if(nargs() > 2){
              .NotYetImplemented()
            }else if(!is.null(localInputDir)){
              inputdir <- .jnew("java.io.File", localInputDir)
              .Object@recommender <- .jnew("net.myrrix.online.ServerRecommender", inputdir)
            }else{
              stop("need to supply either localInputDir or all elements localInputDir, bucket, instanceID, partition, allPartitions")
            }
            .Object
          })