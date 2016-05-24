#' @title Methods to apply on objects of class ClientRecommender and ServerRecommender
#' @description Methods to apply on objects of class ClientRecommender and ServerRecommender.\cr
#' Objects of \code{ClientRecommender} and \code{ServerRecommender} behave similarly for the user. 
#' Both are classes which provide the interface with the recommendation engine of Myrrix, 
#' which is either running locally or in a distributed fashion.\cr
#' The methods which can be applied on this recommendation engine are \code{await}, \code{getAllItemIDs}, \code{getAllUserIDs}, 
#' \code{estimatePreference}, \code{mostPopularItems}, \code{recommend}.\cr
#' If Myrrix is running locally, you can set the hyperparameters of the recommendation engine which are set
#' in java system variables and are used by Myrrix. This can be done by using the provided methods
#' \code{setMyrrixHyperParameters} and \code{getMyrrixHyperParameters}. A full description of 
#' these hyperparamters which influence the model are listed below.
#' 
#' @section Methods: 
#'  \describe{
#'    \item{\code{setMyrrixHyperParameters(list)}:}{Set a list of hyperparameters for building and tuning the recommendation engine}
#'    \item{\code{getMyrrixHyperParameters()}:}{Get a list of hyperparameters which is currently used for building and tuning the recommendation engine}
#'    \item{\code{getMyrrixHyperParameters(parameters)}:}{Get a list of hyperparameters which is currently used for building and tuning the recommendation engine,
#'    limited to the parameters specified}
#'    \item{\code{await(ClientRecommender/ServerRecommender)}:}{Wait until the model is finished}
#'    \item{\code{getAllItemIDs(ClientRecommender/ServerRecommender)}:}{Get all item id's known to the model}
#'    \item{\code{getAllUserIDs(ClientRecommender/ServerRecommender)}:}{Get all user id's known to the model}
#'    \item{\code{estimatePreference(ClientRecommender/ServerRecommender, userID, itemIDs)}:}{Score a user for different items alongside the recommendation engine}
#'    \item{\code{mostPopularItems(ClientRecommender/ServerRecommender, howMany)}:}{Get the most popular items}
#'    \item{\code{recommend(ClientRecommender/ServerRecommender, userID, howMany)}:}{Recommend a number of items to a specific user}
#'  }
#' @section Hyperparameters: 
#'  \describe{
#'    \item{\code{model.iterations.max}:}{A hard limit of the number of iterations that will run in one build. Defaults to 30.}
#'    \item{\code{model.features}:}{Number of features to use when creating the matrix factorization. Defaults to 30.}
#'    \item{\code{model.als.iterations.convergenceThreshold}:}{Estimated strength values in the original matrix change a 
#'    little after each iteration, and less over time. If average absolute change in estimates is below this threshold, iteration will stop. 
#'    Defaults to 0.001.}
#'    \item{\code{model.als.lambda}:}{Controls the lambda overfitting parameter in the ALS algorithm. Defaults to 0.01}
#'    \item{\code{model.als.alpha}:}{Controls the alpha parameter in the ALS algorithm. Defaults to 40}
#'    \item{\code{model.noKnownItems}:}{If true, does not store in memory items that each user is already associated to. 
#'    This saves memory, but means that the recommender does not remember which items the user is already associated to. 
#'    These can't be automatically removed from consideration as recommendations. 
#'    This is desirable behavior in some contexts. To use this, the considerKnownItems argument to recommend must be true. 
#'    mostPopularItems will also not work with this flag enabled. Not recommended in general. Defaults to false.}
#'    \item{\code{model.local.writesBetweenRebuild}:}{Sets the number of new data points written to 
#'    the model that will trigger a model rebuild. Only applies to stand-alone mode. Defaults to 10000.}
#'    \item{\code{model.distributed.writesBetweenUpload}:}{Sets the number of new data points written to the model 
#'    that will trigger an upload of local data to the distributed storage system. Only applies to distributed mode. 
#'    Defaults to 50000}
#'    \item{\code{model.lsh.sampleRatio}:}{Enables locality-sensitive hashing to speed up the /recommend method, 
#'    at the cost of accuracy. Set to a value in (0,1]; LSH is disabled unless set to a value less than 1. 
#'    Recommended values are less than 0.1. This feature is experimental. Defaults to 1.0}  
#'  }
#'
#' @param object An object of class \code{ClientRecommender} or of class \code{ServerRecommender}
#' @param userID a user id for which to make the recommendation
#' @param itemIDs a vector of item id's for which to make the recommendation
#' @param params a list of hyperparameters to set for building the recommendation engine. Where
#' the names of the list elements need to be part of the specified hyperparameters below. See the examples.
#' @param parameters a character vector of names of hyperparameters to obtain the values. See the examples.
#' @param howMany an integer indicating how many popular items you want in the call to \code{mostPopularItems} and \code{recommend}
#' @param ... other arguments passed on to the methods
#' @name RecommenderMethods-methods
#' @rdname RecommenderMethods-methods
#' @aliases await await,ClientRecommender-method await,ServerRecommender-method
#' @seealso \code{\link{ClientRecommender-class}}, \code{\link{ServerRecommender-class}}
#' @exportMethod await
#' @importFrom rJava .jlong
#' @examples
#' ##
#' ## Set Hyperparameters to tune the Myrrix recommendation engine
#' ##
#' x <- getMyrrixHyperParameters()
#' str(x)
#' setMyrrixHyperParameters(
#'  params=list(model.iterations.max = 10, model.features=30, model.als.lambda=0.1))
#' x <- getMyrrixHyperParameters(
#'  parameters=c("model.iterations.max","model.features","model.als.lambda"))
#' str(x)
#' ##
#' ## Build a recommendation model locally
#' ##
#' \dontrun{
#' inputfile <- file.path(tempdir(), "audioscrobbler-data.subset.csv.gz")
#' download.file(
#'  url="http://dom2bevkhhre1.cloudfront.net/audioscrobbler-data.subset.csv.gz", 
#'  destfile = inputfile)
#' ## Set hyperparameters
#' setMyrrixHyperParameters(
#'  params=list(model.iterations.max = 2, model.features=10, model.als.lambda=0.1))
#' x <- getMyrrixHyperParameters(
#'  parameters=c("model.iterations.max","model.features","model.als.lambda"))
#' str(x)
#' ## Build a model which will be stored in getwd() and ingest the data file into it
#' recommendationengine <- new("ServerRecommender", localInputDir=getwd())
#' ingest(recommendationengine, inputfile)
#' await(recommendationengine)
#' ## Get all users/items and score
#' items <- getAllItemIDs(recommendationengine)
#' users <- getAllUserIDs(recommendationengine)
#' estimatePreference(recommendationengine, userID=users[5], itemIDs=items[1:20])
#' mostPopularItems(recommendationengine, howMany=10L)
#' recommend(recommendationengine, userID=users[5], howMany=10L)
#' }
setGeneric("await", function(object, ...) standardGeneric("await"))
setMethod("await", "ClientRecommender", function(object) object@recommender$await())
setMethod("await", "ServerRecommender", function(object) object@recommender$await())

#' @rdname RecommenderMethods-methods
#' @aliases getAllItemIDs getAllItemIDs,ServerRecommender-method getAllItemIDs,ClientRecommender-method
#' @exportMethod getAllItemIDs
setGeneric("getAllItemIDs", function(object, ...) standardGeneric("getAllItemIDs"))
setMethod("getAllItemIDs", "ClientRecommender", function(object) object@recommender$getAllItemIDs()$toArray())
setMethod("getAllItemIDs", "ServerRecommender", function(object) object@recommender$getAllItemIDs()$toArray())

#' @rdname RecommenderMethods-methods
#' @aliases getAllUserIDs getAllUserIDs,ServerRecommender-method getAllUserIDs,ClientRecommender-method
#' @exportMethod getAllUserIDs
setGeneric("getAllUserIDs", function(object, ...) standardGeneric("getAllUserIDs"))
setMethod("getAllUserIDs", "ClientRecommender", function(object) object@recommender$getAllUserIDs()$toArray())
setMethod("getAllUserIDs", "ServerRecommender", function(object) object@recommender$getAllUserIDs()$toArray())

#' @rdname RecommenderMethods-methods
#' @aliases estimatePreference estimatePreference,ServerRecommender,numeric,numeric-method estimatePreference,ClientRecommender,numeric,numeric-method 
#' @exportMethod estimatePreference
setGeneric("estimatePreference", function(object, userID, itemIDs, ...) standardGeneric("estimatePreference"))
.estimatePreference <- function(object, userID, itemIDs){
  stopifnot(length(userID) == 1)
  if(length(itemIDs) == 1){
    return(object@recommender$estimatePreference(.jlong(userID), .jlong(itemIDs)))
  }else{
    return(object@recommender$estimatePreferences(.jlong(userID), .jlong(itemIDs)))
  }
}
setMethod("estimatePreference", signature=signature(object = "ClientRecommender", userID="numeric", itemIDs="numeric"), definition = .estimatePreference)
setMethod("estimatePreference", signature=signature(object = "ServerRecommender", userID="numeric", itemIDs="numeric"), definition = .estimatePreference)

#' @rdname RecommenderMethods-methods
#' @aliases ingest ingest,ServerRecommender,character-method  ingest,ClientRecommender,character-method 
#' @exportMethod ingest
setGeneric("ingest", function(object, file, ...) standardGeneric("ingest"))
.ingest <- function(object, file){
  stopifnot(file.exists(file))
  ingestme <- .jnew("java.io.File", file)
  object@recommender$ingest(ingestme)
}
setMethod("ingest", signature=signature(object = "ClientRecommender", file="character"), definition = .ingest)
setMethod("ingest", signature=signature(object = "ServerRecommender", file="character"), definition = .ingest)

#' @rdname RecommenderMethods-methods
#' @aliases mostPopularItems mostPopularItems,ServerRecommender,integer-method  mostPopularItems,ClientRecommender,integer-method 
#' @exportMethod mostPopularItems
setGeneric("mostPopularItems", function(object, howMany, ...) standardGeneric("mostPopularItems"))
.mostPopularItems <- function(object, howMany){
  x <- object@recommender$mostPopularItems(howMany)
  scores <- .org.apache.mahout.cf.taste.recommender.RecommendedItem_list_to_Rlist(x)
  scores
}
setMethod("mostPopularItems", signature=signature(object = "ClientRecommender", howMany="integer"), definition = .mostPopularItems)
setMethod("mostPopularItems", signature=signature(object = "ServerRecommender", howMany="integer"), definition = .mostPopularItems)

#' @rdname RecommenderMethods-methods
#' @aliases recommend recommend,ServerRecommender,numeric,integer-method  recommend,ClientRecommender,numeric,integer-method 
#' @exportMethod recommend
setGeneric("recommend", function(object, userID, howMany, ...) standardGeneric("recommend"))
.recommend <- function(object, userID, howMany){
  x <- object@recommender$recommend(.jlong(userID), howMany)
  scores <- .org.apache.mahout.cf.taste.recommender.RecommendedItem_list_to_Rlist(x)
  scores
}
setMethod("recommend", signature=signature(object = "ClientRecommender", userID="numeric", howMany="integer"), definition = .recommend)
setMethod("recommend", signature=signature(object = "ServerRecommender", userID="numeric", howMany="integer"), definition = .recommend)
