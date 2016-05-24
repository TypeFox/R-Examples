utils::globalVariables(".")
utils::globalVariables("RefDbcache")
.check.seedset<- function(object)
{
  if(!is(object, "seedset")) stop("object has to be of class \"seedset\" ")
  
  errors <- character()
  if (!is.igraph(object@GsMN))
    errors <- c(errors, "GsMN must be a igrpah object")  
  if(!is.list(object@seeds))
    errors <- c(errors, "seeds must be a list")
  if(length(errors) == 0){
    TRUE
  }else{
    errors
  }
}
#' \code{seedset-class}
#'
#' Object representing the seed sets of a given metabolic network
#'
#'@slot GsMN, a igraph network
#'@slot seeds, a character list represents seeds of a given metabolic 
#'      network which is composed of  the KEGG compound index.
#'      
#'@section method:
#'    \itemize{
#'      \item{getGsMN, \code{signature(object = "seedset")}:
#'        get the genome scale metabolic network whose seed set is caculated}
#'      \item{len, \code{signature(object = "seedset")}:
#'        return the number of source SCC}
#'      \item{seedSize, \code{signature(object = "seedset")}:
#'        returns the sizes of each source SCCs}
#'      \item{nonseed, \code{signature(object = "seedset")}:
#'        the non seeds of the GsMN}
#'      \item{show, \code{signature(object = "seedset")}:
#'        show the short summary of a seedset class}
#'      \item{confidencescore, \code{signature(object = "seedset")}:
#'        confidence score of the seed set}
#'  }
#' @name seedset-class
#' @rdname seedset-class
#' @exportClass seedset
#' @seealso \code{\link{getSeedSets}},\code{\link{getGsMN}},\code{\link{len}},
#' \code{\link{nonseed}},\code{\link{seedSize}},\code{\link{confidencescore}}
#' 
#' @examples 
#' \dontrun{
#' #' ## generate a metabolic network in igraph class and a seed set of this graph
#' annodir <- system.file("extdata","koanno.tab",package = "RevEcoR")
#' metabolic.data <- read.delim2(file=annodir,stringsAsFactors=FALSE)
#' g <- reconstructGsMN(metabolic.data)
#' seeds <- getSeedSets(g)@@seeds
#' seed.set <- new("seedset",GsMN = g, seeds = seeds)
#' }
 
setOldClass("igraph")
setClass("seedset", slot=list(GsMN="igraph",seeds="list"),
         prototype = list(GsMN=NULL,seeds=NULL),
         validity = .check.seedset)

#' Size of the each seed source component
#'
#' Caculate the size of each seed source component.
#'
#'@exportMethod seedSize
#'@rdname seedSize-methods
#'@name seedSize-methods 
#'@param object, \code{seedset} class
#'@aliases seedSize seedSize-methods 
#'@docType methods
#'@seealso \code{\link{seedset-class}}
#'@return a vector represents size of each source seed componet of network
#'@examples
#'\dontrun{
#'seedSize(seed.set)
#'}


setGeneric("seedSize",
  function(object){standardGeneric("seedSize")})
#' @rdname seedSize-methods
#' @aliases seedSize seedSize-methods
setMethod("seedSize",signature="seedset",
  function(object){
    ## listLen(object@seeds)
    lapply(object@seeds,length)
  }
)

#' The genome scale metabolic network 
#'
#' T he genome scale metabolic network (GsMN) whose seed set is caculated.
#'
#'@exportMethod getGsMN
#'@rdname getGsMN-methods
#'@param object, \code{seedset} class
#'@name getGsMN-methods
#'@aliases getGsMN getGsMN-methods 
#'@docType methods
#'@seealso \code{\link{seedset-class}}
#'@return a igraph
#'@examples
#'\dontrun{
#'getGsMN(seed.set)
#'}


setGeneric("getGsMN",
  function(object){standardGeneric("getGsMN")})
#' @rdname getGsMN-methods
#' @aliases getGsMN getGsMN-methods
setMethod("getGsMN",signature="seedset",
  function(object){
    object@GsMN
  }
)

#' Non seed of the network 
#'
#' Non seed of the network.
#'
#'@exportMethod nonseed
#'@rdname nonseed-methods
#'@param object, \code{seedset} class
#'@name nonseed-methods
#'@aliases nonseed nonseed-methods 
#'@docType methods
#'@seealso \code{\link{seedset-class}}
#'@return a vector
#'@examples
#'\dontrun{
#'nonseed(seed.set)
#'}

setGeneric("nonseed",
  function(object){standardGeneric("nonseed")})
#' @rdname nonseed-methods
#' @aliases nonseed nonseed-methods
setMethod("nonseed",signature="seedset",
  function(object){
    non.seed  <- V(object@GsMN)$name
    setdiff(non.seed,unlist(object@seeds))
  }
)


#' Conficence score 
#'
#' Caculate confidence score of seed set 
#'
#'@exportMethod confidencescore 
#'@rdname confidencescore-methods
#'@param object, \code{seedset} class
#'@name confidencescore-methods
#'@aliases confidencescore confidencescore-methods 
#'@docType methods
#'@seealso \code{\link{seedset-class}}
#'@return a list
#'@examples
#'\dontrun{
#'confidencescore(seed.set)
#'}


setGeneric("confidencescore",
  function(object){standardGeneric("confidencescore")})
#' @rdname confidencescore-methods
#' @aliases confidencescore confidencescore-methods
setMethod("confidencescore",signature="seedset",
  function(object){
    size <- unlist(seedSize(object)) 
    confidence.score  <- 1/size 
    map2(confidence.score,size, rep)
  }
)

#' the length of the seed set
#'
#' Caculate the number of the seed source components.
#'
#'@exportMethod len
#'@param object, \code{seed-set} class
#'@name len-methods
#'@rdname len-methods
#'@aliases len len-methods 
#'@docType methods
#'@seealso \code{\link{seedset-class}}
#'@return an interger
#' @examples
#'\dontrun{
#'len(seed.set)
#'}
setGeneric("len",
           function(object)
             standardGeneric("len")
)
#' @rdname len-methods
#' @aliases len len-methods
setMethod("len",valueClass = c("numeric"),signature="seedset",
          function(object){
            return(length(object@seeds))
          }
)

#' The show generic function
#'
#' Show a short summary of seedset object
#'
#'@exportMethod show 
#'@param object, \code{seed-set} class
#'@docType methods
#'@rdname show-methods
#'@aliases show show-methods
#'@seealso \code{\link{seedset-class}}
#'@examples
#'\dontrun{
#'show(seed.set)
#'}
setMethod("show",signature="seedset",
  function(object){
    direction  <- igraph::is.directed(object@GsMN)
    if (direction){
      direction  <- NULL
    }else{
      direction  <- "Directed Network"
    }
    nodes.no  <- length(V(object@GsMN))
    edges.no  <- length(E(object@GsMN))
    cat("Object of class ", class(object), "\n", sep = "")
    cat("  ", "IGRAPH:",direction,"--", nodes.no,"nodes",edges.no,"edges","--","\n")
    cat("  seedset length ",length(object@seeds), "\n")
  }
)

