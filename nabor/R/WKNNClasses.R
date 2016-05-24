#' Wrapper classes for k-NN searches enabling repeated queries of the same tree
#' 
#' @description \code{WKNNF} and \code{WKNND} are reference classes that wrap 
#'   C++ classes of the same name that include a space-efficient k-d tree along 
#'   with the target data points. They have \code{query} methods with exactly 
#'   the same interface as the \code{knn} function. One important point compared
#'   with \code{knn} - they must be intialised with floating point data and you 
#'   are responsible for this - see \code{\link{storage.mode}}) and the example 
#'   below.
#'   
#' @details \code{WKNNF} expects and returns matrices in R's standard (double, 8
#'   bytes) data type but uses floats internally. \code{WKNND} uses doubles 
#'   throughout. When retaining large numbers of points, the \code{WKNNF} 
#'   objects will have a small memory saving, especially if tree building is 
#'   delayed.
#'   
#'   The constructor for WKNN objects includes a logical flag indicating whether
#'   to build the tree immediately (default: \code{TRUE}) or (when \code{FALSE})
#'   to delay building the tree until a query is made (this happens 
#'   automatically when required).
#'   
#' @section Performance: The use of \code{WKNN} objects will incur a performance
#'   penalty for single queries of trees with < ~1000 data points. This is 
#'   because of the overhead associated with the R wrapper class. It therefore 
#'   makes sense to use \code{knn} in these circumstances.
#'   
#'   If you wish to make repeated queries of the same target data, then using 
#'   \code{WKNN} objects can give significant advantages. If you are going to 
#'   make repeated queries with the same set of query points (presumably against
#'   different target data), you can obtain benefits in some cases by converting
#'   the query points into WKNNF objects without building the trees.
#' @name WKNNF-class
#' @aliases WKNN
#' @rdname WKNN-class
#' @exportClass WKNNF
#' @export WKNNF
#' @seealso \code{\link{knn}}
#' @examples
#' ## Basic usage
#' # load sample data consisting of list of 3 separate 3d pointets
#' data(kcpoints)
#' # build a tree and query it with two different sets of points
#' w1 <- WKNNF(kcpoints[[1]])
#' w1q2 <- w1$query(kcpoints[[2]], k=5, eps=0)
#' str(w1q2)
#' w1q3 <- w1$query(kcpoints[[3]], k=5, eps=0)
#' # note that there will be small difference between WKNNF and knn due to loss 
#' # of precision in the double to float conversion when a WKNNF tree is 
#' # built and queried.
#' stopifnot(all.equal(knn(data=kcpoints[[1]], query=kcpoints[[2]], k=5, eps=0),
#'  w1q2, tolerance=1e-6))
#'  
#' ## storage mode: must be double
#' m=matrix(1:24, ncol=3)
#' storage.mode(m)
#' # this will generate an error unless we change to a double
#' w=tools::assertCondition(WKNND(m), "error")
#' storage.mode(m)="double"
#' w=WKNND(matrix(m, ncol=3))
#' 
#' ## construct wrapper objects but delay tree construction
#' w1 <- WKNNF(kcpoints[[1]], FALSE)
#' # query triggers tree construction
#' w1q2 <- w1$query(kcpoints[[2]], k=5, eps=0)
#' str(w1q2)
#' 
#' ## queries using wrapper objects
#' wkcpoints <-lapply(kcpoints, WKNNF, FALSE)
#' # query all 3 point sets against first
#' # this will trigger tree construction only for pointset 1
#' qall <- lapply(wkcpoints, function(x) wkcpoints[[1]]$queryWKNN(x$.CppObject, k=5, eps=0))
#' str(qall)
WKNNF <- setRcppClass("WKNNF")

#' @name WKNND-class
#' @rdname WKNN-class
#' @exportClass WKNND
#' @export WKNND
WKNND <- setRcppClass("WKNND")
