
#'Unwrap a clustered object
#'
#'\code{unwrap} is a utility function that reformats a CBData or CMData object so
#'that each row is one observation (instead of one or more clusters). A new
#'`ID' variable is added to indicate clusters. This form can be useful for
#'setting up the data for a different package.
#'
#'@aliases unwrap unwrap.CBData
#'@export
#'@param object a \code{\link{CBData}} object
#'@param \dots other potential arguments; not currently used
#'@return For \code{unwrap.CBData}: a data frame with one row for each cluster element (having a binary
#'outcome) with the following standardized column names
#'@return \item{Trt}{factor, the treatment group}
#'@return \item{ClusterSize}{numeric, the cluster size}
#'@return \item{ID}{factor, each level representing a different cluster}
#'@return \item{Resp}{numeric with 0/1 values, giving the response of the cluster
#'element}
#'@author Aniko Szabo
#'@keywords manip
#'@examples
#'
#'data(shelltox)
#'ush <- unwrap(shelltox)
#'head(ush)
#'

unwrap <- function(object,...) UseMethod("unwrap")
