# Data documentation for Roxygen

#' Bell Labs 1router data from Cao et al. (2000)
#'
#' Data from 4-node network with star topology collected from Bell Labs; used in
#' Cao et al. (2000).
#'
#' @section Objects:
#' The list bell.labs, which contains several objects:
#' \itemize{
#'  \item \code{A}, the routing matrix for this network (truncated for full row
#'      rank)
#'  \item \code{df}, a data.frame with all data
#'  \item \code{X}, a matrix of origin-destination flows formatted for analysis
#'  \item \code{Y}, a matrix of link loads formatted for analysis
#'  \item \code{tvec}, a vector of times
#' }
#' In this data, we have \code{A \%*\% t(X) == t(Y)}.
#'
#' @section Variables:
#' The list bell.labs contains the following:
#' \itemize{
#'  \item The routing matrix \code{A}. The columns of this matrix correspond to
#'  individual OD flows (the columns of X), and its rows correspond to individual
#'  link loads (the columns of Y).
#'  \item The data.frame \code{df}, containing
#'  \itemize{
#'   \item value, level of traffic recorded
#'   \item nme, name of flow or load
#'   \item method, whether flow was directly observered or inferred
#'       (all observed)
#'   \item time, time of observation
#'   \item od, flag for origin-destination vs. link loads
#'   \item orig, origin of flow or load
#'   \item dest, destination of flow or load
#'   \item node, node involved in flow or load
#'  }
#'  \item The OD matrix X. Columns correspond to individual OD flows, and the rows
#'  correspond to observations.
#'  \item The link load matrix Y. Columns of the Y matrix correspond to individual
#'  link loads, and the rows correspond to observations.
#'  \item The vector tvec, containing the time in decimal hours since midnight for
#'  each observation. 
#' }
#'
#' @docType data
#' @name bell.labs
#' @usage bell.labs
#' @references J. Cao, D. Davis, S. Van Der Viel, and B. Yu.
#' Time-varying network tomography: router link data.
#' Journal of the American Statistical Association, 95:1063-75, 2000.
#' @keywords datasets
#' @family bell.labs
NULL

#' CMU data from Blocker & Airoldi (2011)
#'
#' Data from the 12 node CMU network used in Blocker & Airoldi (2011). The OD
#' flows are actual, observed traffic from a CMU network. The topology does not,
#' however, correspond to the original network due to security considerations.
#'
#' @section Objects:
#' The list cmu, which contains several objects:
#' \itemize{
#'  \item \code{A}, the routing matrix for this network (truncated for full row
#'      rank)
#'  \item \code{X}, a matrix of origin-destination flows formatted for analysis
#'  \item \code{Y}, a matrix of link loads formatted for analysis
#'  \item \code{A.full}, the routing matrix for this network without
#'  truncatation for full row rank)
#'  \item \code{Y.full}, a matrix of link loads corresponding to code{A.full}
#' }
#' In this data, we have \code{A \%*\% t(X) == t(Y)} and
#' \code{A.full \%*\% t(X) == t(Y.full)}
#'
#' @section Variables:
#' The list cmu contains the following:
#' \itemize{
#'  \item The routing matrix \code{A}. The columns of this matrix correspond to
#'  individual OD flows (the columns of X), and its rows correspond to individual
#'  link loads (the columns of Y).
#'  \item The OD matrix X. Columns correspond to individual OD flows, and the rows
#'  correspond to observations.
#'  \item The link load matrix Y. Columns of the Y matrix correspond to individual
#'  link loads, and the rows correspond to observations.
#'  \item The routing matrix \code{A.full}. This is the complete routing matrix
#'  before reduction for full row-rank.
#'  \item The link load matrix Y.full, corresponding to A.full.
#' }
#'
#' @docType data
#' @name cmu
#' @usage cmu
#' @references A.W. Blocker and E.M. Airoldi. Deconvolution of mixing
#' time series on a graph. Proceedings of the Twenty-Seventh Conference Annual
#' Conference on Uncertainty in Artificial Intelligence (UAI-11) 51-60, 2011.
#' @keywords datasets
#' @family cmu
NULL

#' Abilene data from Fang et al. (2007)
#'
#' Data from the 12 node Abilene network from Fang et al. (2007). Both the OD
#' flows and the topology correspond to the actual network. This is the X1
#' dataset from the given paper.
#'
#' @section Objects:
#' The list abilene, which contains several objects:
#' \itemize{
#'  \item \code{A}, the routing matrix for this network (truncated for full row
#'      rank)
#'  \item \code{X}, a matrix of origin-destination flows formatted for analysis
#'  \item \code{Y}, a matrix of link loads formatted for analysis
#'  \item \code{A.full}, the routing matrix for this network without
#'  truncatation for full row rank)
#'  \item \code{Y.full}, a matrix of link loads corresponding to code{A.full}
#' }
#' In this data, we have \code{A \%*\% t(X) == t(Y)} and
#' \code{A.full \%*\% t(X) == t(Y.full)}
#'
#' @section Variables:
#' The list abilene contains the following:
#' \itemize{
#'  \item The routing matrix \code{A}. The columns of this matrix correspond to
#'  individual OD flows (the columns of X), and its rows correspond to individual
#'  link loads (the columns of Y).
#'  \item The OD matrix X. Columns correspond to individual OD flows, and the rows
#'  correspond to observations.
#'  \item The link load matrix Y. Columns of the Y matrix correspond to individual
#'  link loads, and the rows correspond to observations.
#'  \item The routing matrix \code{A.full}. This is the complete routing matrix
#'  before reduction for full row-rank.
#'  \item The link load matrix Y.full, corresponding to A.full.
#' }
#'
#' @docType data
#' @name abilene
#' @usage abilene
#' @references J. Fang, Y. Vardi, and C.-H. Zhang. An iterative tomogravity
#' algorithm for the estimation of network traffic. In R. Liu, W. Strawderman,
#' and C.-H. Zhang, editors, Complex Datasets and Inverse Problems: Tomography,
#' Networks and Beyond, volume 54 of Lecture Notes-Monograph Series. IMS, 2007.
#' @keywords datasets
#' @family abilene
NULL
