## Copyright 2013-2015 Stefan Widgren and Maria Noremark,
## National Veterinary Institute, Sweden
##
## Licensed under the EUPL, Version 1.1 or - as soon they
## will be approved by the European Commission - subsequent
## versions of the EUPL (the "Licence");
## You may not use this work except in compliance with the
## Licence.
## You may obtain a copy of the Licence at:
##
## http://ec.europa.eu/idabc/eupl
##
## Unless required by applicable law or agreed to in
## writing, software distributed under the Licence is
## distributed on an "AS IS" basis,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
## express or implied.
## See the Licence for the specific language governing
## permissions and limitations under the Licence.

##' \code{OutDegree}
##'
##' The number of herds with direct movements of animals from the root
##' herd during the defined time window used for tracing
##'
##'
##' The time period used for \code{OutDegree} can either be specified
##' using \code{tEnd} and \code{days} or \code{outBegin} and
##' \code{outEnd}.
##'
##' If using \code{tEnd} and \code{days}, the time period for outgoing
##' contacts ends at \code{tEnd} and starts at \code{days} prior to
##' \code{tEnd}. The outdegree will be calculated for each combination
##' of \code{root}, \code{tEnd} and \code{days}.
##'
##' An alternative way is to use \code{outBegin} and \code{outEnd}.
##' The time period for outgoing contacts starts at outBegin and ends
##' at outEndDate. The vectors \code{root} \code{outBegin},
##' \code{outEnd} must have the same lengths and the outdegree will be
##' calculated for each index of them.
##'
##' The movements in \code{OutDegree} is a \code{data.frame}
##' with the following columns:
##' \describe{
##'
##'   \item{source}{
##'     an integer or character identifier of the source holding.
##'   }
##'
##'   \item{destination}{
##'     an integer or character identifier of the destination holding.
##'   }
##'
##'   \item{t}{
##'     the Date of the transfer
##'   }
##'
##'   \item{id}{
##'     an optional character vector with the identity of the animal.
##'   }
##'
##'   \item{n}{
##'     an optional numeric vector with the number of animals moved.
##'   }
##'
##'   \item{category}{
##'     an optional character or factor with category of the animal e.g. Cattle.
##'   }
##' }
##'
##' @rdname OutDegree-methods
##' @docType methods
##' @keywords methods
##' @include Contacts.r
##' @include ContactTrace.r
##' @seealso \code{\link{NetworkSummary}}
##' @param x a ContactTrace object, or a list of ContactTrace objects
##' or a \code{data.frame} with movements of animals between holdings,
##' see \code{\link{Trace}} for details.
##' @param ... Additional arguments to the method
##' @param root vector of roots to calculate outdegree for.
##' @param tEnd the last date to include outgoing movements. Defaults
##' to \code{NULL}
##' @param days the number of previous days before tEnd to include
##' outgoing movements. Defaults to \code{NULL}
##' @param outBegin the first date to include outgoing
##' movements. Defaults to \code{NULL}
##' @param outEnd the last date to include outgoing movements. Defaults
##' to \code{NULL}
##' @return A \code{data.frame} with the following columns:
##' \describe{
##'   \item{root}{
##'     The root of the contact tracing
##'   }
##'
##'   \item{outBegin}{
##'     The first date to include outgoing movements
##'   }
##'
##'   \item{outEnd}{
##'     The last date to include outgoing movements
##'   }
##'
##'   \item{outDays}{
##'     The number of days in the interval outBegin to outEnd
##'   }
##'
##'   \item{outDegree}{
##'     The \code{\link{OutDegree}} of the root within the time-interval
##'   }
##' }
##' @section Methods:
##' \describe{
##'   \item{\code{signature(x = "ContactTrace")}}{
##'     Get the OutDegree of a \code{ContactTrace} object.
##'   }
##'
##'   \item{\code{signature(x = "data.frame")}}{
##'     Get the OutDegree for a data.frame with movements, see details and examples.
##'   }
##' }
##' @references \itemize{
##'   \item Dube, C., et al., A review of network analysis terminology
##'     and its application to foot-and-mouth disease modelling and policy
##'     development. Transbound Emerg Dis 56 (2009) 73-85, doi:
##'     10.1111/j.1865-1682.2008.01064.x
##'
##'   \item Noremark, M., et al., Network analysis
##'     of cattle and pig movements in Sweden: Measures relevant for
##'     disease control and riskbased surveillance.  Preventive Veterinary
##'     Medicine 99 (2011) 78-90, doi: 10.1016/j.prevetmed.2010.12.009
##' }
##' @examples
##' \dontrun{
##'
##' ## Load data
##' data(transfers)
##'
##' ## Perform contact tracing using tEnd and days
##' contactTrace <- Trace(movements=transfers,
##'                       root=2645,
##'                       tEnd='2005-10-31',
##'                       days=91)
##'
##' ## Calculate outdegree from a ContactTrace object
##' od.1 <- OutDegree(contactTrace)
##'
##' ## Calculate outdegree using tEnd and days
##' od.2 <- OutDegree(transfers,
##'                   root=2645,
##'                   tEnd='2005-10-31',
##'                   days=91)
##'
##' ## Check that the result is identical
##' identical(od.1, od.2)
##'
##' ## Calculate outdegree for all included herds
##' ## First extract all source and destination from the dataset
##' root <- sort(unique(c(transfers$source,
##'                       transfers$destination)))
##'
##' ## Calculate outdegree
##' result <- OutDegree(transfers,
##'                     root=root,
##'                     tEnd='2005-10-31',
##'                     days=91)
##' }
setGeneric("OutDegree",
           signature = "x",
           function(x, ...) standardGeneric("OutDegree"))

##' @rdname OutDegree-methods
##' @export
setMethod("OutDegree",
          signature(x = "Contacts"),
          function (x)
      {
          if(!identical(x@direction, "out")) {
              stop("Unable to determine OutDegree for ingoing contacts")
          }

          return(length(unique(x@destination[x@source==x@root])))
      }
)

##' @rdname OutDegree-methods
##' @export
setMethod("OutDegree",
          signature(x = "ContactTrace"),
          function (x)
      {
          return(NetworkSummary(x)[, c("root",
                                       "outBegin",
                                       "outEnd",
                                       "outDays",
                                       "outDegree")])
      }
)

##' @rdname OutDegree-methods
##' @export
setMethod("OutDegree",
          signature(x = "data.frame"),
          function(x,
                   root,
                   tEnd = NULL,
                   days = NULL,
                   outBegin = NULL,
                   outEnd = NULL)
      {
          if(missing(root)) {
              stop("Missing parameters in call to OutDegree")
          }

          if(all(is.null(tEnd), is.null(days))) {
              inBegin <- outBegin
              inEnd <- outBegin
          } else {
              inBegin <- NULL
              inEnd <- NULL
          }

          return(NetworkSummary(x,
                                root,
                                tEnd,
                                days,
                                inBegin,
                                inEnd,
                                outBegin,
                                outEnd)[, c("root",
                                            "outBegin",
                                            "outEnd",
                                            "outDays",
                                            "outDegree")])
      }
)
