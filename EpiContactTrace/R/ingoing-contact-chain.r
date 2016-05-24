## Copyright 2013 Stefan Widgren and Maria Noremark,
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

##' \code{IngoingContactChain}
##'
##' The ingoing contact chain is the number of holdings in the network
##' of direct and indirect contacts to the root holding, with regard
##' to temporal and order of the contacts during the defined time
##' window used for contact tracing.
##'
##'
##' The time period used for \code{IngoingContactChain} can either be
##' specified using \code{tEnd} and \code{days} or \code{inBegin} and
##' \code{inEnd}.
##'
##' If using \code{tEnd} and \code{days}, the time period for ingoing
##' contacts ends at \code{tEnd} and starts at \code{days} prior to
##' \code{tEnd}. The indegree will be calculated for each combination
##' of \code{root}, \code{tEnd} and \code{days}.
##'
##' An alternative way is to use \code{inBegin} and \code{inEnd}.  The
##' time period for ingoing contacts starts at inBegin and ends at
##' inEndDate. The vectors \code{root} \code{inBegin}, \code{inEnd}
##' must have the same lengths and the indegree will be calculated for
##' each index of them.
##'
##' The movements in \code{IngoingContactChain} is a \code{data.frame}
##' with the following columns: \describe{
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
##' @rdname IngoingContactChain-methods
##' @docType methods
##' @keywords methods
##' @include Contacts.r
##' @include ContactTrace.r
##' @seealso \code{\link{NetworkSummary}}
##' @param x a ContactTrace object, or a list of ContactTrace objects
##' or a \code{data.frame} with movements of animals between holdings,
##' see \code{\link{Trace}} for details.
##' @param ... Additional arguments to the method
##' @param root vector of roots to calculate ingoing contact chain
##' for.
##' @param tEnd the last date to include ingoing movements. Defaults
##' to \code{NULL}
##' @param days the number of previous days before tEnd to include
##' ingoing movements. Defaults to \code{NULL}
##' @param inBegin the first date to include ingoing
##' movements. Defaults to \code{NULL}
##' @param inEnd the last date to include ingoing movements. Defaults
##' to \code{NULL}
##' @return A \code{data.frame} with the following columns:
##' \describe{
##'   \item{root}{
##'     The root of the contact tracing
##'   }
##'
##'   \item{inBegin}{
##'     The first date to include ingoing movements
##'   }
##'
##'   \item{inEnd}{
##'     The last date to include ingoing movements
##'   }
##'
##'   \item{inDays}{
##'     The number of days in the interval inBegin to inEnd
##'   }
##'
##'   \item{ingoingContactChain}{
##'     The \code{\link{IngoingContactChain}} of the root within the time-interval
##'   }
##' }
##'
##' @section Methods:
##' \describe{
##'   \item{\code{signature(x = "ContactTrace")}}{
##'     Get the IngoingContactChain of a \code{ContactTrace} object.
##'   }
##'
##'   \item{\code{signature(x = "data.frame")}}{
##'     Get the IngoingContactChain for a data.frame with movements,
##'     see details and examples.
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
##' ## Calculate ingoing contact chain from a ContactTrace object
##' ic.1 <- IngoingContactChain(contactTrace)
##'
##' ## Calculate ingoing contact chain using tEnd and days
##' ic.2 <- IngoingContactChain(transfers,
##'                             root=2645,
##'                             tEnd='2005-10-31',
##'                             days=91)
##'
##' ## Check that the result is identical
##' identical(ic.1, ic.2)
##'
##' ## Calculate ingoing contact chain for all included herds
##' ## First extract all source and destination from the dataset
##' root <- sort(unique(c(transfers$source,
##'                       transfers$destination)))
##'
##' ## Calculate ingoing contact chain
##' result <- IngoingContactChain(transfers,
##'                               root=root,
##'                               tEnd='2005-10-31',
##'                               days=91)
##' }
setGeneric("IngoingContactChain",
           signature = "x",
           function(x, ...) standardGeneric("IngoingContactChain"))

##' @rdname IngoingContactChain-methods
##' @export
setMethod("IngoingContactChain",
          signature(x = "Contacts"),
          function (x)
      {
          if(!identical(x@direction, "in")) {
              stop("Unable to determine IngoingContactChain for outgoing contacts")
          }

          return(length(setdiff(x@source,x@root)))
      }
)

##' @rdname IngoingContactChain-methods
##' @export
setMethod("IngoingContactChain",
          signature(x = "ContactTrace"),
          function (x)
      {
          return(NetworkSummary(x)[, c("root",
                                       "inBegin",
                                       "inEnd",
                                       "inDays",
                                       "ingoingContactChain")])
      }
)

##' @rdname IngoingContactChain-methods
##' @export
setMethod("IngoingContactChain",
          signature(x = "data.frame"),
          function(x,
                   root,
                   tEnd = NULL,
                   days = NULL,
                   inBegin = NULL,
                   inEnd = NULL)
      {
          if(missing(root)) {
              stop("Missing parameters in call to IngoingContactChain")
          }

          if(all(is.null(tEnd), is.null(days))) {
              outBegin <- inBegin
              outEnd <- inBegin
          } else {
              outBegin <- NULL
              outEnd <- NULL
          }

          return(NetworkSummary(x,
                                root,
                                tEnd,
                                days,
                                inBegin,
                                inEnd,
                                outBegin,
                                outEnd)[, c("root",
                                            "inBegin",
                                            "inEnd",
                                            "inDays",
                                            "ingoingContactChain")])
      }
)
