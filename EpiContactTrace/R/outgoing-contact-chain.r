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

##' \code{OutgoingContactChain}
##'
##' The outgoing contact chain is the number of holdings in the network of
##' direct and indirect contacts from the root holding, with regard to temporal
##' and order of the contacts during the defined time window used for contact
##' tracing.
##'
##' @rdname OutgoingContactChain-methods
##' @docType methods
##' @keywords methods
##' @include Contacts.r
##' @include ContactTrace.r
##' @seealso \code{\link{NetworkSummary}}
##' @param x a ContactTrace object, or a list of ContactTrace objects
##' or a \code{data.frame} with movements of animals between holdings,
##' see \code{\link{Trace}} for details.
##' @param ... Additional arguments to the method
##' @param root vector of roots to calculate outgoing contact chain
##' for.
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
##'     The \code{\link{OutgoingContactChain}} of the root within the time-interval
##'   }
##' }
##' @section Methods:
##' \describe{
##'   \item{\code{signature(x = "ContactTrace")}}{
##'     Get the OutgoingContactChain of a \code{ContactTrace} object.
##'   }
##'
##'   \item{\code{signature(x = "data.frame")}}{
##'     Get the OutgoingContactChain for a data.frame with movements, see examples.
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
##' ## Calculate outgoing contact chain from a ContactTrace object
##' oc.1 <- OutgoingContactChain(contactTrace)
##'
##' ## Calculate outgoing contact chain using tEnd and days
##' oc.2 <- OutgoingContactChain(transfers,
##'                             root=2645,
##'                             tEnd='2005-10-31',
##'                             days=91)
##'
##' ## Check that the result is identical
##' identical(oc.1, oc.2)
##'
##' ## Calculate outgoing contact chain for all included herds
##' ## First extract all source and destination from the dataset
##' root <- sort(unique(c(transfers$source,
##'                       transfers$destination)))
##'
##' ## Calculate outgoing contact chain
##' result <- OutgoingContactChain(transfers,
##'                                root=root,
##'                                tEnd='2005-10-31',
##'                                days=91)
##' }
##'
setGeneric("OutgoingContactChain",
           signature = "x",
           function(x, ...) standardGeneric("OutgoingContactChain"))

##' @rdname OutgoingContactChain-methods
##' @export
setMethod("OutgoingContactChain",
          signature(x = "Contacts"),
          function (x)
      {
          if(!identical(x@direction, "out")) {
              stop("Unable to determine OutgoingContactChain for ingoing contacts")
          }

          return(length(setdiff(x@destination,x@root)))
      }
)

##' @rdname OutgoingContactChain-methods
##' @export
setMethod("OutgoingContactChain",
          signature(x = "ContactTrace"),
          function (x)
      {
          return(NetworkSummary(x)[, c("root",
                                       "outBegin",
                                       "outEnd",
                                       "outDays",
                                       "outgoingContactChain")])
      }
)

##' @rdname OutgoingContactChain-methods
##' @export
setMethod("OutgoingContactChain",
          signature(x = "data.frame"),
          function(x,
                   root,
                   tEnd = NULL,
                   days = NULL,
                   outBegin = NULL,
                   outEnd = NULL)
      {
          if(missing(root)) {
              stop("Missing parameters in call to OutgoingContactChain")
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
                                            "outgoingContactChain")])
      }
)
