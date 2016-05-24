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

##' \code{NetworkSummary}
##'
##' \code{NetworkSummary} gives a summary of the contact tracing including the
##' time-window, \code{\link{InDegree}}, \code{\link{OutDegree}},
##' \code{\link{IngoingContactChain}} and \code{\link{OutgoingContactChain}}.
##'
##'
##' The time period used for \code{NetworkSummary} can either be specified
##' using \code{tEnd} and \code{days} or \code{inBegin}, \code{inEnd},
##' \code{outBegin} and \code{outEnd}.
##'
##' If using \code{tEnd} and \code{days}, the time period for ingoing
##' and outgoing contacts ends at \code{tEnd} and starts at
##' \code{days} prior to \code{tEnd}. The network summary will be
##' calculated for each combination of \code{root}, \code{tEnd} and
##' \code{days}.
##'
##' An alternative way is to use \code{inBegin}, \code{inEnd},
##' \code{outBegin} and \code{outEnd}. The time period for ingoing
##' contacts starts at inBegin and ends at inEndDate.  For outgoing
##' contacts the time period starts at outBegin and ends at outEnd.
##' The vectors \code{root} \code{inBegin}, \code{inEnd},
##' \code{outBegin} and \code{outEnd} must have the same lengths and
##' the network summary will be calculated for each index of them.
##'
##' The movements in \code{NetworkSummary} is a \code{data.frame}
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
##' @rdname NetworkSummary-methods
##' @docType methods
##' @param x a ContactTrace object or a \code{data.frame} with
##' movements of animals between holdings, see \code{\link{Trace}} for
##' details.
##' @param ... Additional arguments to the method
##' @param root vector of roots to calculate network summary for.
##' @param tEnd the last date to include ingoing movements. Defaults
##' to \code{NULL}
##' @param days the number of previous days before tEnd to include
##' ingoing movements. Defaults to \code{NULL}
##' @param inBegin the first date to include ingoing
##' movements. Defaults to \code{NULL}
##' @param inEnd the last date to include ingoing movements. Defaults
##' to \code{NULL}
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
##'   \item{inBegin}{
##'     Equals inBegin in \code{\link{Trace}}
##'   }
##'
##'   \item{inEnd}{
##'     Equals inEnd in \code{\link{Trace}}
##'   }
##'
##'   \item{outBegin}{
##'     Equals outBegin in \code{\link{Trace}}
##'   }
##'
##'   \item{outEnd}{
##'     Equals outEnd in \code{\link{Trace}}
##'   }
##'
##'   \item{inDegree}{
##'     The \code{\link{InDegree}} of the contact tracing
##'   }
##'
##'   \item{outDegree}{
##'     The \code{\link{OutDegree}} of the contact tracing
##'   }
##'
##'   \item{ingoingContactChain}{
##'     The \code{\link{IngoingContactChain}} of the contact tracing
##'   }
##'
##'   \item{outgoingContactChain}{
##'     The \code{\link{OutgoingContactChain}} of the contact tracing
##'   }
##' }
##'
##' @section Methods:
##' \describe{
##'   \item{\code{signature(x = "ContactTrace")}}{
##'     Get the network summary for the ingoing and outgoing
##'     \code{Contacts} of a ContactTrace object.
##'   }
##'
##'   \item{\code{signature(x = "data.frame")}}{
##'     Get the network summary for a data.frame with movements,
##'     see details and examples.
##'   }
##' }
##'
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
##' @keywords methods
##' @useDynLib EpiContactTrace
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
##' ## Calculate network summary from a ContactTrace object
##' ns.1 <- NetworkSummary(contactTrace)
##'
##' ## Calculate network summary using tEnd and days
##' ns.2 <- NetworkSummary(transfers,
##'                        root=2645,
##'                        tEnd='2005-10-31',
##'                        days=91)
##'
##' ## Check that the result is identical
##' identical(ns.1, ns.2)
##'
##' ## Calculate network summary using inBegin, inEnd
##' ## outBegin and outEnd
##' ns.3 <- NetworkSummary(transfers,
##'                        root=2645,
##'                        inBegin='2005-08-01',
##'                        inEnd='2005-10-31',
##'                        outBegin='2005-08-01',
##'                        outEnd='2005-10-31')
##'
##' ## Check that the result is identical
##' identical(ns.2, ns.3)
##'
##' ## When calculating the network summary for a data.frame of movements
##' ## a data.frame for each combination of root, tEnd and days are returned.
##' root <- c(1,2,3)
##' tEnd <- c("2005-09-01", "2005-10-01")
##' days <- c(30, 45)
##'
##' ## The network summary are calculated at the following
##' ## 12 combinations.
##' ## root = 1, tEnd = "2005-09-01", days = 30
##' ## root = 1, tEnd = "2005-09-01", days = 45
##' ## root = 1, tEnd = "2005-10-01", days = 30
##' ## root = 1, tEnd = "2005-10-01", days = 45
##' ## root = 2, tEnd = "2005-09-01", days = 30
##' ## root = 2, tEnd = "2005-09-01", days = 45
##' ## root = 2, tEnd = "2005-10-01", days = 30
##' ## root = 2, tEnd = "2005-10-01", days = 45
##' ## root = 3, tEnd = "2005-09-01", days = 30
##' ## root = 3, tEnd = "2005-09-01", days = 45
##' ## root = 3, tEnd = "2005-10-01", days = 30
##' ## root = 3, tEnd = "2005-10-01", days = 45
##' NetworkSummary(transfers, root, tEnd, days)
##'
##' ## Create a network summary for all included herds
##' ## First extract all source and destination from the dataset
##' root <- sort(unique(c(transfers$source,
##'                       transfers$destination)))
##'
##' ## Perform contact tracing using tEnd and days
##' result.1 <- NetworkSummary(transfers,
##'                            root=root,
##'                            tEnd='2005-10-31',
##'                            days=90)
##'
##' ## Perform contact tracing using inBegin, inEnd, outBegin and outEnd.
##' result.2 <- NetworkSummary(transfers,
##'                            root=root,
##'                            inBegin=rep('2005-08-02', length(root)),
##'                            inEnd=rep('2005-10-31', length(root)),
##'                            outBegin=rep('2005-08-02', length(root)),
##'                            outEnd=rep('2005-10-31', length(root)))
##' }
setGeneric('NetworkSummary',
           signature = 'x',
           function(x, ...) standardGeneric('NetworkSummary'))

##' @rdname NetworkSummary-methods
##' @export
setMethod('NetworkSummary',
          signature(x = 'ContactTrace'),
          function(x)
      {
          data.frame(root=x@root,
                     inBegin=x@ingoingContacts@tBegin,
                     inEnd=x@ingoingContacts@tEnd,
                     inDays=as.integer(x@ingoingContacts@tEnd - x@ingoingContacts@tBegin),
                     outBegin=x@outgoingContacts@tBegin,
                     outEnd=x@outgoingContacts@tEnd,
                     outDays=as.integer(x@outgoingContacts@tEnd - x@outgoingContacts@tBegin),
                     inDegree=InDegree(x@ingoingContacts),
                     outDegree=OutDegree(x@outgoingContacts),
                     ingoingContactChain=IngoingContactChain(x@ingoingContacts),
                     outgoingContactChain=OutgoingContactChain(x@outgoingContacts))
      }
)

##' @rdname NetworkSummary-methods
##' @export
setMethod('NetworkSummary',
          signature(x = 'data.frame'),
          function(x,
                   root,
                   tEnd = NULL,
                   days = NULL,
                   inBegin = NULL,
                   inEnd = NULL,
                   outBegin = NULL,
                   outEnd = NULL)
      {
          ## Check that arguments are ok from various perspectives...

          ## Check the data.frame x with movements
          if(!all(c('source', 'destination', 't') %in% names(x))) {
              stop('x must contain the columns source, destination and t.')
          }

          if(any(is.factor(x$source), is.integer(x$source))) {
              x$source <- as.character(x$source)
          } else if(!is.character(x$source)) {
              stop('invalid class of column source in x')
          }

          if(any(is.factor(x$destination), is.integer(x$destination))) {
              x$destination <- as.character(x$destination)
          } else if(!is.character(x$destination)) {
              stop('invalid class of column destination in x')
          }

          if(any(is.character(x$t), is.factor(x$t))) {
              x$t <- as.Date(x$t)
          }

          if(!identical(class(x$t), 'Date')) {
              stop('invalid class of column t in x')
          }

          if(any(is.na(x$t))) {
              stop("t in x contains NA")
          }

          ## Make sure the columns are in expected order and remove
          ## non-unique observations
          x <- unique(x[, c('source', 'destination', 't')])

          ## Check root
          if(missing(root)) {
              stop('Missing root in call to NetworkSummary')
          }

          if(any(is.factor(root), is.integer(root))) {
              root <- as.character(root)
          } else if(is.numeric(root)) {
              ## root is supposed to be a character or integer identifier
              ## so test that root is a integer the same way as binom.test test x
              rootr <- round(root)
              if(any(max(abs(root - rootr) > 1e-07))) {
                  stop("'root' must be an integer or character")
              }

              root <- as.character(rootr)
          } else if(!is.character(root)) {
              stop('invalid class of root')
          }

          ## Check if we are using the combination of tEnd and days or
          ## specify inBegin, inEnd, outBegin and outEnd
          if(all(!is.null(tEnd), !is.null(days))) {
              ## Using tEnd and days...check that
              ## inBegin, inEnd, outBegin and outEnd is NULL
              if(!all(is.null(inBegin), is.null(inEnd), is.null(outBegin), is.null(outEnd))) {
                  stop('Use either tEnd and days or inBegin, inEnd, outBegin and outEnd in call to NetworkSummary')
              }

              if(any(is.character(tEnd), is.factor(tEnd))) {
                  tEnd <- as.Date(tEnd)
              }

              if(!identical(class(tEnd), 'Date')) {
                  stop("'tEnd' must be a Date vector")
              }

              ## Test that days is a nonnegative integer the same way as binom.test test x
              daysr <- round(days)
              if (any(is.na(days) | (days < 0)) || max(abs(days - daysr)) > 1e-07) {
                  stop("'days' must be nonnegative and integer")
              }
              days <- daysr

              ## Make sure root, tEnd and days are unique
              root <- unique(root)
              tEnd <- unique(tEnd)
              days <- unique(days)

              n.root <- length(root)
              n.tEnd <- length(tEnd)
              n.days <- length(days)
              n <- n.root * n.tEnd * n.days

              root <- rep(root, each=n.tEnd*n.days, length.out=n)
              inEnd <- rep(tEnd, each=n.days, length.out=n)
              inBegin <- inEnd - rep(days, each=1, length.out=n)
              outEnd <- inEnd
              outBegin <- inBegin
          } else if(all(!is.null(inBegin), !is.null(inEnd), !is.null(outBegin), !is.null(outEnd))) {
              ## Using tEnd and days...check that
              ## Using inBegin, inEnd, outBegin and outEnd...check that
              ## tEnd and days are NULL
              if(!all(is.null(tEnd), is.null(days))) {
                  stop('Use either tEnd and days or inBegin, inEnd, outBegin and outEnd in call to NetworkSummary')
              }
          } else {
              stop('Use either tEnd and days or inBegin, inEnd, outBegin and outEnd in call to NetworkSummary')
          }

          ##
          ## Check inBegin
          ##
          if(any(is.character(inBegin), is.factor(inBegin))) {
              inBegin <- as.Date(inBegin)
          }

          if(!identical(class(inBegin), 'Date')) {
              stop("'inBegin' must be a Date vector")
          }

          if(any(is.na(inBegin))) {
              stop('inBegin contains NA')
          }

          ##
          ## Check inEnd
          ##
          if(any(is.character(inEnd), is.factor(inEnd))) {
              inEnd <- as.Date(inEnd)
          }

          if(!identical(class(inEnd), 'Date')) {
              stop("'inEnd' must be a Date vector")
          }

          if(any(is.na(inEnd))) {
              stop('inEnd contains NA')
          }

          ##
          ## Check outBegin
          ##
          if(any(is.character(outBegin), is.factor(outBegin))) {
              outBegin <- as.Date(outBegin)
          }

          if(!identical(class(outBegin), 'Date')) {
              stop("'outBegin' must be a Date vector")
          }

          if(any(is.na(outBegin))) {
              stop('outBegin contains NA')
          }

          ##
          ## Check outEnd
          ##
          if(any(is.character(outEnd), is.factor(outEnd))) {
              outEnd <- as.Date(outEnd)
          }

          if(!identical(class(outEnd), 'Date')) {
              stop("'outEnd' must be a Date vector")
          }

          if(any(is.na(outEnd))) {
              stop('outEnd contains NA')
          }

          ##
          ## Check ranges of dates
          ##
          if(any(inEnd < inBegin)) {
              stop('inEnd < inBegin')
          }

          if(any(outEnd < outBegin)) {
              stop('outEnd < outBegin')
          }

          ##
          ## Check length of vectors
          ##
          if(!identical(length(unique(c(length(root),
                                        length(inBegin),
                                        length(inEnd),
                                        length(outBegin),
                                        length(outEnd)))),
                        1L)) {
              stop('root, inBegin, inEnd, outBegin and outEnd must have equal length')
          }

          ## Arguments seems ok...go on with calculations

          ## Make sure all nodes have a valid variable name by making
          ## a factor of source and destination
          nodes <- as.factor(unique(c(x$source,
                                      x$destination,
                                      root)))

          ## Call networkSummary in EpiContactTrace.dll
          contact_chain<- .Call("networkSummary",
                                as.integer(factor(x$source, levels=levels(nodes))),
                                as.integer(factor(x$destination, levels=levels(nodes))),
                                as.integer(julian(x$t)),
                                as.integer(factor(root, levels=levels(nodes))),
                                as.integer(julian(inBegin)),
                                as.integer(julian(inEnd)),
                                as.integer(julian(outBegin)),
                                as.integer(julian(outEnd)),
                                length(nodes),
                                PACKAGE = "EpiContactTrace")

          return(data.frame(root=root,
                            inBegin=inBegin,
                            inEnd=inEnd,
                            inDays=as.integer(inEnd - inBegin),
                            outBegin=outBegin,
                            outEnd=outEnd,
                            outDays=as.integer(outEnd - outBegin),
                            inDegree=contact_chain[['inDegree']],
                            outDegree=contact_chain[['outDegree']],
                            ingoingContactChain=contact_chain[['ingoingContactChain']],
                            outgoingContactChain=contact_chain[['outgoingContactChain']]))
     }
)
