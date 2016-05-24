## Copyright 2013-2014 Stefan Widgren and Maria Noremark,
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

##' \code{ShortestPaths}
##'
##' Methods for function \code{ShortestPaths} in package \pkg{EpiContactTrace}
##' to get the shortest distance from/to the root given by the contact tracing.
##'
##' The contact tracing performs a depth first search starting at the root. The
##' \code{ShortestPaths} gives the shortest distance from root at each node.
##' The network tree structure given by the depth first search is shown by
##' \code{\link{show}}.
##'
##' @rdname ShortestPaths-methods
##' @docType methods
##' @keywords methods
##' @include ContactTrace.r
##' @param x a \code{\linkS4class{ContactTrace}} object, or a
##' \code{data.frame} with movements of animals between holdings, see
##' \code{\link{Trace}} for details.
##' @param ... Additional arguments to the method
##' @param root vector of roots to calculate shortest path for.
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
##'     If the direction is ingoing, then inBegin equals inBegin in
##'     \code{\link{Trace}} else NA.
##'   }
##'
##'   \item{inEnd}{
##'     If the direction is ingoing, then inEnd equals inEnd in
##'     \code{\link{Trace}} else NA.
##'   }
##'
##'   \item{outBegin}{
##'     If the direction is outgoing, then outBegin equals outBegin in
##'     \code{\link{Trace}} else NA.
##'   }
##'
##'   \item{outEnd}{
##'     If the direction is outgoing, then outEnd equals outEnd in
##'     \code{\link{Trace}} else NA.
##'   }
##'
##'   \item{direction}{
##'     If the direction is ingoing, then direction equals 'in' else 'out'
##'   }
##'
##'   \item{source}{
##'     The source of the contact at distance from root
##'   }
##'
##'   \item{destination}{
##'     The destination of the contact at distance from root
##'   }
##'
##'   \item{distance}{
##'     The shortest distance from/to root in the depth first search
##'   }
##' }
##' @section Methods: \describe{
##'
##'   \item{\code{signature(object = "ContactTrace")}}{
##'     Get the shortest paths for the ingoing and outgoing
##'     \code{Contacts} of a \code{ContactTrace} object.
##'   }
##'
##'   \item{\code{signature(x = "data.frame")}}{
##'     Get the shortest paths for a data.frame with movements,
##'     see details and examples.
##'   }
##' }
##' @seealso \code{\link{show}} and \code{\link{NetworkStructure}}.
##' @examples
##' \dontrun{
##'
##' ## Load data
##' data(transfers)
##'
##' ## Perform contact tracing
##' contactTrace <- Trace(movements=transfers,
##'                       root=2645,
##'                       tEnd='2005-10-31',
##'                       days=90)
##'
##' ShortestPaths(contactTrace)
##'
##' ## Calculate shortest paths for all included herds
##' ## First extract all source and destination from the dataset
##' root <- sort(unique(c(transfers$source, transfers$destination)))
##'
##' sp <- ShortestPaths(transfers, root=root, tEnd='2005-10-31', days=90)
##'
##' }
setGeneric('ShortestPaths',
           signature = 'x',
           function(x, ...) standardGeneric('ShortestPaths'))

##' @rdname ShortestPaths-methods
##' @export
setMethod('ShortestPaths',
          signature(x = 'ContactTrace'),
          function(x)
      {
          ns <- NetworkStructure(x)
          ns.in <- ns[ns$direction == "in",]
          ns.out <- ns[ns$direction == "out",]

          result <- NULL
          if(nrow(ns.in)) {
              ns.in <- ns.in[order(ns.in$distance, ns.in$source),]
              ns.in <- ns.in[!duplicated(ns.in$source),]
              ns.in$destination <- NA_character_
              result <- ns.in
          }

          if(nrow(ns.out)) {
              ns.out <- ns.out[order(ns.out$distance, ns.out$destination),]
              ns.out <- ns.out[!duplicated(ns.out$destination),]
              ns.out$source <- NA_character_
              result <- rbind(result, ns.out)
          }

          if(is.null(result)) {
              result <- ns
          } else {
              rownames(result) <- NULL
          }

          return(result)
      }
)

##' @rdname ShortestPaths-methods
##' @export
setMethod('ShortestPaths',
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
              stop('Missing root in call to ShortestPaths')
          }

          if(any(is.factor(root), is.integer(root))) {
              root <- as.character(root)
          } else if(is.numeric(root)) {
              ## root is supposed to be a character or integer
              ## identifier so test that root is a integer the same
              ## way as binom.test test x
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
                  stop('Use either tEnd and days or inBegin, inEnd, outBegin and outEnd in call to ShortestPaths')
              }

              if(any(is.character(tEnd), is.factor(tEnd))) {
                  tEnd <- as.Date(tEnd)
              }

              if(!identical(class(tEnd), 'Date')) {
                  stop("'tEnd' must be a Date vector")
              }

              ## Test that days is a nonnegative integer the same way
              ## as binom.test test x
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
                  stop('Use either tEnd and days or inBegin, inEnd, outBegin and outEnd in call to ShortestPaths')
              }
          } else {
              stop('Use either tEnd and days or inBegin, inEnd, outBegin and outEnd in call to ShortestPaths')
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

          sp <- .Call("shortestPaths",
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

          result <- NULL
          if(length(sp$inIndex)) {
              result <- data.frame(root        = root[sp$inIndex],
                                   inBegin     = inBegin[sp$inIndex],
                                   inEnd       = inEnd[sp$inIndex],
                                   outBegin    = as.Date(NA_character_),
                                   outEnd      = as.Date(NA_character_),
                                   direction   = 'in',
                                   source      = x$source[sp$inRowid],
                                   destination = NA_character_,
                                   distance    = sp$inDistance,
                                   stringsAsFactors=FALSE)
          }

          if(length(sp$outIndex)) {
              result <- rbind(result,
                              data.frame(root        = root[sp$outIndex],
                                         inBegin     = as.Date(NA_character_),
                                         inEnd       = as.Date(NA_character_),
                                         outBegin    = outBegin[sp$outIndex],
                                         outEnd      = outEnd[sp$outIndex],
                                         direction   = 'out',
                                         source      = NA_character_,
                                         destination = x$destination[sp$outRowid],
                                         distance    = sp$outDistance,
                                         stringsAsFactors=FALSE))
          }

          if(is.null(result)) {
              result <- data.frame(root        = character(0),
                                   inBegin     = as.Date(character(0)),
                                   inEnd       = as.Date(character(0)),
                                   outBegin    = as.Date(character(0)),
                                   outEnd      = as.Date(character(0)),
                                   direction   = character(0),
                                   source      = character(0),
                                   destination = character(0),
                                   distance    = integer(0),
                                   stringsAsFactors=FALSE)
          } else {
              rownames(result) <- NULL
          }

          return(result)
     }
)
