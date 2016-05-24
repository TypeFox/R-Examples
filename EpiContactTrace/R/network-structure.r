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

##' \code{NetworkStructure}
##'
##' Methods for function \code{NetworkStructure} in package
##' \pkg{EpiContactTrace} to get the network tree structure from the contact
##' tracing.
##'
##' The contact tracing performs a depth first search starting at the root. The
##' \code{NetworkStructure} gives the distance from root at each node. The
##' network tree structure given by the depth first search is shown by
##' \code{\link{show}}.
##'
##' @rdname NetworkStructure-methods
##' @docType methods
##' @keywords methods
##' @include Contacts.r
##' @include ContactTrace.r
##' @param object A \code{\linkS4class{Contacts}} or
##' \code{linkS4class{ContactTrace}} object.
##' @return A \code{data.frame} with the following columns:
##' \describe{
##'   \item{root}{The root of the contact tracing}
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
##'     If the direction is outgoing, then outBegin equals
##'     outBegin in \code{\link{Trace}} else NA.
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
##'     The source of the contacts in the depth first search
##'   }
##'
##'   \item{destination}{
##'     The destination of the contacts in the depth first search
##'   }
##'
##'   \item{distance}{
##'     The distance from the destination to root in the depth first search
##'   }
##' }
##' @section Methods:
##' \describe{
##'   \item{\code{signature(object = "Contacts")}}{
##'     Get the network structure for the Contacts object.
##'   }
##'
##'   \item{\code{signature(object = "ContactTrace")}}{
##'     Get the network structure for the ingoing and outgoing
##'     \code{Contacts} of a \code{ContactTrace} object.
##'   }
##'
##'   \item{\code{signature(object = "list")}}{
##'     Get the network structure for a list of \code{ContactTrace}
##'     objects. Each item in the list must be a \code{ContactTrace} object.
##'   }
##' }
##' @seealso \code{\link{show}}.
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
##' NetworkStructure(contactTrace)
##' }
setGeneric("NetworkStructure",
           signature = "object",
           function(object) standardGeneric("NetworkStructure"))

##' @rdname NetworkStructure-methods
##' @export
setMethod("NetworkStructure",
          signature(object = "Contacts"),
          function(object)
      {
          if(length(object@source) > 0L) {
              ## Create a matrix with source, destination and distance
              m <- cbind(object@source[object@index],
                         object@destination[object@index],
                         object@distance,
                         deparse.level=0)

              ## To be able to identify duplicate rows, create strings from rows
              tmp <- apply(m, 1, function(x) paste(x, collapse="\r"))

              ## Identify which rows are not identical to previous rows.
              ## row[i] != row[i-1] for all i > 1
              i <- tmp[seq_len(length(tmp)-1)] != tmp[seq_len(length(tmp))[-1]]

              ## Select the i rows, including first row
              m <- as.data.frame(m[c(TRUE, i), , drop=FALSE], stringsAsFactors=FALSE)
              names(m) <- c("source", "destination", "distance")

              ## Convert distance from character to integer
              m$distance <- as.integer(m$distance)

              if(identical(object@direction, "in")) {
                  result <- data.frame(root=object@root,
                                       inBegin=object@tBegin,
                                       inEnd=object@tEnd,
                                       outBegin=as.Date(as.character(NA)),
                                       outEnd=as.Date(as.character(NA)),
                                       direction="in",
                                       stringsAsFactors=FALSE)
              } else {
                  result <- data.frame(root=object@root,
                                       inBegin=as.Date(as.character(NA)),
                                       inEnd=as.Date(as.character(NA)),
                                       outBegin=object@tBegin,
                                       outEnd=object@tEnd,
                                       direction="out",
                                       stringsAsFactors=FALSE)
              }

              return(cbind(result, m))
          } else {
              ## No contacts, return a zero row data.frame
              return(data.frame(root=character(0),
                                inBegin=as.Date(character(0)),
                                inEnd=as.Date(character(0)),
                                outBegin=as.Date(character(0)),
                                outEnd=as.Date(character(0)),
                                direction=character(0),
                                source=character(0),
                                destination=character(0),
                                distance=integer(0),
                                stringsAsFactors=FALSE))
          }
      }
)

##' @rdname NetworkStructure-methods
##' @export
setMethod("NetworkStructure",
          signature(object="ContactTrace"),
          function(object)
      {
          return(rbind(NetworkStructure(object@ingoingContacts),
                       NetworkStructure(object@outgoingContacts)))
      }
)
