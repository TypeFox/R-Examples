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

##' Trace Contacts.
##'
##' Contact tracing for a specied node(s) (root) during a specfied time period.
##' The time period is divided into two parts, one for ingoing contacts and one
##' for outgoing contacts.
##'
##'
##' The time period used for \code{Trace} can either be specified
##' using \code{tEnd} and \code{days} or \code{inBegin}, \code{inEnd},
##' \code{outBegin} and \code{outEnd}.
##'
##' If using \code{tEnd} and \code{days}, the time period for ingoing
##' and outgoing contacts ends at \code{tEnd} and starts at
##' \code{days} prior to \code{tEnd}. The tracing will be performed
##' for each combination of \code{root}, \code{tEnd} and \code{days}.
##'
##' An alternative way is to use \code{inBegin}, \code{inEnd},
##' \code{outBegin} and \code{outEnd}. The time period for ingoing
##' contacts starts at inBegin and ends at inEndDate.  For outgoing
##' contacts the time period starts at outBegin and ends at outEnd.
##' The vectors \code{root} \code{inBegin}, \code{inEnd},
##' \code{outBegin} and \code{outEnd} must have the same lengths and
##' the tracing will be performed for each index of them.
##'
##' The argument movements in Trace is a \code{data.frame}
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
##' @param movements a \code{data.frame} data.frame with movements,
##' see details.
##' @param root vector of roots to perform contact tracing for.
##' @param tEnd the last date to include ingoing and outgoing
##' movements. Defaults to \code{NULL}
##' @param days the number of previous days before tEnd to include
##' ingoing and outgoing movements. Defaults to \code{NULL}
##' @param inBegin the first date to include ingoing
##' movements. Defaults to \code{NULL}
##' @param inEnd the last date to include ingoing movements. Defaults
##' to \code{NULL}
##' @param outBegin the first date to include outgoing
##' movements. Defaults to \code{NULL}
##' @param outEnd the last date to include outgoing
##' movements. Defaults to \code{NULL}
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
##' @export
##' @examples
##' \dontrun{
##'
##' ## Load data
##' data(transfers)
##'
##' ## Perform contact tracing using tEnd and days
##' trace.1 <- Trace(movements=transfers,
##'                  root=2645,
##'                  tEnd='2005-10-31',
##'                  days=91)
##'
##' ## Perform contact tracing using inBegin, inEnd
##' ## outBegin and outEnd
##' trace.2 <- Trace(movements=transfers,
##'                  root=2645,
##'                  inBegin='2005-08-01',
##'                  inEnd='2005-10-31',
##'                  outBegin='2005-08-01',
##'                  outEnd='2005-10-31')
##'
##' ## Check that the result is identical
##' identical(trace.1, trace.2)
##'
##' ## Show result of contact tracing
##' show(trace.1)
##'
##' ## Create a network summary for all included herds
##' ## First extract all source and destination from the dataset
##' root <- sort(unique(c(transfers$source,
##'                       transfers$destination)))
##'
##' ## Perform contact tracing using tEnd and days.
##' trace.3 <- Trace(movements=transfers,
##'                  root=root,
##'                  tEnd='2005-10-31',
##'                  days=91)
##'
##' ## Perform contact tracing using inBegin, inEnd
##' ## outBegin and outEnd
##' trace.4 <- Trace(movements=transfers,
##'                  root=root,
##'                  inBegin=rep('2005-08-01', length(root)),
##'                  inEnd=rep('2005-10-31', length(root)),
##'                  outBegin=rep('2005-08-01', length(root)),
##'                  outEnd=rep('2005-10-31', length(root)))
##'
##' ## Check that the result is identical
##' identical(trace.3, trace.4)
##'
##' NetworkSummary(trace.3)
##' }
##'
Trace <- function(movements,
                  root,
                  tEnd = NULL,
                  days = NULL,
                  inBegin = NULL,
                  inEnd = NULL,
                  outBegin = NULL,
                  outEnd = NULL)
{
    ## Before doing any contact tracing check that arguments are ok
    ## from various perspectives.
    if(any(missing(movements),
           missing(root))) {
        stop('Missing parameters in call to Trace')
    }

    if(!is.data.frame(movements)) {
        stop('movements must be a data.frame')
    }

    if(!all(c('source', 'destination', 't') %in% names(movements))) {
        stop('movements must contain the columns source, destination and t.')
    }

    ##
    ## Check movements$source
    ##
    if(any(is.factor(movements$source), is.integer(movements$source))) {
        movements$source <- as.character(movements$source)
    } else if(!is.character(movements$source)) {
        stop('invalid class of column source in movements')
    }

    if(any(is.na(movements$source))) {
        stop('source in movements contains NA')
    }

    ##
    ## Check movements$destination
    ##
    if(any(is.factor(movements$destination), is.integer(movements$destination))) {
        movements$destination <- as.character(movements$destination)
    } else if(!is.character(movements$destination)) {
        stop('invalid class of column destination in movements')
    }

    if(any(is.na(movements$destination))) {
        stop('destination in movements contains NA')
    }

    ##
    ## Check movements$t
    ##
    if(any(is.character(movements$t), is.factor(movements$t))) {
        movements$t <- as.Date(movements$t)
    }
    if(!identical(class(movements$t), 'Date')) {
        stop('invalid class of column t in movements')
    }

    if(any(is.na(movements$t))) {
        stop('t in movements contains NA')
    }

    if('n' %in% names(movements)) {
        if(is.integer(movements$n)) {
            movements$n <- as.numeric(movements$n)
        } else if(!is.numeric(movements$n)) {
            stop('invalid class of column n in movements')
        }
    } else {
        movements$n <- as.numeric(NA)
    }

    if('id' %in% names(movements)) {
        if(any(is.factor(movements$id), is.integer(movements$id))) {
            movements$id <- as.character(movements$id)
        } else if(!is.character(movements$id)) {
            stop('invalid class of column id in movements')
        }
    } else {
        movements$id <- as.character(NA)
    }

    if('category' %in% names(movements)) {
        if(any(is.factor(movements$category), is.integer(movements$category))) {
            movements$category <- as.character(movements$category)
        } else if(!is.character(movements$category)) {
            stop('invalid class of column category in movements')
        }
    } else {
        movements$category <- as.character(NA)
    }

    ## Make sure the columns are in expected order
    if(!identical(names(movements), c('source',
                                      'destination',
                                      't',
                                      'id',
                                      'n',
                                      'category'))) {
        movements <- movements[, c('source',
                                   'destination',
                                   't',
                                   'id',
                                   'n',
                                   'category')]
    }

    ## Make sure that no duplicate movements exists
    movements <- unique(movements)

    ##
    ## Check root
    ##
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

    if(any(is.na(root))) {
        stop('root contains NA')
    }

    ## Check if we are using the combination of tEnd and days or
    ## specify inBegin, inEnd, outBegin and outEnd
    if(all(!is.null(tEnd), !is.null(days))) {
        ## Using tEnd and days...check that
        ## inBegin, inEnd, outBegin and outEnd is NULL
        if(!all(is.null(inBegin), is.null(inEnd), is.null(outBegin), is.null(outEnd))) {
            stop('Use either tEnd and days or inBegin, inEnd, outBegin and outEnd in call to Trace')
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
            stop('Use either tEnd and days or inBegin, inEnd, outBegin and outEnd in call to Trace')
        }
    } else {
        stop('Use either tEnd and days or inBegin, inEnd, outBegin and outEnd in call to Trace')
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

    ## Arguments seems ok...go on with contact tracing

    ## Make sure all nodes have a valid variable name by making
    ## a factor of source and destination
    nodes <- as.factor(unique(c(movements$source,
                                movements$destination,
                                root)))


    trace_contacts <- .Call("traceContacts",
                            as.integer(factor(movements$source, levels=levels(nodes))),
                            as.integer(factor(movements$destination, levels=levels(nodes))),
                            as.integer(julian(movements$t)),
                            as.integer(factor(root, levels=levels(nodes))),
                            as.integer(julian(inBegin)),
                            as.integer(julian(inEnd)),
                            as.integer(julian(outBegin)),
                            as.integer(julian(outEnd)),
                            length(nodes),
                            PACKAGE = "EpiContactTrace")

    result <- lapply(seq_len(length(root)), function(i) {
        j <- (i-1) * 4

        ## Extract data from contact tracing
        contacts_all <- movements[trace_contacts[[j + 1]], ]
        distance <- trace_contacts[[j + 2]]

        ## Since the algorithm might visit the same node more than once
        ## make sure we have unique contacts
        contacts <- unique(contacts_all)

        ## Create an index to contacts, so that the result matrix can be reconstructed
        ## from the contacts, combined with index and distance
        ## contacts_all <- cbind(contacts[index,], distance)
        index <- match(apply(contacts_all, 1, function(x) paste(x, collapse="\r")),
                       apply(contacts, 1, function(x) paste(x, collapse="\r")))

        ingoingContacts <- new('Contacts',
                               root = root[i],
                               tBegin = inBegin[i],
                               tEnd = inEnd[i],
                               source = contacts[,1],
                               destination = contacts[,2],
                               t = contacts[,3],
                               id = contacts[,4],
                               n = contacts[,5],
                               category = contacts[,6],
                               index = index,
                               distance = distance,
                               direction = 'in')

        ## Extract data from contact tracing
        contacts_all <- movements[trace_contacts[[j + 3]], ]
        distance <- trace_contacts[[j + 4]]

        ## Since the algorithm might visit the same node more than once
        ## make sure we have unique contacts
        contacts <- unique(contacts_all)

        ## Create an index to contacts, so that the result matrix can be reconstructed
        ## from the contacts, combined with index and distance
        ## contacts_all <- cbind(contacts[index,], distance)
        index <- match(apply(contacts_all, 1, function(x) paste(x, collapse="\r")),
                       apply(contacts, 1, function(x) paste(x, collapse="\r")))

        outgoingContacts <- new('Contacts',
                                root = root[i],
                                tBegin = outBegin[i],
                                tEnd = outEnd[i],
                                source = contacts[,1],
                                destination = contacts[,2],
                                t = contacts[,3],
                                id = contacts[,4],
                                n = contacts[,5],
                                category = contacts[,6],
                                index = index,
                                distance = distance,
                                direction = 'out')

        return(new('ContactTrace',
                   root = root[i],
                   ingoingContacts = ingoingContacts,
                   outgoingContacts = outgoingContacts))
    })

    ## Name each list item with ContactTrace objects to the name of the ContactTrace root.
    names(result) <-  sapply(result, function(listItem) listItem@ingoingContacts@root)

    if(identical(length(result), 1L))
      return(result[[1]])

    return(result)
}
