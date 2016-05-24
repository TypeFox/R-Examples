###############################################################################
# R (http://r-project.org/) Instrument Class Model
#
# Copyright (c) 2009-2012
# Peter Carl, Dirk Eddelbuettel, Jeffrey Ryan, 
# Joshua Ulrich, Brian G. Peterson, and Garrett See
#
# This library is distributed under the terms of the GNU Public License (GPL)
# for full details see the file COPYING
#
# $Id: synthetic.R 1498 2013-08-25 00:26:39Z gsee $
#
###############################################################################

#' @export
#' @rdname synthetic.instrument
synthetic <- function(primary_id=NULL, currency=NULL, multiplier=1, 
                      identifiers=NULL, assign_i=TRUE, overwrite=TRUE, ..., 
                      members=NULL, type="synthetic")
{
    if (missing(primary_id) || (is.null(primary_id))) primary_id <- make_spread_id(members)
    if (missing(currency) || (is.null(currency))) {
        if (is.null(members)) {
            stop("'currency' is a required argument") 
        } else {
            instr <- try(getInstrument(members[[1]],silent=TRUE))
            if (is.instrument(instr)) currency <- instr$currency
        }
    }
    instrument(primary_id=primary_id , currency=currency , 
               multiplier=multiplier , identifiers = identifiers, 
               assign_i=assign_i, overwrite=overwrite, ...=..., type=type, 
               members=members)
}


#' synthetic instrument constructors
#' 
#' define spreads, guaranteed_spreads, butterflies, and other synthetic instruments
#'
#' Simple derivatives like \code{\link{option}} or \code{\link{future}} contracts typically have one underlying instrument.  
#' While properties like strike and expiration vary for these derivative contracts or series, the underlying is well understood.
#'
#' More complex derivatives are typically modeled as baskets of underlying products, and are typically traded over-the-counter or as proprietary in-house products.
#'
#' The general \code{synthetic} function is intended to be extended to support these arbitrary baskets of assets.  
#' 
#' \code{spread} \code{guaranteed_spread} and \code{butterfly} are wrappers for \code{synthetic.instrument}. \code{synthetic.instrument} will make a call to synthetic to create the final instrument.
#' 
#' The \code{suffix_id} parameter of wrapper functions such as  \code{guaranteed_spread} is presumed to 
#' be a string describing the \code{members}. 
#' It will be \code{\link{strsplit}} using the regex "[-;:_,\\.]" to create the \code{members} vector,
#' and potentially combined with a \code{root_id}.
#'
#' Most wrappers will build \code{primary_id} if it is NULL, either by combining \code{root_id} and \code{suffix_id}, or
#' by passing \code{members} in a call to \code{\link{make_spread_id}}
#'
#' \code{ICS} will build an Intercommodity Spread.  Although the expiration date and ratio may change, 
#' the members of a given ICS will not change.  Therefore, \code{ICS_root} can be used to hold the 
#' members of an Intercommodity Spread.  If an \code{ICS_root} has not been defined, then \code{members}
#' will be a required argument for \code{ICS}
#'
#' We welcome assistance from others to model more complex OTC derivatives such as swap products.
#'
#' @aliases synthetic.instrument synthetic spread guaranteed_spread butterfly
#' @param primary_id chr string of primary identifier of instrument to be defined.
#' @param currency chr string name of currency denomination
#' @param members vector of primary_ids of member instruments
#' @param memberratio vector of weights for each leg. negative numbers for selling.
#' @param \dots any other passthrough parameters
#' @param multiplier multiplier of the spread (1 / divisor for price weighted baskets)
#' @param tick_size minimum price change of the spread
#' @param identifiers identifiers
#' @param assign_i TRUE/FALSE. Should the instrument be assigned in the \code{.instrument} environment?
#' @param type type of instrument; wrappers do not require this.
#' @param root_id instrument identifier for the root contract, default NULL
#' @param suffix_id identifiers for the member contract suffixes, default NULL, 
#'   will be split as \code{members}, see Details
#' @param overwrite if FALSE and an instrument with the same \code{primary_id}
#'   is already defined, an error will be thrown and no instruments will be 
#'   created.
#' @return called for side effect. stores an instrument in .instrument environment
#' @author Brian Peterson, Garrett See
#' @seealso instrument, future, option_series.yahoo
#' @examples
#'
#' \dontrun{
#' stock('SPY','USD',1)
#' stock('DIA','USD',1)
#' spread('SPY.DIA','USD',c('SPY','DIA'),c(1,-1))
#' }
#' @export
synthetic.instrument <- function (primary_id, currency, members, 
                                  memberratio, ..., multiplier = 1, 
                                  tick_size=NULL, identifiers = NULL, 
                                  assign_i=TRUE, 
                                  type = c("synthetic.instrument", "synthetic")) 
{
    dargs <- list(...)
    if (!is.list(members)) {
        if (length(members) != length(memberratio) || length(members) < 2) {
            stop("length of members and memberratio must be equal, and contain two or more instruments")
        }
        memberlist <- list(members = members, memberratio = memberratio, 
                            currencies = vector(), memberpositions = NULL)
        for (member in members) {
            tmp_instr <- try(getInstrument(member, silent=TRUE))
            if (inherits(tmp_instr, "try-error") | !is.instrument(tmp_instr)) {
                if(missing(currency) || is.null(currency)) {
                    stop("'currency' must be provided if member instruments are not defined")
                    warning(paste("Instrument", member, "not found, using currency of", currency))
                } 
                memberlist$currencies[member] <- currency
            }
            else {
                memberlist$currencies[member] <- tmp_instr$currency
            }
        }

        # expires will be whichever member expires first (unless it was passed through dots)
        if (is.character(members) && is.null(dargs$expires)) {
            ids <- sort_ids(members) #sort chronologically by expiry
            expires <- NULL
            tmpinstr <- try(getInstrument(ids[1], silent=TRUE))
            if (is.instrument(tmpinstr)) expires <- tmpinstr$expires
            if (!is.null(expires) && 
                !inherits(expires, "try-error") && 
                is.null(dargs$expires)) {
                    dargs$expires <- expires
            }
        }
    } else {
        warning("passing in members as a list not fully tested")
        if (all(do.call(c, lapply(members, is.instrument)))) { #if members is a list of instruments
            instrlist <- members
            members <- do.call(c, lapply(instrlist, FUN=function(x) x$primary_id))
            memberlist <- list(members = members, memberratio = memberratio, 
                            currencies = vector(), memberpositions = NULL)
            for (i in 1:length(members)) {
                tmp_instr <- instrlist[[i]]
                memberlist$currencies[members[i]] <- tmp_instr$currency
            }
        } else {
            memberlist = members
            members <- memberlist$members
        }    
    }

    names(memberlist$members) <- memberlist$members
    names(memberlist$memberratio) <- memberlist$members
    names(memberlist$currencies) <- memberlist$members

    if (missing(primary_id) || is.null(primary_id)) 
        primary_id <- make_spread_id(members)
    if (missing(currency) || is.null(currency)) 
        currency <- as.character(memberlist$currencies[1])
    
    synthetic(primary_id = primary_id, currency = currency, multiplier = multiplier, 
        identifiers = identifiers, assign_i=assign_i, memberlist = memberlist, 
              memberratio = memberratio, tick_size=tick_size,
              ... = dargs, members = members, type = type)
}


#' @export
#' @rdname synthetic.instrument
spread <- function (primary_id = NULL, currency = NULL, members, memberratio, 
                    tick_size=NULL, ..., multiplier = 1, identifiers = NULL, 
                    assign_i=TRUE) {
    synthetic.instrument(primary_id = primary_id, currency = currency, 
      members = members, memberratio = memberratio, ...=..., tick_size=tick_size,
      multiplier = multiplier, identifiers = identifiers, assign_i=assign_i,
      type = c("spread", "synthetic.instrument", "synthetic", "instrument"))
}


#' @export
#' @rdname synthetic.instrument
butterfly <- function(primary_id = NULL, currency=NULL, members,tick_size=NULL, identifiers=NULL, assign_i=TRUE, ...)
{
##TODO: butterfly can refer to expirations (futures) or strikes (options)
##TODO: A butterfly could either have 3 members that are outrights, or 2 members that are spreads
    if (length(primary_id) > 1L) stop("primary_id must be length 1")
    if (missing(members)) {  
        pid <- parse_id(primary_id)
        root_id <- pid$root
        suffix_id <- pid$suffix
        #synthetic flies will have a root_id that looks like "
        if (suffix_id == "")
            members <- unlist(strsplit(root_id, "[-;:,\\.]"))
        else members <- paste(root_id, unlist(strsplit(suffix_id, "[-;:_,\\.]")), sep="_")
    }

    if (length(members) == 3) {
        if (is.null(currency)) {
            m1 <- getInstrument(members[[1]],silent=TRUE)
            if (!is.instrument(m1)) m1 <- getInstrument(members[[2]], silent=TRUE)
            if (!is.instrument(m1)) m1 <- getInstrument(members[[3]], silent=TRUE)
            if (is.instrument(m1)) {
                currency <- m1$currency
            } else {
                warning('currency is NULL and spread members are not defined. Using currency of USD')
                currency <- 'USD'
            }
        }
        synthetic.instrument(primary_id=primary_id,currency=currency,members=members,
            memberratio=c(1,-2,1), multiplier=1, tick_size=tick_size,
            identifiers=NULL, assign_i=assign_i, ...=..., type = c('butterfly','spread','synthetic.instrument',
            'synthetic','instrument'))
    } else if (length(members) == 2) {
        stop('butterfly currently only supports 3 leg spreads (i.e. no spread of spreads yet.)')
    } else stop('A butterfly must either have 3 outright legs or 2 spread legs')
}


#' @export
#' @rdname synthetic.instrument
guaranteed_spread <- calendar_spread <- 
    function (primary_id=NULL, currency=NULL, root_id=NULL, suffix_id=NULL, 
              members=NULL, memberratio=c(1,-1), ..., multiplier=NULL, 
              identifiers = NULL, assign_i=TRUE, 
              tick_size=NULL) {

    if (!is.null(suffix_id)) {
        if(!is.null(root_id)) {
            id<- paste(root_id,suffix_id,sep="_")
        } else {
            id <- paste(primary_id, suffix_id, sep = "_")
        }
    } else id <- primary_id

    if (is.null(id) && !is.null(members)) id <- make_spread_id(members, root=root_id)

    id<-make.names(id) #force syntactically valid primary_id

    if (is.null(suffix_id)) suffix_id <- parse_id(id)$suffix
    if (is.null(root_id)) root_id <- parse_id(id)$root

    if (is.null(members)) {
        #construct members from suffix_id and either primary_id or root_id
        members <- unlist(strsplit(suffix_id, "[-;:_,\\.]"))
        members <- paste(root_id,members, sep ="_")
    }
    
    # go get other instrument quantities from the root contract
    root_contract<-try(getInstrument(root_id,silent=TRUE,type='future'))
    if (inherits(root_contract, 'try-error')) 
        root_contract <-try(getInstrument(root_id,silent=TRUE,type='option'))
    if(is.instrument(root_contract)){
        if(is.null(currency)) currency <- root_contract$currency
        if(is.null(multiplier)) multiplier <- root_contract$multiplier
        if(is.null(tick_size)) tick_size <-  root_contract$tick_size
    } else {
        if (is.null(multiplier)) {
            message(paste(root_id, 'is not defined, using multiplier of 1'))
            multiplier <- 1
        }
        if (is.null(currency)) {
            m1 <- getInstrument(members[[1]],silent=TRUE)
            if (!is.instrument(m1)) m1 <- getInstrument(members[[2]], silent=TRUE)
            if (is.instrument(m1)) {
                currency <- m1$currency
            } else {
                warning('currency is NULL and spread members are not defined. Using currency of USD')
                currency <- 'USD'
            }
        }
    }
    
    synthetic.instrument(primary_id = id, currency = currency, members = members, 
    memberratio = memberratio, multiplier = multiplier, identifiers = identifiers, assign_i=assign_i,
    tick_size=tick_size, ... = ..., type = c("guaranteed_spread", "spread", 
    "synthetic.instrument", "synthetic", 'instrument'))
}



#' @export
#' @rdname synthetic.instrument
ICS_root <- function(primary_id, currency = NULL, members, multiplier=NULL, 
                    identifiers=NULL, assign_i=TRUE, overwrite=TRUE, 
                     tick_size=NULL, ...) {
    if (length(primary_id) < 1) stop("primary_id must be length 1")
    if (!isTRUE(overwrite) && isTRUE(assign_i) && 
            is.instrument.name(primary_id)) {
        stop("overwrite is FALSE and primary_id ", sQuote(primary_id), 
             "is already in use.")
    }

    # future roots may begin with a dot; make sure we've got the primary_ids
    members <- do.call(c, lapply(members, function(x) {
        instr <- try(getInstrument(x, type='future', silent=TRUE))
        if (is.instrument(instr)) 
            instr$primary_id
        else {
            warning(x, ' is not defined.') 
            x
        }
    }))
    getfirst <- function(chr) { # value of 'chr' field of the first of "members" that has a field named "chr"
        tmp <- suppressWarnings(try(na.omit(as.data.frame(
                buildHierarchy(members, chr), stringsAsFactors=FALSE)[[chr]][[1]])))
        if (identical(character(0), as.vector(tmp))) stop(chr, ' is required if no members are defined')
        tmp
    }

    # If currency was not given, use the currency of the first 'member' that is defined
    if (is.null(currency)) currency <- getfirst('currency')
    # do the same with multiplier and tick_size
    if (is.null(multiplier)) multiplier <- getfirst('multiplier')
    if (is.null(tick_size)) tick_size <- getfirst('tick_size')

    synthetic(primary_id, currency, multiplier, 
        identifiers=identifiers, assign_i=assign_i, 
        tick_size=tick_size, ..., members=members, type='ICS_root')
}    

#' @export
#' @rdname synthetic.instrument
ICS <- function(primary_id, assign_i=TRUE, identifiers = NULL, ...)
{ #author gsee
    pid <- parse_id(primary_id)
    if (!"ICS" %in% pid$type) stop("suffix of primary_id should look like 'H2.0302'")
 
    dargs <- list(...)
    root <- getInstrument(pid$root, silent=TRUE, type=c('ICS_root', 'spread', 'synthetic'))
    # look in dots for arguments that you can use to call ICS_root if there isn't
    #       an ICS_root already defined.
    if (!is.instrument(root)) {
        if (is.null(dargs$members)) stop(paste('Please provide "members" or define ICS_root', 
            pid$root))
        # See if we can create a temporary ICS_root with args in dots
        icsra <- list() #ICS_root args
        icsra$primary_id <- pid$root
        if (!is.null(dargs$currency)) {
            icsra$currency <- dargs$currency
            dargs$currency <- NULL
        }
        if (!is.null(dargs$multiplier)) {
            icsra$multiplier <- dargs$multiplier
            dargs$multiplier <- NULL
        }
        if (!is.null(dargs$members)) {
            icsra$members <- dargs$members
            dargs$members <- NULL
        }
        if (!is.null(dargs$tick_size)) {
            icsra$tick_size <- dargs$tick_size
            dargs$tick_size <- NULL
        }
        icsra$assign_i <- FALSE
        root <- do.call(ICS_root, icsra)
    } else {
        dargs$currency <- NULL
        dargs$multiplier <- NULL
        dargs$members <- NULL
    }
    if (!is.instrument(root)) stop("'ICS_root' must be defined first") 
    members <- root$members 
    #split the suffix in half. 1st half is CY, 2nd half is ratio string
    suff.1 <- strsplit(pid$suffix, "\\.")[[1]][1]
    suff.2 <- strsplit(pid$suffix, "\\.")[[1]][2]

    # if members are futures (roots) change them to the future_series
    # get a list of member instruments    
    memlist <- lapply(members, getInstrument, type=c('future_series', 'future'))
    #memtypes <- do.call(c, lapply(memlist, "[[", "type"))

    # if any members are future, create a future_series id, else don't change member primary_id
    members <- sapply(memlist, function(x) {
        if (x$type[1] == 'future') {
            if (is.null(x$root)) {
                paste(gsub("\\.", "", x$primary_id), suff.1, sep="_")
            } else paste(x$root, suff.1, sep="_")
        } else x$primary_id
    })

    # Check to make sure members exist in instrument envir.  Warn if not.
    defined <- sapply(members, exists, where=.instrument)
    if (any(defined == FALSE)) warning("No instrument definition found for ", 
                                       paste(members[!defined], collapse=" "))
    memberratio <- suff.2
    if (is.character(memberratio) && length(memberratio == 1)) { 
        # "0503" means c(5, -3).  "010201" is c(1,-2,1)
        memberratio <- do.call(c, lapply(seq(2, nchar(memberratio), 2), 
                    function(i) as.numeric(substr(memberratio, i-1, i))))
        # every other weight will be negative -- i.e. every other position is short
        if (length(memberratio) > 1) memberratio <- suppressWarnings(memberratio * c(1,-1)) 
    }
    #paste(sub("\\.\\.", "", members)
    if (length(dargs) == 0) dargs <- NULL
    siargs <- list() #synthetic.instrument arguments
    siargs$primary_id <- primary_id
    siargs$currency <- root$currency
    siargs$members <- members
    siargs$memberratio <- memberratio
    siargs$multiplier <- root$multiplier
    siargs$identifiers <- identifiers
    siargs$assign_i <- assign_i
    siargs$tick_size <- root$tick_size
    siargs$type <- c('ICS', 'guaranteed_spread', 'spread', 'synthetic.instrument', 'synthetic', 'instrument')
    siargs <- c(siargs, dargs)
    do.call(synthetic.instrument, siargs)
}


