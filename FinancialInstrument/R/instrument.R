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
# $Id: instrument.R 1655 2014-11-23 22:53:26Z gsee $
#
###############################################################################

#.onLoad <- function(lib, pkg) {
#    if(!exists('.instrument'))
#        .instrument <<- new.env(hash=TRUE)
#}

.instrument <- new.env(parent=emptyenv())


#' class test for object supposedly of type 'instrument'
#' @param x object to test for type
#' @export
is.instrument <- function( x ) {
  inherits( x, "instrument" )
}


#' check each element of a character vector to see if it is either the 
#' primary_id or an identifier of an \code{\link{instrument}}
#' @param x character vector
#' @return logical vector
#' @export
is.instrument.name <- function(x) {
  if (!is.character(x)) return(FALSE)
  sapply(lapply(x, getInstrument, silent=TRUE), inherits, "instrument")
}


#' class test for object supposedly of type 'currency'
#' @param x object to test for type
#' @export
is.currency <- function( x ) {
#  x<-getInstrument(x, silent=TRUE) # Please use is.currency.name if x is character
  inherits( x, "currency" )
}


#' check each element of a character vector to see if it is either the 
#' primary_id or an identifier of a \code{\link{currency}}
#' @param x character vector
#' @export
is.currency.name <- function( x ) {
  if (!is.character(x)) return(FALSE)
  sapply(lapply(x, getInstrument, type='currency', silent=TRUE), inherits, 
         "currency")
}


#' instrument class constructors
#' 
#' All 'currency' instruments must be defined before instruments of other types 
#' may be defined.
#' 
#' In \dots you may pass any other arbitrary instrument fields that will be used 
#' to create 'custom' fields.  S3 classes in \R are basically lists with a class 
#' attribute.  We use this to our advantage to allow us to set arbitrary fields.
#' 
#' \code{identifiers} should be a named list to specify other identifiers beyond 
#' the \code{primary_id}.  Please note that whenever possible, these should 
#' still be unique.  Perhaps Bloomberg, Reuters-X.RIC, CUSIP, etc.
#' \code{\link{getInstrument}} will return the first (and only the first) match 
#' that it finds, starting with the primary_id, and then searching the 
#' primary_ids of all instruments for each of the \code{identifiers}.  Note that
#' when a large number of instruments are defined, it is faster to find 
#' instruments by \code{primary_id} than by \code{identifiers} because it looks
#' for \code{primary_id}s first.
#' 
#' The \code{primary_id} will be coerced within reason to a valid \R variable 
#' name by using \code{\link{make.names}}. We also remove any leading '1' digit 
#' (a simple workaround to account for issues with the Reuters API).  If you are
#' defining an instrument that is not a \code{currency}, with a primary_id that
#' already belongs to a \code{currency}, a new primary_id will be create using
#' \code{make.names}.  For example, \code{stock("USD", currency("USD"))}, would
#' create a stock with a primary_id of \dQuote{USD.1} instead of overwritting
#' the \code{currency}.
#'
#' Please use some care to choose your primary identifiers so that R won't 
#' complain.  If you have better regular expression code, we'd be happy to 
#' include it.   
#' 
#' Identifiers will also try to be discovered as regular named arguments passed 
#' in via \code{...}.  We currently match any of the following: 
#' \code{"CUSIP","SEDOL","ISIN","OSI","Bloomberg","Reuters","X.RIC","CQG","TT","Yahoo","Google"}
#' Others may be specified using a named list of identifiers, as described above.
#' 
#' \code{assign_i} will use \code{\link{assign}} to place the constructed 
#' instrument class object into the \code{.instrument} environment.  Most of the 
#' special type-specific constructors will use \code{assign_i=TRUE} internally. 
#' Calling with \code{assign_i=FALSE}, or not specifying it, will return an 
#' object and will \emph{not} store it.  Use this option ether to wrap calls to 
#' \code{instrument} prior to further processing (and presumably assignment) or 
#' to test your parameters before assignment.
#' 
#' If \code{overwrite=FALSE} is used, an error will be thrown if any 
#' \code{primary_id}s are already in use.
#' 
#' As of version 0.10.0, the .instrument environment is located at the top level
#' of the package. i.e. \code{.instrument}.
#' 
#' \code{future} and \code{option} are used to define the contract specs of a 
#' series of instruments.  The \code{primary_id} for these can begin with 1 or 
#' 2 dots if you need to avoid overwriting another instrument.
#' For example, if you have a \code{stock} with \sQuote{SPY} as the 
#' \code{primary_id}, you could use \sQuote{.SPY} as the \code{primary_id} of 
#' the \code{option} specs, and \sQuote{..SPY} as the \code{primary_id} of the 
#' single stock \code{future} specs. (or vice versa)
#'
#' You can (optionally) provide a \code{src} argument in which case, it will be 
#' used in a call to \code{\link[quantmod]{setSymbolLookup}}.
#' @param primary_id String describing the unique ID for the instrument. Most 
#'   of the wrappers allow this to be a vector.
#' @param ... Any other passthru parameters, including 
#' @param underlying_id For derivatives, the identifier of the instrument that 
#'   this one is derived from, may be \code{NULL} for cash settled instruments
#' @param currency String describing the currency ID of an object of type 
#'   \code{\link{currency}}
#' @param multiplier Numeric multiplier to apply to the price in the instrument 
#'   to get to notional value.
#' @param tick_size The tick increment of the instrument price in it's 
#'   trading venue, as numeric quantity (e.g. 1/8 is .125)
#' @param identifiers Named list of any other identifiers that should also be 
#'   stored for this instrument
#' @param type instrument type to be appended to the class definition, typically 
#'   not set by user
#' @param assign_i TRUE/FALSE. Should the instrument be assigned to the 
#'   \code{.instrument} environment?  Default is FALSE for \code{instrument}, 
#'   TRUE for wrappers.
#' @param overwrite TRUE/FALSE. Should existing instruments with the same
#'   primary_id be overwritten? Default is TRUE. If FALSE, an error will be 
#'   thrown and the instrument will not be created.
#' @aliases 
#' stock
#' bond
#' future
#' option
#' currency
#' instrument
#' fund
#' @seealso 
#' \code{\link{currency}},
#' \code{\link{exchange_rate}},
#' \code{\link{option_series}},
#' \code{\link{future_series}},
#' \code{\link{spread}},
#' \code{\link{load.instruments}}
#' @export
instrument <- function(primary_id , ..., currency , multiplier , tick_size=NULL, 
                     identifiers = NULL, type=NULL, assign_i=FALSE, 
                     overwrite=TRUE) {
  if(is.null(primary_id)) {
      stop("you must specify a primary_id for the instrument")
  }
  raw_id <- primary_id
  #deal with leading digits or illegal characters
  if(substr(primary_id,1,1)==1) {
      primary_id <- substr(primary_id,2,nchar(primary_id))
  }
  primary_id<-make.names(primary_id)
  
  if(missing(currency) || is.null(currency) || 
        (!missing(currency) && !is.currency.name(currency))) {
      stop("currency ",currency," must be defined first")
  }
  if(!hasArg(identifiers) || is.null(identifiers)) identifiers = list()
  if(!is.list(identifiers)) {
      warning("identifiers",identifiers,"do not appear to be a named list")
  } 

  if (raw_id != primary_id) {
      identifiers <- c(identifiers, raw_id=raw_id)
  }
  arg<-list(...)
  if(is.list(arg[['...']])){
      if(length(arg)==1) arg <- arg[['...']]
      else {
          targ<-arg[['...']]
          arg[['...']]<-NULL
          arg<-c(arg,targ)
      }
  }
  if (!is.null(arg$src)) {
      sarg <- list()
      sarg[[primary_id]] <- arg$src
      setSymbolLookup(sarg)
      #arg[["src"]]<-NULL
  }
  #check for identifiers we recognize 
  ident_str<-tolower(c("X.RIC", "RIC", "CUSIP", "SEDOL", "OSI", "Bloomberg", 
                       "Reuters", "ISIN", "CQG", "TT", "Yahoo", "Google")) #converted to lowercase for easier case-insensitive matching
  lnarg <- tolower(names(arg)) #lower case names of arguments
  pos_arg <- which(lnarg %in% ident_str)
  identifiers <- c(identifiers, arg[pos_arg])
  arg[pos_arg] <- NULL
  
  ## TODO note that multiplier could be a time series, probably add code here to check
  if(!is.numeric(multiplier) || length(multiplier) > 1) {
      stop("multiplier must be a single number")
  }
  if(!is.null(tick_size) && (!is.numeric(tick_size) | length(tick_size) > 1)) {
      stop("tick_size must be NULL or a single number")
  }  
  if(is.null(type)) {
      tclass="instrument" 
  } else tclass = unique(c(type,"instrument"))
  
  if (is.currency.name(primary_id) && 
          !inherits(getInstrument(primary_id, type="currency"), 
                    "exchange_rate")) {
      oid <- primary_id
      primary_id <- tail(make.names(c(ls_instruments(), oid), unique=TRUE), 1)
      warning(paste(oid, "is the name of a currency. Using", primary_id, 
                    "for the primary_id of this", type))
      identifiers <- c(identifiers, ticker=oid)
  } else if ((primary_id %in% ls_instruments()) && !overwrite && 
                 isTRUE(assign_i)) {
      # primary_id already exists and we are not overwriting
	  stop(paste("an instrument with primary_id", primary_id, 
                 "already exists in the .instrument environment.",
                 "Set overwrite=TRUE to overwrite."))
  }
  tmpinstr <- list(primary_id = primary_id,
                   currency = currency,
                   multiplier = multiplier,
                   tick_size=tick_size,
                   identifiers = identifiers,
                   type = type)
  if(length(arg)>=1) {
      tmpinstr <- c(tmpinstr,arg)   
  }
  class(tmpinstr)<-tclass
  
  if(assign_i)  {
      assign(primary_id, tmpinstr, 
             envir=as.environment(.instrument) )
      return(primary_id)  
  } else return(tmpinstr) 
}

#' @export
#' @rdname instrument
stock <- function(primary_id , currency=NULL , multiplier=1 , tick_size=.01, 
                  identifiers = NULL, assign_i=TRUE, overwrite=TRUE, ...){
    if (is.null(currency)) stop ("'currency' is a required argument")
    if (!isTRUE(overwrite) && isTRUE(assign_i) &&
        any(in.use <- primary_id %in% (li <- ls_instruments()))) {
        stop(paste(paste("In stock(...) : ",
                          "overwrite is FALSE and primary_id", 
                          if (sum(in.use) > 1) "s are" else " is", 
                          " already in use:\n", sep=""),
                   paste(intersect(primary_id, li), collapse=", ")), 
             call.=FALSE)
    }
    if (length(primary_id) > 1) {
        out <- sapply(primary_id, stock, currency=currency, 
                      multiplier=multiplier, tick_size=tick_size, 
                      identifiers=identifiers, assign_i=assign_i,
                      ...=..., simplify=assign_i)
        return(if (assign_i) unname(out) else out)
    }
    instrument(primary_id=primary_id, currency=currency, multiplier=multiplier, 
               tick_size=tick_size, identifiers = identifiers, ..., 
               type="stock", assign_i=assign_i)
}

#' @export
#' @rdname instrument
fund <- function(primary_id , currency=NULL , multiplier=1 , tick_size=.01, 
                 identifiers = NULL, assign_i=TRUE, overwrite=TRUE, ...){
    if (is.null(currency)) stop ("'currency' is a required argument")
    if (!isTRUE(overwrite) && isTRUE(assign_i) &&
        any(in.use <- primary_id %in% (li <- ls_instruments()))) {
        stop(paste(paste("In fund(...) : ",
                          "overwrite is FALSE and primary_id", 
                          if (sum(in.use) > 1) "s are" else " is", 
                          " already in use:\n", sep=""),
                   paste(intersect(primary_id, li), collapse=", ")), 
             call.=FALSE)
    }
    if (length(primary_id) > 1) {
        out <- sapply(primary_id, fund, currency=currency, 
                      multiplier=multiplier, tick_size=tick_size, 
                      identifiers=identifiers, assign_i=assign_i, ...=..., 
                      simplify=assign_i)
        return(if (assign_i) unname(out) else out)
    }
    instrument(primary_id=primary_id, currency=currency, 
               multiplier=multiplier, tick_size=tick_size, 
               identifiers=identifiers, ..., type="fund", assign_i=assign_i)
}

#' @export
#' @rdname instrument
future <- function(primary_id , currency , multiplier , tick_size=NULL, 
                   identifiers = NULL, assign_i=TRUE, overwrite=TRUE, ..., 
                   underlying_id=NULL){
    if(missing(primary_id)) primary_id <- paste("..",underlying_id,sep="")
    if (length(primary_id) > 1) stop('primary_id must be of length 1')
    if (!isTRUE(overwrite) && assign_i==TRUE && 
            primary_id %in% ls_instruments()) {
        stop(sQuote(primary_id), " already in use and overwrite=FALSE")
    }
    if (missing(currency) && !is.null(underlying_id)) {
        uinstr <- getInstrument(underlying_id,silent=TRUE)
        if (is.instrument(uinstr)) {
            currency <- uinstr$currency
        } else stop("'currency' is a required argument")
    }
    if(is.null(underlying_id)) {
        warning("underlying_id should only be NULL for cash-settled futures")
    } else {
        if(!exists(underlying_id, where=.instrument,
                   inherits=TRUE)) {
            warning("underlying_id not found") # assumes that we know where to look
        }
        if (primary_id == underlying_id) {
            primary_id <- paste("..",primary_id,sep="")
            warning(paste('primary_id is the same as underlying_id,',
                'the instrument will be given a primary_id of', primary_id))
        }  
    }

    instrument(primary_id=primary_id, currency=currency, multiplier=multiplier, 
               tick_size=tick_size, identifiers = identifiers, ... , 
               type="future", underlying_id=underlying_id, assign_i=assign_i)
}


#' Constructors for series contracts
#' 
#' Constructors for series contracts on instruments such as options and futures
#'
#' The root \code{instrument} (e.g. the \code{future} or \code{option}) must be
#' defined first.
#'
#' In custom parameters for these series contracts, we have often found it
#' useful to store attributes such as local roll-on and roll-off dates
#' (rolling not on the \code{first_listed} or \code{expires}.  
#'
#' For \code{future_series} and \code{option_series} you may either provide a 
#'   \code{primary_id} (or vector of \code{primary_id}s), 
#' OR both a \code{root_id} and \code{suffix_id}.
#'
#' Note that the code for \code{bond} and \code{bond_series} has not been 
#' updated recently and may not support all the features supported for
#' \code{option_series} and \code{future_series}.  Patches welcome.
#'
#' @param primary_id String describing the unique ID for the instrument. May be 
#'   a vector for \code{future_series} and \code{option_series}
#' @param root_id String product code or underlying_id, usually something like 
#'   'ES' or 'CL' for futures, or the underlying stock symbol (maybe preceded 
#'   with a dot) for equity options.
#' @param suffix_id String suffix that should be associated with the series, 
#'   usually something like 'Z9' or 'Mar10' denoting expiration and year.
#' @param first_traded String coercible to Date for first trading day.
#' @param expires String coercible to Date for expiration date
#' @param maturity String coercible to Date for maturity date of bond series.
#' @param callput Right of option; call or put
#' @param strike Strike price of option
#' @param payment_schedule Not currently being implemented
#' @param identifiers Named list of any other identifiers that should also be 
#'   stored for this instrument.
#' @param assign_i TRUE/FALSE. Should the instrument be assigned in the 
#'   \code{.instrument} environment?
#' @param overwrite TRUE/FALSE. If FALSE, only \code{first_traded} and 
#'   \code{expires} will be updated.
#' @param ... any other passthru parameters
#' @aliases 
#' option_series
#' future_series
#' bond_series
#' @examples
#' \dontrun{
#' currency("USD")
#' future("ES","USD",multiplier=50, tick_size=0.25)
#' future_series('ES_U1')
#' future_series(root_id='ES',suffix_id='Z11')
#' stock('SPY','USD')
#' option('.SPY','USD',multiplier=100,underlying_id='SPY')
#' #can use either .SPY or SPY for the root_id. 
#' #it will find the one that is option specs.
#' option_series('SPY_110917C125', expires='2011-09-16')
#' option_series(root_id='SPY',suffix_id='111022P125')
#' option_series(root_id='.SPY',suffix_id='111119C130')
#' #multiple series instruments at once.
#' future_series(c("ES_H12","ES_M12"))
#' option_series(c("SPY_110917C115","SPY_110917P115"))
#' }
#' @export
#' @rdname series_instrument
future_series <- function(primary_id, root_id=NULL, suffix_id=NULL, 
                          first_traded=NULL, expires=NULL, identifiers = NULL, 
                          assign_i=TRUE, overwrite=TRUE, ...){
  # if overwrite==FALSE and assign_i==TRUE, we'll need to know what instruments
  # are already defined.  Don't bother doing this if we're overwriting anyway
  if (!isTRUE(overwrite) && isTRUE(assign_i)) li <- ls_instruments()
  if (missing(primary_id)) {
      if (all(is.null(c(root_id,suffix_id)))) {
          stop(paste('must provide either a primary_id or',
                     'both a root_id and a suffix_id'))
      } else {
          if (is.null(suffix_id)) {
              sdate <- gsub("-","",expires)
              if (is.null(expires) || nchar(sdate) < 6) {
                  stop("must provide either 'expires' or 'suffix_id'")
              }
              suffix_id <- paste(M2C()[as.numeric(substr(sdate,5,6))], 
                                                  substr(sdate,3,4),sep="")              
          }
          primary_id <- paste(gsub("\\.","",root_id), suffix_id, sep="_")
      }
  } else if (length(primary_id) > 1) {
      if (!is.null(expires) || !is.null(first_traded)) {
          stop(paste("'first_traded' and 'expires' must be NULL",
                     "if calling with multiple primary_ids"))
      }
      if (!isTRUE(overwrite) && isTRUE(assign_i) &&
          any(in.use <- primary_id %in% li)) {
          stop(paste(paste("In future_series(...) : ",
                            "overwrite is FALSE and primary_id", 
                            if (sum(in.use) > 1) "s are" else " is", 
                            " already in use:\n", sep=""),
                   paste(intersect(primary_id, li), collapse=", ")), 
               call.=FALSE)
      }
      out <- sapply(primary_id, future_series, root_id=root_id, 
                    suffix_id=suffix_id, first_traded=first_traded, 
                    expires=expires, identifiers = identifiers, 
                    assign_i=assign_i, ...=..., simplify=assign_i)
      return(if (assign_i) unname(out) else out)
  } else if (is.null(root_id) && !is.null(suffix_id) && 
             parse_id(primary_id)$type == 'root') {
      #if we have primary_id, but primary_id looks like a root_id, and we have 
      #suffix_id and don't have root_id then primary_id is really root_id and we 
      #need to replace primary_id
      root_id <- primary_id
      primary_id <- paste(root_id, suffix_id, sep="_")
  }
  if (!isTRUE(overwrite) && isTRUE(assign_i) && primary_id %in% li) {
      stop(sQuote(primary_id), " already in use and overwrite=FALSE")
  }
  pid <- parse_id(primary_id)
  if (is.null(root_id)) root_id <- pid$root
  if (is.null(suffix_id)) suffix_id <- pid$suffix
  if (is.null(expires)) {
    expires <- paste(pid$year, 
                     sprintf("%02d", match(pid$month, toupper(month.abb))),
                     sep='-') 
    #if expires now has an NA in it, set it back to NULL
    if (!identical(integer(0), grep("NA",expires))) expires <- NULL
  }

  contract<-getInstrument(root_id,type='future')
  
  # TODO add check for Date equivalent in first_traded and expires

  ## with futures series we probably need to be more sophisticated,
  ## and find the existing series from prior periods (probably years or months)
  ## and then add the first_traded and expires to the time series bu splicing
  #temp_series<-try(getInstrument(primary_id, silent=TRUE),silent=TRUE)
  if (!isTRUE(overwrite)) {
      temp_series<-try(getInstrument(primary_id, silent=TRUE),silent=TRUE)
      if(inherits(temp_series,"future_series")) {
          message("updating existing first_traded and expires for ",primary_id)
          temp_series$first_traded <- unique(c(temp_series$first_traded,
                                               first_traded))
          temp_series$expires<-unique(c(temp_series$expires,expires))
          assign(primary_id, temp_series, 
                 envir=as.environment(.instrument))
          return(primary_id)
      } else warning("No contract found to update. A new one will be created.")
  }
  args <- list()
  args$primary_id <- primary_id
  args$root_id <- root_id
  args$suffix_id=suffix_id
  args$currency = contract$currency
  args$multiplier = contract$multiplier
  args$tick_size=contract$tick_size
  args$identifiers = identifiers
  args$first_traded = first_traded
  args$type=c("future_series", "future")
  args$expires = expires
  if (!is.null(contract$exchange)) {
      args$exchange <- contract$exchange
  }
  args$underlying_id = contract$underlying_id
  if (!is.null(contract$marketName)) {
      args$marketName <- contract$marketName
  }
  if (!is.null(contract$exchange_id)) {
      args$exchange_id <- contract$exchange_id
  }
  if (!is.null(contract$description)) {
      args$series_description <- paste(contract$description, expires)
  }
  args$assign_i=assign_i
  dargs<-list(...)
  dargs$currency=NULL
  dargs$multiplier=NULL
  dargs$type=NULL
  if (is.null(dargs$src) && !is.null(contract$src)){
      dargs$src <- contract$src
  }
  args <- c(args, dargs)

  do.call(instrument, args)
}

#' @export
#' @rdname instrument
option <- function(primary_id , currency , multiplier , tick_size=NULL, 
                   identifiers = NULL, assign_i=TRUE, overwrite=TRUE,
                   ..., underlying_id=NULL){
  if (missing(primary_id)) primary_id <- paste(".",underlying_id,sep="")
  if (length(primary_id) > 1) stop("'primary_id' must be of length 1")
  if (!isTRUE(overwrite) && assign_i==TRUE && 
          primary_id %in% ls_instruments()) {
      stop(sQuote(primary_id), " already in use and overwrite=FALSE")
  }
  if (missing(currency) && !is.null(underlying_id)) {
        uinstr <- getInstrument(underlying_id,silent=TRUE)
        if (is.instrument(uinstr)) {
            currency <- uinstr$currency
        } else stop("'currency' is a required argument")
  }
  if(is.null(underlying_id)) {
      warning("underlying_id should only be NULL for cash-settled options")
  } else {
      if(!exists(underlying_id, where=.instrument,
                 inherits=TRUE)) {
          warning("underlying_id not found") # assumes that we know where to look
      }
      if (primary_id == underlying_id) {
          primary_id <- paste(".",primary_id,sep="")
          warning(paste('primary_id is the same as underlying_id,',
                'the instrument will be given a primary_id of', primary_id))
      }  
  }
  ## now structure and return
  instrument(primary_id=primary_id, currency=currency, multiplier=multiplier, 
             tick_size=tick_size, identifiers = identifiers, ..., type="option", 
             underlying_id=underlying_id, assign_i=assign_i )
}

#' @export
#' @rdname series_instrument
option_series <- function(primary_id , root_id = NULL, suffix_id = NULL, 
                          first_traded=NULL, expires=NULL, 
                          callput=c("call","put"), strike=NULL, 
                          identifiers=NULL, assign_i=TRUE, overwrite=TRUE, ...){
    if (!isTRUE(overwrite) && isTRUE(assign_i)) li <- ls_instruments()
    if (missing(primary_id) ) {
        if (all(is.null(c(root_id,suffix_id)))) {
            stop(paste('must provide either a primary_id or',
                       'both a root_id and a suffix_id'))
        } else { #if you give it only a root_id it will make the suffix_id 
                 #using expires, callput, and strike
            if (is.null(suffix_id)) {
                sdate <- if (nchar(expires) == 8) { 
                    try(as.Date(expires, format='%Y%m%d'), silent=TRUE) 
                } else try(as.Date(expires),silent=TRUE)
                if (inherits(sdate,'try-error')) {
                    stop("expires is missing or of incorrect format")
                }
                sright <- try(switch(callput, C=,c=,call="C", P=,p=,put="P"),
                              silent=TRUE)
                if (inherits(sright,'try-error')) {
                    stop(paste("must provide 'callput' or a 'suffix_id'",
                               "from which 'callput' can be inferred."))
                }
                if (is.null(strike)) {
                    stop(paste("must provide 'strike' or a 'suffix_id'",
                               "from which 'strike' can be inferred."))
                }
                suffix_id <- paste(format(sdate,'%y%m%d'), sright, strike, 
                                   sep="")
            }
            primary_id <- paste(gsub("\\.","",root_id), suffix_id, sep="_")
        }
    } else if (length(primary_id) > 1) {
       if (!isTRUE(overwrite) && isTRUE(assign_i) && 
               any(in.use <- primary_id %in% li)) {
          stop(paste(paste("In option_series(...) : ",
                            "overwrite is FALSE and primary_id", 
                            if (sum(in.use) > 1) "s are" else " is", 
                            " already in use:\n", sep=""),
                     paste(intersect(primary_id, li), collapse=", ")), 
               call.=FALSE)
      }
      if (!is.null(expires) || !is.null(first_traded)) {
          stop(paste("'first_traded' and 'expires' must be NULL",
                     "if calling with multiple primary_ids"))
      }
      out <- sapply(primary_id, option_series, root_id=root_id, 
                    suffix_id=suffix_id, first_traded=first_traded, 
                    expires=expires, callput=callput, strike=strike, 
                    identifiers=identifiers, assign_i=assign_i, ...=..., 
                    simplify=assign_i)
      return(if (assign_i) unname(out) else out)
    } else if (is.null(root_id) && !is.null(suffix_id) && 
               parse_id(primary_id)$type == 'root') {
          #if we have primary_id, but primary_id looks like a root_id, and we 
          #have suffix_id and don't have root_id then primary_id is really 
          #root_id and we need to replace primary_id
          root_id <- primary_id
          primary_id <- paste(root_id, suffix_id, sep="_")
    }
    if (!isTRUE(overwrite) && isTRUE(assign_i) && primary_id %in% li) {
        stop(sQuote(primary_id), " already in use and overwrite=FALSE")
    }
    pid <- parse_id(primary_id)
    if (is.null(root_id)) root_id <- pid$root
    if (is.null(suffix_id)) suffix_id <- pid$suffix
    if (is.null(strike)) {
        if (is.na(pid$strike)) stop('strike must be provided.')
        strike <- pid$strike
    }
    if (is.null(expires)) {
        #TODO: option_series ids contain the entire date.  
        #      Don't have to settle for YYYY-MM
        expires <- paste(pid$year, 
                         sprintf("%02d",match(pid$month, 
                                              toupper(month.abb))),sep='-') 
        if (!identical(integer(0), grep("NA",expires))) {
            stop(paste("must provide 'expires' formatted '%Y-%m-%d',",
                       "or a 'suffix_id' from which to infer 'expires'"))
        }
    }
    contract<-getInstrument(root_id, type='option')
    ## with options series we probably need to be more sophisticated,
    ## and find the existing series from prior periods (probably years)
    ## and then add the first_traded and expires to the time series
    if(length(callput)==2) callput <- switch(pid$right, C='call', P='put')
    if (is.null(callput)) {
        stop("value of callput must be specified as 'call' or 'put'")
    }
    if (!isTRUE(overwrite)) {
        temp_series<-try(getInstrument(primary_id, silent=TRUE),silent=TRUE)
        if(inherits(temp_series,"option_series")) {
            message("updating existing first_traded and expires for ", 
                    primary_id)
            temp_series$first_traded <- unique(c(temp_series$first_traded,
                                                 first_traded))
            temp_series$expires<-unique(c(temp_series$expires,expires))
            assign(primary_id, temp_series, 
                   envir=as.environment(.instrument))
            return(primary_id)
        } else {
            warning("No contract found to update.  A new one will be created.")
        }
    } else {
        dargs <- list(...)
        if (is.null(dargs$src) && !is.null(contract$src)) {
            dargs$src <- contract$src
        }
        instrument( primary_id = primary_id,
                    root_id = root_id,
                    suffix_id = suffix_id,
                    currency = contract$currency,
                    multiplier = contract$multiplier,
                    tick_size=contract$tick_size,
                    first_traded = first_traded,
                    expires = expires,
                    identifiers = identifiers,
                    callput = callput,
                    strike = strike,
                    underlying_id = contract$underlying_id,
                    ...=dargs,
                    type=c("option_series", "option"),
                    assign_i=assign_i
                  ) 
    }
}

#' constructor for series of options using yahoo data
#'
#' Defines a chain or several chains of options by looking up necessary info 
#' from yahoo.  
#'
#' If \code{Exp} is missing it will define only the nearby options. 
#' If \code{Exp} is NULL it will define all options
#' 
#' If \code{first_traded} and/or \code{tick_size} should not be the same for all 
#' options being defined, they should be left NULL and defined outside of this 
#' function.
#' @param symbol character vector of ticker symbols of the underlying 
#'   instruments (Currently, should only be stock tickers)
#' @param Exp Expiration date or dates to be passed to getOptionChain
#' @param currency currency of underlying and options
#' @param multiplier contract multiplier. Usually 100 for stock options
#' @param first_traded first date that contracts are tradeable. Probably not 
#'   applicable if defining several chains.
#' @param tick_size minimum price change of options.
#' @param overwrite if an instrument already exists, should it be overwritten?
#' @return Called for side-effect. The instrument that is created and stored 
#'   will inherit option_series, option, and instrument classes. 
#' @references Yahoo \url{http://finance.yahoo.com}
#' @author Garrett See
#' @note Has only been tested with stock options.
#' The options' currency should be the same as the underlying's.
#' @seealso \code{\link{option_series}}, \code{\link{option}}, 
#'   \code{\link{instrument}}, \code{\link[quantmod]{getOptionChain}}
#' @examples
#' \dontrun{
#' option_series.yahoo('SPY') #only nearby calls and puts
#' option_series.yahoo('DIA', Exp=NULL) #all chains
#' ls_instruments()
#' }
#' @export
option_series.yahoo <- function(symbol, Exp, currency="USD", multiplier=100, 
                                first_traded=NULL, tick_size=NULL, overwrite=TRUE) {
    #FIXME: identifiers?

    opts <- getOptionChain(Symbols=symbol,Exp=Exp, src="yahoo")

    locals <- function(x) c(grep(symbol, rownames(x$puts), value=TRUE),
                            grep(symbol, rownames(x$calls), value=TRUE))
    if (is.null(opts$calls)) { #if is.null(Exp) we'll get back all chains
        led <- (lapply(opts, locals))  
        optnames <- unname(do.call(c, led)) #FIXME: Is this a reasonable way to get rownames?
    } else optnames <- locals(opts) #c(rownames(opts$calls),rownames(opts$puts))


    CleanID <- function(x, symbol) {
        si <- gsub(symbol, "", x) #suffix_id        
        out <- list(root_id = symbol,
                    expiry = substr(si, 1, 6),
                    right = substr(si, 7, 7),
                    strike = as.numeric(substr(si, 8, 15))/1000)
        clean.si <- with(out, paste(expiry, right, strike, sep=""))
        c(out, list(clean.si=clean.si, 
                    primary_id = paste(symbol, "_", clean.si, sep="")))
    }

    id.list <- lapply(optnames, CleanID, symbol)
    
    if (!isTRUE(overwrite)) {
        new.ids <- unname((u <- unlist(id.list))[grep("primary_id", names(u))])
        if (any(in.use <- new.ids %in% (li <- ls_instruments()))) {
            stop(paste(paste("In option_series.yahoo(...) : ",
                              "overwrite is FALSE and primary_id", 
                              if (sum(in.use) > 1) "s are" else " is", 
                              " already in use:\n", sep=""),
                       paste(intersect(new.ids, li), collapse=", ")), 
                 call.=FALSE)
        }
    }
    
    idout <- NULL
    for (ID in id.list) {
        #create currency if it doesn't exist #?? Any reason not to ??
        tmpccy <- try(getInstrument(currency, silent=TRUE), silent=TRUE)
        if (!inherits(tmpccy, "currency")) {
            warning(paste("Created currency", currency, 
                          "because it did not exist."))
            currency(currency) #create it
        }
        #create option spec if we need to.
        tmpInstr <- try(getInstrument(paste('.',symbol,sep=""), silent=TRUE),
                        silent=TRUE)
        if (!inherits(tmpInstr, "option")) {
            warning(paste('Created option specs for root',
                              paste('.', symbol, sep="")))
            option(primary_id=paste('.',symbol,sep=""), currency=currency,
                multiplier=multiplier, tick_size=tick_size, 
                underlying_id=symbol)
        }
        #create specific option
        tempseries = instrument(primary_id=ID[["primary_id"]], 
                                suffix_id=ID[["clean.si"]], 
                                first_traded=first_traded, 
                                currency=currency, 
                                multiplier=multiplier, 
                                tick_size=tick_size, 
                                expires=as.Date(paste(paste('20', 
                                                            substr(ID[["expiry"]], 1, 2),
                                                            sep=""), 
                                                      substr(ID[["expiry"]], 3, 4), 
                                                      substr(ID[["expiry"]], 5, 6),
                                                      sep="-")), 
                                callput=switch(ID[["right"]], C="call", P="put"), #to be consistent with the other option_series function
                                strike=ID[["strike"]], 
                                underlying_id=symbol, 
                                type = c("option_series","option"), 
                                defined.by='yahoo', assign_i=TRUE, overwrite=overwrite
                                )    
        idout <- c(idout, ID[["primary_id"]])
    }
    idout
}

#' @export
#' @rdname instrument
currency <- function(primary_id, identifiers = NULL, assign_i=TRUE, ...){
    if (hasArg("overwrite")) {
        if (!list(...)$overwrite && isTRUE(assign_i) &&
            any(in.use <- primary_id %in% (li <- ls_instruments()))) {
            stop(paste(paste("In currency(...) : ",
                              "overwrite is FALSE and primary_id", 
                              if (sum(in.use) > 1) "s are" else " is", 
                              " already in use:\n", sep=""),
                       paste(intersect(primary_id, li), collapse=", ")), 
                call.=FALSE)
        }
    }
    if (length(primary_id) > 1) {
        out <- sapply(primary_id, currency, identifiers=identifiers, 
                      assign_i=assign_i, ...=..., simplify=assign_i)
        return(if (assign_i) unname(out) else out)
    }
    if (is.null(identifiers)) identifiers <- list()
    ccy <- try(getInstrument(primary_id,type='currency',silent=TRUE))
    if (is.instrument(ccy)) {
        if (length(identifiers) > 0) {
            if (!is.list(identifiers)) identifiers <- list(identifiers)
            for (nm in names(ccy$identifiers)[names(ccy$identifiers) %in% 
                                              names(identifiers)]) {
                ccy$identifiers[[nm]] <- identifiers[[nm]]
            }
            identifiers <- identifiers[names(identifiers)[!names(identifiers) 
                               %in% names(ccy$identifiers)]]
            ccy$identifiers <- c(identifiers, ccy$identifiers)
        }
    } else ccy <- list(primary_id = primary_id,
                        currency = primary_id,
                        multiplier = 1,
                        tick_size= .01,
                        identifiers = identifiers,
                        type = "currency")
    dargs <- list(...)
    if (!is.null(dargs)) {
        for (nm in names(ccy)[names(ccy) %in% names(dargs)]) {
            ccy[[nm]] <- dargs[[nm]]
        }
        dargs <- dargs[names(dargs)[!names(dargs) %in% names(ccy)]]
        ccy <- c(ccy,dargs)
    }        
    class(ccy)<-c("currency","instrument")
    if (assign_i) {
        assign(primary_id, ccy, 
               pos=as.environment(.instrument) )
        return(primary_id)
    }
    ccy
}


#' constructor for spot exchange rate instruments
#' 
#' Currency symbols (like any symbol) may be any combination of alphanumeric 
#' characters, but the FX market has a convention that says that the first 
#' currency in a currency pair is the 'target'  and the second currency in the 
#' symbol pair is the currency the rate ticks in.  So 'EURUSD' can be read as 
#' 'USD per 1 EUR'.
#' 
#' In \code{FinancialInstrument} the \code{currency} of the instrument should 
#' be the currency that the spot rate ticks in, so it will typically be the 
#' second currency listed in the symbol. 
#' 
#' Thanks to Garrett See for helping sort out the inconsistencies in different 
#' naming and calculating conventions. 
#' @param primary_id string identifier, usually expressed as a currency pair 
#'   'USDYEN' or 'EURGBP'
#' @param currency string identifying the currency the exchange rate ticks in
#' @param counter_currency string identifying the currency which the rate uses 
#'   as the base 'per 1' multiplier
#' @param tick_size minimum price change
#' @param identifiers named list of any other identifiers that should also be 
#'   stored for this instrument
#' @param assign_i TRUE/FALSE. Should the instrument be assigned in the 
#'   \code{.instrument} environment? (Default TRUE)
#' @param overwrite \code{TRUE} by default.  If \code{FALSE}, an error will
#'   be thrown if there is already an instrument defined with the same 
#'   \code{primary_id}.
#' @param ... any other passthru parameters
#' @references http://financial-dictionary.thefreedictionary.com/Base+Currency
#' @export
exchange_rate <- function (primary_id = NULL, currency = NULL, 
                           counter_currency = NULL, tick_size=0.01, 
                           identifiers = NULL, assign_i=TRUE, overwrite=TRUE, 
                           ...){
  if (is.null(primary_id) && !is.null(currency) && !is.null(counter_currency)) {
    primary_id <- c(outer(counter_currency,currency,paste,sep=""))
    same.same <- function(x) substr(x,1,3) == substr(x,4,6)
    primary_id <- primary_id[!same.same(primary_id)]
  } else if (is.null(primary_id) && (is.null(currency) || 
             is.null(counter_currency))) {
    stop(paste("Must provide either 'primary_id' or both",
               "'currency' and 'counter_currency'"))
  }
  if (!isTRUE(overwrite) && isTRUE(assign_i) &&
        any(in.use <- primary_id %in% (li <- ls_instruments()))) {
        stop(paste(paste("In exchange_rate(...) : ",
                          "overwrite is FALSE and primary_id", 
                          if (sum(in.use) > 1) "s are" else " is", 
                          " already in use:\n", sep=""),
                   paste(intersect(primary_id, li), collapse=", ")), 
             call.=FALSE)
  }

  if (length(primary_id) > 1) {
    out <- sapply(primary_id, exchange_rate, identifiers=identifiers, 
                         assign_i=assign_i, ...=..., simplify=assign_i)
    return(if (assign_i) unname(out) else out)
  }
  if (is.null(currency)) currency <- substr(primary_id,4,6)
  if (is.null(counter_currency)) counter_currency <- substr(primary_id,1,3)
  if(!exists(currency, where=.instrument,inherits=TRUE)) {
    warning(paste("currency",currency,"not found")) # assumes that we know where to look
  }
  if(!exists(counter_currency, 
    where=.instrument,inherits=TRUE)) {
        warning(paste("counter_currency",counter_currency,"not found")) # assumes that we know where to look
  }

  ## now structure and return
  instrument(primary_id=primary_id , currency=currency , multiplier=1, 
             tick_size=tick_size, identifiers = identifiers, ..., 
             counter_currency=counter_currency, 
             type=c("exchange_rate","currency"), assign_i=assign_i)
}

#TODO  auction dates, coupons, etc for govmt. bonds
#' @export
#' @rdname instrument
bond <- function(primary_id, currency, multiplier, tick_size=NULL, 
                 identifiers = NULL, assign_i=TRUE, overwrite=TRUE, ...){
    if (missing(currency)) stop ("'currency' is a required argument")
    if (length(primary_id) > 1) stop("'primary_id' must be of length 1 for this function")
    if (!isTRUE(overwrite) && isTRUE(assign_i) && primary_id %in% ls_instruments()) {
        stop("overwrite is FALSE and the primary_id ", sQuote(primary_id), 
             " is already in use.")
    }
    instrument(primary_id=primary_id, currency=currency, multiplier=multiplier, 
               tick_size=tick_size, identifiers = identifiers, ..., type="bond", 
               assign_i=assign_i )
}

#' @export
#' @rdname series_instrument
bond_series <- function(primary_id , suffix_id, ..., first_traded=NULL, 
                        maturity=NULL, identifiers = NULL, 
                        payment_schedule=NULL, assign_i=TRUE){
    contract<-try(getInstrument(primary_id))
    if(!inherits(contract,"bond")) {
        stop("bonds contract spec must be defined first")
    }
    
    # TODO add check for Date equivalent in first_traded and expires
    
    ## with bond series we probably need to be more sophisticated,
    ## and find the existing series from prior periods (probably years or months)
    ## and then add the first_traded and expires to the time series by splicing
    id<-paste(primary_id, suffix_id,sep="_")
    temp_series<-try(getInstrument(id),silent=TRUE)
    if(inherits(temp_series,"bond_series")) {
        message("updating existing first_traded and maturity for ",id)
        temp_series$first_traded<-c(temp_series$first_traded,first_traded)
        temp_series$maturity<-c(temp_series$maturity,maturity)
        assign(id, temp_series, 
               envir=as.environment(.instrument))
    } else {
        dargs<-list(...)
        dargs$currency=NULL
        dargs$multiplier=NULL
        dargs$type=NULL
        temp_series = instrument( primary_id = id,
                suffix_id=suffix_id,
                currency = contract$currency,
                multiplier = contract$multiplier,
                tick_size=contract$tick_size,
                first_traded = first_traded,
                maturity = maturity,
                identifiers = identifiers,
                type=c("bond_series", "bond"),
                ...=dargs,
                assign_i=assign_i
        ) 
    }
}

#' Create an instrument based on name alone
#'
#' Given a name, this function will attempt to create
#' an instrument of the appropriate type.
#'
#' If \code{currency} is not already defined, it will be defined (unless it is 
#' not 3 uppercase characters).  The default value for \code{currency} is 
#' \dQuote{USD}.  If you do not provide a value for \code{currency}, 
#' \dQuote{USD} will be defined and used to create the instrument.  
#'
#' If \code{primary_id} is 6 uppercase letters and \code{default_type} is not 
#' provided, it will be assumed that it is the primary_id of an 
#' \code{\link{exchange_rate}}, in which case, the 1st and 2nd half of 
#' \code{primary_id} will be defined as \code{\link{currency}}s if not 
#' the names of already defined \code{\link{instrument}}s.
#' If the \code{primary_id} begins with a \dQuote{^} it will be assumed that it 
#' is a yahoo symbol and that the instrument is an index (synthetic), and the 
#' \sQuote{src} will be set to \dQuote{yahoo}. 
#' (see \code{\link{setSymbolLookup}})
#'
#' If it is not clear from the \code{primary_id} what type of instrument to 
#' create, an instrument of type \code{default_type} will be created (which is 
#' 'NULL' by default).  This will happen when \code{primary_id} is that of a 
#' \code{\link{stock}}, \code{\link{future}}, \code{\link{option}}, or 
#' \code{\link{bond}}.  This may also happen if \code{primary_id} is that of a 
#' \code{\link{future_series}} or \code{\link{option_series}} but the 
#' corresponding \code{future} or \code{option} cannot be found.  In this case, 
#' the instrument type would be \code{default_type}, but a lot of things would 
#' be filled in as if it were a valid series instrument (e.g. \sQuote{expires}, 
#' \sQuote{strike}, \sQuote{suffix_id}, etc.)
#' @param primary_id charater primary identifier of instrument to be created
#' @param currency character name of currency that instrument will be 
#'   denominated it. Default=\dQuote{USD}
#' @param multiplier numeric product multiplier
#' @param silent TRUE/FALSE. silence warnings?
#' @param default_type What type of instrument to make if it is not clear from 
#'   the primary_id. ("stock", "future", etc.) Default is NULL.
#' @param root character string to pass to \code{\link{parse_id}} to be used as 
#'   the root_id for easier/more accurate parsing.
#' @param assign_i TRUE/FALSE. Should the \code{instrument} be assigned in the 
#'   \code{.instrument} environment?
#' @param ... other passthrough parameters
#' @return Primarily called for its side-effect, but will return the name of the 
#'   instrument that was created
#' @note This is not intended to be used to create instruments of type 
#'   \code{stock}, \code{future}, \code{option},
#' or \code{bond} although it may be updated in the future.
#' @author Garrett See
#' @examples
#' \dontrun{
#' instrument.auto("CL_H1.U1")
#' getInstrument("CL_H1.U1") #guaranteed_spread
#' 
#' instrument.auto("ES_H1.YM_H1")
#' getInstrument("ES_H1.YM_H1") #synthetic
#' 
#' currency(c("USD","EUR"))
#' instrument.auto("EURUSD")
#' getInstrument("EURUSD") #made an exchange_rate
#' 
#' instrument.auto("VX_H11") #no root future defined yet!
#' getInstrument("VX_H11") #couldn't find future, didnt make future_series
#' future("VX","USD",1000,underlying_id=synthetic("SPX","USD")) #make the root 
#' instrument.auto("VX_H11") #and try again
#' getInstrument("VX_H11") #made a future_series
#' }
#' @export
instrument.auto <- function(primary_id, currency=NULL, multiplier=1, silent=FALSE, 
                            default_type='unknown', root=NULL, assign_i=TRUE, 
                            ...) {
##TODO: check formals against dots and remove duplicates from dots before calling constructors to avoid
# 'formal argument "multiplier" matched by multiple actual arguments'
    if (!is.null(currency) && !is.currency.name(currency)) {
        if (nchar(currency) != 3 || currency != toupper(currency)) {
            stop(paste(currency, "is not defined, and it will not be auto",
                       "defined because it does not appear to be valid."))
        }
        currency(currency)
        if (!silent) message(paste('Created currency', currency,
                                   'because it was not defined.\n'))
    } 
    warned <- FALSE
    dargs <- list(...)
    primary_id <- make.names(primary_id)
    pid <- parse_id(primary_id, root=root)
    type <- NULL
    if (any(pid$type == 'calendar')) {
        return(guaranteed_spread(primary_id, currency=currency, 
                                defined.by='auto', multiplier=multiplier, 
                                assign_i=assign_i, ...))
    } 
    if (any(pid$type == 'butterfly')) {
        return(butterfly(primary_id, currency=currency, defined.by='auto', 
                         assign_i=assign_i, ...))
    }
    if (any(pid$type == 'ICS')) {
        root <- getInstrument(pid$root, type='ICS_root', silent=TRUE)
        if (is.instrument(root)) {
            return(ICS(primary_id, assign_i=assign_i, ...))
        } else {
            #TODO: look for members in dots
            if (!silent) {
                warning(paste(primary_id, " appears to be an ICS, ", 
                        "but its ICS_root cannot be found. ",
                        "Creating _", default_type, "_ instrument instead.", 
                        sep=""))
                warned <- TRUE
            }
            dargs$root_id <- pid$root
            dargs$suffix_id <- pid$suffix
            dargs$expires <- paste(pid$year, 
                                   sprintf("%02d", 
                                           month_cycle2numeric(pid$month)), 
                                   sep="-")
        }
    }
    if (any(pid$type == 'future') || any(pid$type == 'SSF')) {
        root <- getInstrument(pid$root,silent=TRUE,type='future')
        if (is.instrument(root) && !inherits(root, 'future_series')) {
            if (is.null(currency) && is.null(root[['currency']])) {
                if (!isTRUE(silent)) warning('using USD as the currency')
                currency <- currency("USD")
            }
            return(future_series(primary_id, currency=currency, 
                                 defined.by='auto', assign_i=assign_i,...))
        } else {
            if (!isTRUE(silent)) {
                warning(paste(primary_id," appears to be a future_series, ", 
                              "but its root cannot be found. ", 
                              "Creating _", default_type, 
                              "_ instrument instead.", sep=""))
                if (is.null(currency)){
                    warning('using USD as the currency')
                    currency <- currency("USD")
                }
                warned <- TRUE
            }
            dargs$root_id <- pid$root
            dargs$suffix_id <- pid$suffix
            dargs$expires <- paste(pid$year, 
                                   sprintf("%02d", 
                                           month_cycle2numeric(pid$month)), 
                                   sep="-")
        }
    }
    if (any(pid$type == 'option')) {
        root <- getInstrument(pid$root,silent=TRUE,type='option')
        if (is.instrument(root) && !inherits(root, 'option_series')) {
            if (is.null(currency) && is.null(root[['currency']])) {
                if (!isTRUE(silent)) warning('using USD as the currency')
                currency <- currency("USD")
            }
            return(option_series(primary_id, currency=currency, 
                                 defined.by='auto', assign_i=assign_i, ...))
        } else {
            if (!isTRUE(silent)) {
                warning(paste(primary_id," appears to be an option_series, ", 
                             "but its root cannot be found. Creating _", 
                              default_type, "_ instrument instead.", sep=""))
                if (is.null(currency)){
                    warning('using USD as the currency')
                    currency <- currency("USD")
                }
                warned <- TRUE
            }
            dargs$root_id <- pid$root
            dargs$suffix_id <- pid$suffix
            dargs$expires <- if(pid$format == 'opt2') {
                    as.Date(substr(pid$suffix,1,6),format='%y%m%d')
                } else if (pid$format == 'opt4') {
                    as.Date(substr(pid$suffix,1,8),format='%Y%m%d')
                } else paste(pid$year, sprintf("%02d", 
                                               month_cycle2numeric(pid$month)), 
                             sep="-")
            dargs$multiplier=100
            dargs$callput <- switch(pid$right, C='call', P='put')
            dargs$strike <- pid$strike
        }
    } 
    if (any(pid$type == 'exchange_rate'))
        return(exchange_rate(primary_id, defined.by='auto', assign_i=assign_i, 
                             ...))
    #if we weren't given a default_type, then if it's 6 uppercase letters, make an exchange rate    
    if (default_type == 'unknown' && nchar(primary_id) == 6 && 
            sum(attr(gregexpr("[A-Z]",primary_id)[[1]],"match.length")) == 6) {
        if (!is.instrument(getInstrument(substr(primary_id,1,3), silent=TRUE))) {
            ccy.st <- currency(substr(primary_id,1,3), defined.by='auto') 
            if (!silent) { 
                message(paste("Created currency", ccy.st, 
                              "because it was not defined."))
            } 
        }
        if (!is.instrument(getInstrument(substr(primary_id,4,6), silent=TRUE))) { 
            ccy.st <- currency(substr(primary_id,4,6), defined.by='auto') 
            if (!silent) { message(paste("Created currency", ccy.st, 
                                         "because it was not defined.")) 
            }
        }
        return(exchange_rate(primary_id, defined.by='auto', assign_i=assign_i, ...))
    }
    if (any(pid$type == 'synthetic')) {
        if (!is.na(pid$format) && pid$format == 'yahooIndex') {
            if (is.null(currency)) {
                if (!silent) {
                    warning(paste('currency will be assumed to be USD', 
                                  'because NULL is not a currency.'))
                }
                currency <- currency('USD')
            }
            return(synthetic(gsub("\\^","",primary_id), currency=currency, 
                            multiplier=multiplier, 
                            identifiers=list(yahoo=primary_id), 
                            src=list(src='yahoo',name=primary_id),
                            defined_by='auto', assign_i=assign_i, ...))
        } else {
            members <- strsplit(primary_id,"\\.")[[1]]
            if (!all(is.instrument.name(members))) {
                # at least 1 member is not defined, so we have to assume this 
                # is an index (e.g. TICK-NYSE)
                if (is.null(currency)) {
                    if (!silent) {
                        warning(paste('currency will be assumed to be USD',
                                      'because NULL is not a currency.'))
                    }
                    currency <- currency('USD')
                }
                return(synthetic(primary_id, currency=currency, 
                                multiplier=multiplier, defined.by='auto', 
                                assign_i=assign_i, ...))
            }
            # if all members are already defined, we'll pass those to synthetic 
            # and let it figure out currency.
            return(synthetic(members=members, currency=currency, 
                            multiplier=multiplier, defined.by='auto', 
                            assign_i=assign_i, ...) )
        }
    } 
    if (any(pid$type == 'root')) {
        if (primary_id %in% c(paste(pid$root, "O", sep="."), 
                              paste(pid$root, "OQ", sep="."))) { 
            #X.RIC for NASDAQ stock.  e.g. AAPL.O, MSFT.OQ
            return(stock(pid$root, currency=currency('USD'), 
                         multiplier=multiplier, 
                         identifiers=list(X.RIC=primary_id), defined.by='auto', 
                         assign_i=assign_i, ...))
            #update_instruments.yahoo(pid$root) 
        }
    }
    ss <- strsplit(primary_id," ")[[1]]  #take out spaces (OSI uses spaces, but makenames would turn them into dots)
    ss <- ss[!ss %in% ""]
    if (length(ss) == 2) primary_id <- paste(ss,collapse="_")
    dargs$primary_id <- primary_id
    dargs$currency <- if (!is.null(currency)) { 
        currency 
    } else { 
        if (!silent) warning('using USD as the currency'); currency("USD") 
    } 
    dargs$multiplier <- multiplier
    dargs$defined.by='auto'
    dargs$assign_i <- assign_i
    if(is.function(try(match.fun(default_type),silent=TRUE))) {
        if (!silent && !warned) 
            warning('Creating a _', default_type, '_ instrument because ', 
                    primary_id, ' is of an ambiguous format.') 
         return(do.call(default_type, dargs))
    }
    if (!silent && !warned) 
        warning(paste(primary_id, 'is not of an unambiguous format.', 
                'Creating _unknown_ instrument with multiplier 1.'))
    dargs$type <- 'unknown'
    do.call(instrument, dargs)
}

  
#' Primary accessor function for getting objects of class 'instrument'
#' 
#' This function will search the \code{.instrument} environment for objects of
#' class \code{type}, using first the \code{primary_id} and then any 
#' \code{identifiers} to locate the instrument.  Finally, it will try adding 1 
#' and then 2 dots to the beginning of the \code{primary_id} to see if an 
#' instrument was stored there to avoid naming conflicts.
#' 
#' \code{\link{future}} and \code{\link{option}} objects may have a primary_id 
#' that begins with 1 or 2 dots (in order to avoid naming conflics).  For 
#' example, the root specs for options (or futures) on the stock with ticker 
#' "SPY" may be stored with a primary_id of "SPY", ".SPY", or "..SPY".  
#' \code{getInstrument} will try using each possible \code{primary_id}
#' until it finds an instrument of the appropriate \code{type}
#' @param x String identifier of instrument to retrieve
#' @param Dates date range to retrieve 'as of', may not currently be implemented
#' @param silent if TRUE, will not warn on failure, default FALSE
#' @param type class of object to look for. See Details
#' @examples
#' \dontrun{
#' option('..VX', multiplier=100, 
#'   underlying_id=future('.VX',multiplier=1000, 
#'     underlying_id=synthetic('VIX', currency("USD"))))
#'
#' getInstrument("VIX")
#' getInstrument('VX') #returns the future
#' getInstrument("VX",type='option')
#' getInstrument('..VX') #finds the option
#' }
#' @export
#' @rdname getInstrument
getInstrument <- function(x, Dates=NULL, silent=FALSE, type='instrument'){
    tmp_instr <- try(get(x,pos=.instrument),silent=TRUE)
    if(inherits(tmp_instr,"try-error") || !inherits(tmp_instr, type)){
        xx <- make.names(x)
        ## First, look to see if x matches any identifiers.
        # unlist all instruments into a big named vector
        ul.instr <- unlist(as.list(.instrument, 
                                   all.names=TRUE))
        # subset by names that include "identifiers"
        ul.ident <- ul.instr[grep('identifiers', names(ul.instr))]
        # if x (or make.names(x)) is in the identifiers subset, extract the 
        # primary_id from the name
        tmpname <- ul.ident[ul.ident %in% unique(c(x, xx))]
        # if x was not in ul.ident, tmpname will == named character(0)
        if (length(tmpname) > 0) {
            #primary_id is everything before .identifiers
            id <- gsub("\\.identifiers.*", "", names(tmpname))
            tmp_instr <- try(get(id, pos=.instrument), 
                             silent=TRUE)
            if (inherits(tmp_instr, type)) {
                #&& (x %in% tmp_instr$identifiers || x %in% make.names(tmp_instr$identifiers))
                return(tmp_instr) 
            }
        }
        #If not found, see if it begins with dots (future or option root)
        #Remove any dots at beginning of string and add them back 1 at a time 
        # to the beginning of id.
        char.x <- strsplit(x, "")[[1]] # split x into vector of characters
        x <- substr(x, grep("[^\\.]", char.x)[1], length(char.x)) # excluding leading dots
        tmp_instr<-try(get(x,pos=.instrument),silent=TRUE)
        if(!inherits(tmp_instr,type)) {
            tmp_instr<-try(get(paste(".",x,sep=""),
                               pos=.instrument),
                           silent=TRUE)
            if(!inherits(tmp_instr,type)) {
                tmp_instr<-try(get(paste("..",x,sep=""),
                                   pos=.instrument),
                               silent=TRUE)
            }
        }
        if (inherits(tmp_instr, type)) return(tmp_instr)
        if(!silent) warning(paste(type,x,"not found, please create it first."))
        return(FALSE)
    } else{
        return(tmp_instr)
    }
    #TODO add Date support to instrument, to get the proper value given a specific date
}


#' Add or change an attribute of an instrument
#' 
#' This function will add or overwrite the data stored in the specified slot of 
#' the specified instrument.
#'
#' If the \code{attr} you are trying to change is the \dQuote{primary_id,} the 
#' instrument will be renamed. (A copy of the instrument will be stored by the 
#' name of \code{value} and the old instrument will be removed.)
#' If the \code{attr} you are changing is \dQuote{type}, the instrument will be 
#' reclassed with that type. If \code{attr} is \dQuote{src}, \code{value} will 
#' be used in a call to \code{setSymbolLookup}.  Other checks are in place to 
#' make sure that \dQuote{currency} remains a \code{\link{currency}} object and 
#' that \dQuote{multiplier} and \dQuote{tick_size} can only be changed to 
#' reasonable values.
#' 
#' If \code{attr} is \dQuote{identifiers} and \code{value} is \code{NULL}, 
#' \code{identifiers} will be set to \code{list()}.  If \code{value} is not a 
#' list, \code{\link{add.identifier}} will be called with \code{value}.
#' \code{add.identifier} will convert \code{value} to a list and append it to
#' the current \code{identifiers}
#' @param primary_id primary_id of the instrument that will be updated
#' @param attr Name of the slot that will be added or changed
#' @param value What to assign to the \code{attr} slot of the \code{primary_id} 
#'   instrument
#' @param ... arguments to pass to \code{getInstrument}. For example,
#'   \code{type} could be provided to allow for \code{primary_id} to be an
#'   identifier that is shared by more that one instrument (of different types)
#' @return called for side-effect
#' @note You can remove an attribute/level from an instrument by calling this 
#'   function with \code{value=NULL}
#' @examples
#' \dontrun{
#' currency("USD")
#' stock("SPY","USD")
#' instrument_attr("USD","description","U.S. Dollar")
#' instrument_attr("SPY", "description", "An ETF")
#' getInstrument("USD")
#' getInstrument("SPY")
#' 
#' #Call with value=NULL to remove an attribute
#' instrument_attr("SPY", "description", NULL)
#' getInstrument("SPY")
#'
#' instrument_attr("SPY","primary_id","SPX") #move/rename it
#' instrument_attr("SPX","type","synthetic") #re-class
#' instrument_attr("SPX","src",list(src='yahoo',name='^GSPC')) #setSymbolLookup
#' getSymbols("SPX") #knows where to look because the last line setSymbolLookup
#' getInstrument("SPX")
#' }
#' @export
instrument_attr <- function(primary_id, attr, value, ...) {
    instr <- try(getInstrument(primary_id, silent=TRUE, ...))
    if (inherits(instr, 'try-error') || !is.instrument(instr))
        stop(paste('instrument ',primary_id,' must be defined first.',sep=''))
    if (attr == 'primary_id') {
        rm(list = primary_id, pos = .instrument)
    } else if (attr == 'currency') {
        if (!is.currency.name(value)) {
            stop("currency ", value, " must be an object of type 'currency'")
        }
    } else if (attr == 'multiplier') {
        if (!is.numeric(value) || length(value) > 1) {
            stop("multiplier must be a single number")
        }
    } else if (attr == 'tick_size') {
        if (!is.null(value) && (!is.numeric(value) || length(value) > 1)) {
            stop("tick_size must be NULL or a single number")
        }
    } else if (attr == 'type') {
        tclass <- unique(c(value, "instrument"))
        class(instr) <- tclass
    } else if (attr == 'IB') {
        if (inherits(value, 'twsContract')) {
            class(instr) <- unique(c(class(instr)[1], 'twsInstrument', 
                                     class(instr)[-1]))
        } else {
            warning('non-twsContract assigned to $IB')
            class(instr) <- class(instr)[!class(instr) %in% 'twsInstrument']
        }
    } else if (attr == 'src') {
        sarg <- list()
        sarg[[instr$primary_id]] <- value
        setSymbolLookup(sarg)
    } else if (attr == 'identifiers') {
        if (length(value) == 0) {
            value <- list()
        } else if (!is.list(value)) {
            #warning("identifiers must be a list. Appending current identifiers.")
            # add.identifier will convert to list
            return(add.identifier(primary_id, value))
        }
    }
    instr[[attr]] <- value
    assign(instr$primary_id, instr, pos=.instrument)
}


#' Add an identifier to an \code{instrument}
#' 
#' Add an identifier to an \code{\link{instrument}} unless the instrument 
#' already has that identifier.
#' @param primary_id primary_id of an \code{\link{instrument}}
#' @param ... identifiers passed as regular named arguments.
#' @return called for side-effect
#' @author Garrett See
#' @seealso \code{\link{instrument_attr}}
#' @examples
#' \dontrun{
#' stock("XXX", currency("USD"))
#' add.identifier("XXX", yahoo="^XXX") 
#' getInstrument("^XXX")
#' add.identifier("^XXX", "x3")
#' all.equal(getInstrument("x3"), getInstrument("XXX")) #TRUE
#' }
#' @export
add.identifier <- function(primary_id, ...) {
    new.ids <- as.list(unlist(list(...)))
    instr <- getInstrument(primary_id)
    if (!inherits(instr, "instrument")) {
        stop(paste(primary_id, "is not a defined instrument"))
    }
    ids <- c(instr[["identifiers"]], new.ids)
    if (all(is.null(names(ids)))) {
        instrument_attr(primary_id, "identifiers", unique(ids))
    } else instrument_attr(primary_id, "identifiers", 
                ids[!(duplicated(unlist(ids)) & duplicated(names(ids)))])
}


#' Add a source to the defined.by field of an \code{instrument}
#' 
#' Concatenate a string or strings (passed through dots) to the defined.by 
#' field of an instrument (separated by semi-colons).  Any duplicates will be
#' removed.  See Details.
#' 
#' If there is already a value for the \code{defined.by} attribute of the 
#' \code{primary_id} instrument, that string will be split on semi-colons and
#' converted to a character vector.  That will be \code{c}ombined with any new 
#' strings (in \code{...}).  The unique value of this new vector will then
#' be converted into a semi-colon delimited string that will be assigned to 
#' the \code{defined.by} attribute of the \code{primary_ids}' instruments
#' 
#' Many functions that create or update instrument definitions will also add or
#' update the value of the defined.by attribute of that instrument.  If an 
#' instrument has been updated by more than one function, it's \code{defined.by}
#' attribute will likely be a semi-colon delimited string (e.g. 
#' \dQuote{TTR;yahoo}).  
#' @param primary_ids character vector of primary_ids of 
#'   \code{\link{instrument}}s
#' @param ... strings, or character vector, or semi-colon delimited string.
#' @return called for side-effect
#' @author Garrett See
#' @seealso \code{\link{add.identifier}}, \code{\link{instrument_attr}}
#' @examples
#' \dontrun{
#' update_instruments.TTR("GS")
#' getInstrument("GS")$defined.by #TTR
#' add.defined.by("GS", "gsee", "demo")
#' add.defined.by("GS", "gsee;demo") #same
#' }
#' @export
add.defined.by <- function(primary_ids, ...) {
    for(id in primary_ids) {
        db <- getInstrument(id)[["defined.by"]]
        instrument_attr(id, "defined.by", 
            paste(unique(c(unlist(strsplit(db, ";")), 
                unlist(strsplit(unlist(list(...)), ";")))), collapse=";"))
    }
}


#' instrument class print method
#' 
#' @author Joshua Ulrich, Garrett See
#' @keywords internal
#' @export
print.instrument <- function(x, ...) {
  str(unclass(x), comp.str="", no.list=TRUE, give.head=FALSE,
    give.length=FALSE, give.attr=FALSE, nest.lev=-1, indent.str="")
  invisible(x)
}

#' instrument class sort method
#' 
#' @author Garrett See
#' @keywords internal
#' @export
sort.instrument <- function(x, decreasing=FALSE, na.last=NA, ...) {
    anchored <- x[c("primary_id", "currency", "multiplier", "tick_size", 
                  "identifiers", "type")] 
    sortable <- x[setdiff(names(x), names(anchored))]
    out <- c(anchored, sortable[order(names(sortable), decreasing=decreasing, 
                                      na.last=na.last, ...)])
    class(out) <- class(x)
    out
}
