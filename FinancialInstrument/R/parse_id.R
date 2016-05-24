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
# $Id: parse_id.R 1655 2014-11-23 22:53:26Z gsee $
#
###############################################################################


#' Parse a primary_id
#'
#' Extract/infer descriptive information about an instrument from its name.
#'
#' This function is primarily intended to be used on the names of \code{\link{future_series}} 
#' and \code{\link{option_series}} instruments, and it will work best if the id has an 
#' underscore in it that separates the root_id from the suffix_id.  (However, it should be able 
#' to handle most ids even if the underscore is missing). 
#' After splitting \code{x} into a root_id and suffix_id, the suffix_id is 
#' passed to \code{\link{parse_suffix}} (see also) for further processing.
#' 
#' TODO: add support for bond_series.
#' @param x the id to be parsed (e.g. \sQuote{ES_U11}, \sQuote{SPY_111217C130})
#' @param silent silence warnings?
#' @param root character name of instrument root_id.  Optionally provide this to make parsing easier.
#' @return a list of class \sQuote{id.list} containing \sQuote{root} and \sQuote{suffix} as well as 
#' what is returned from \code{\link{parse_suffix}} (type, month, year, strike, right, cm, cc, format)
#' @author Garrett See
#' @note this function will identify \code{x} as an \code{\link{exchange_rate}} only if it is 
#' 6 characters long and made up of 2 previously defined \code{\link{currency}} instruments.
#' @seealso \code{\link{parse_suffix}}
#' @examples
#' parse_id("ES_Z11")
#' parse_id("CLZ1")
#' parse_id("SPY_111217C130")
#' @export
parse_id <- function(x, silent=TRUE, root=NULL) {
    sufftype <- TRUE #will we use the type given by parse_suffix, or overwrite it with e.g. 'exchange_rate', or 'synthetic'
    suffformat <- TRUE #If x begins with "^" this will be set to FALSE, and we'll overwrite parse_suffix(...)$format with yahooIndex"
    x <- gsub("-", ".", x)
    all.numeric <- as.logical(!is.na(suppressWarnings(as.numeric(x))))
    if (!is.null(root)) {
        suffix <- sub(root,"",x) #turns ESU1 into U1, or ES_U11 into _U11 
        suffix <- gsub("_","",suffix) #take out the underscore if there is one
    } else if (identical(integer(0), grep("[0-9]",x))) { 
        #if there are no numbers in the id, then it has no year, so it is not a recognized future or option        
        if (substr(x, nchar(x)-1, nchar(x)) %in% c(".O", ".K")) {
            #only one dot, and the id ends with ".O" -- it looks like an X.RIC for a NASDAQ stock
            suffix <- ""
            root <- substr(x, 1, nchar(x)-2)
            type <- 'root'
        } else if (substr(x, nchar(x)-2, nchar(x)) == ".OQ") {
            # id ends with ".OQ" -- looks like an X.RIC for NMS stock
            suffix <- ""
            root <- substr(x, 1, nchar(x)-3)
            type <- 'root'
        } else if (identical(all.equal(nchar(x) - nchar( gsub("\\.","",x)),1), TRUE)) { 
            #only 1 dot, so it's not a fly
            #SPY.DIA, EUR.USD, SPY110917C122.5, T2010917P25
            if (suppressWarnings(!is.na(as.numeric(strsplit(x,"\\.")[[1]][2])))) { #probably an option with a decimal in the strike
                #if we take out all the numbers, periods, and dashes, 
                #we should be left with the ticker and either "C" or "P"                
                root <- gsub("[0-9.-]","",x) #now it looks like SPYC or TP
                root <- substr(root, 1,nchar(root)-1)
                suffix <- gsub(root,"",x) #whatever isn't the root
            } else { #probably a synthetic: SPY.DIA, GLD.EUR
                suffix <- ""
                root <- x
                type <- 'synthetic'
                sufftype <- FALSE
            }
        } else if (identical(all.equal(nchar(x) - nchar( gsub("\\.","",x)),2), TRUE)) { 
            #2 dots, so we'll treat it as a fly, although it could be a basket
            #SPY.DIA.QQQ, 
            suffix <- ""
            root <- x
            type <- 'synthetic'
            sufftype <- FALSE
        } else { #no dots, no numbers
            root <- x
            suffix <- ""
            if (substr(x,1,1) == "^") {
                type <- c('synthetic','root')
                format <- 'yahooIndex'
                sufftype <- FALSE
                suffformat <- FALSE
            } else if (nchar(x) == 6) {
                if (is.instrument(getInstrument(substr(x,1,3),silent=TRUE)) 
                    && is.instrument(getInstrument(substr(x,4,6),silent=TRUE))) {
                    type <- c('exchange_rate', 'root')
                    sufftype <- FALSE
                }
            }
        }
    } else if (!all.numeric && identical(x, gsub('_','',x))) { #no underscore; have to guess what is root and what is suffix       
        hasdot <- !identical(integer(0),grep("\\.",x))        
        if (!silent && !hasdot) 
            warning("id of future_series should have an underscore in it. Trying to parse anyway.")
        #if (nchar(x) < 9 && !hasdot) { #assume it's a future like ESU1 or ESU11
        if (!hasdot) { #assume it's a future like ESU1 or ESU11
            if (suppressWarnings(!is.null(parse_suffix(substr(x,3,nchar(x)))) && 
                    !is.na(parse_suffix(substr(x,3,nchar(x)))$format))) {
                root <- substr(x,1,2)
                suffix <- substr(x,3,nchar(x))
            } else if (suppressWarnings(!is.null(parse_suffix(substr(x,4,nchar(x)))) && 
                    !is.na(parse_suffix(substr(x,4,nchar(x)))$format))) {
                root <- substr(x,1,3)        
                suffix <- substr(x,4,nchar(x))
            } else if (suppressWarnings(!is.null(parse_suffix(substr(x,5,nchar(x)))) && 
                    !is.na(parse_suffix(substr(x,5,nchar(x)))$format))) {
                root <- substr(x,1,4)        
                suffix <- substr(x,5,nchar(x))
            } else if (suppressWarnings(!is.null(parse_suffix(substr(x,5,nchar(x)))) && 
                    !is.na(parse_suffix(substr(x,6,nchar(x)))$format))) {
                root <- substr(x,1,5)        
                suffix <- substr(x,6,nchar(x))
            } else if (suppressWarnings(!is.null(parse_suffix(substr(x,2,nchar(x)))) && 
                    !is.na(parse_suffix(substr(x,2,nchar(x)))$format))) {
                root <- substr(x,1,1)
                suffix <- substr(x,2,nchar(x))
            } else if (nchar(x) <= 3) {
                root <- substr(x, 1, nchar(x))
                suffix <- "" 
            } else { 
                root <- substr(x,1,6)
                suffix <- substr(x,7,nchar(x))
            }
        } else {
            # Look for spread (e.g. CLF2.F3). Split by the dot, then try to extract suffix_id from the part before
            # the dot. If that suffix is the same length as the part after the dot, we'll assume it's a spread
            # (it could also be a fly or condor, but only the first 2 parts are checked)
            #x <- "CLF2.G2"
            ss <- strsplit(x, "\\.")[[1]]
            if (length(ss) > 1 && !any(ss == "")) {
                pidhalf <- parse_id(ss[1]) #parse_id on part before dot (looking for spread)
                if (nchar(pidhalf$suffix) == nchar(strsplit(x, "\\.")[[1]][2])) {
                    root <- pidhalf$root
                    suffix <- gsub(root, "", x)
                } else if ("outright" %in% pidhalf$type) {
                    #x <- "DIA111230P139.75"
                    root <- gsub("[0-9.-]","",x) #now it looks like DIAP, SPYC or TP
                    root <- substr(root, 1,nchar(root)-1)
                    suffix <- gsub(root,"",x) #whatever isn't the root
                } else {
                    #x <- "USD.1" #a stock whose symbols was created by 'make.names'
                    suffix <- ""
                    root <- x
                    type <- 'root'
                }
            } else {
                #has a dot, has a number, and a non-number, and no underscore.  Probably a strange root (e.g. "..BL2")
                root <- x
                suffix <- ""
                type <- 'root'
                format <- NA
                sufftype <- FALSE
                suffformat <- FALSE
            } 
        }
    } else { #there _is_ an underscore and at least 1 number.     
        #if there are dots then maybe it is a spread of futures?
        #e.g. "CL_N1.HO_M1"
        ss <- strsplit(x,"\\.")[[1]]
        has.und <- function(x) { #TRUE if it has an underscore
            sapply(x, FUN=function(x) !identical(x, gsub('_','',x)))
        }
        if (length(ss) > 1) {    
            if (all(has.und(ss))) { #all parts have an underscore. 
                #e.g. CL_N1.HO_M1 --> "CL_N1" "HO_M1"
                #or CL_N1.CL_M1           
                ssu <- strsplit(ss,"_")
                tmprt <- ssu[[1]][1]
                if (all(sapply(ssu, FUN=function(x) x[1] == tmprt))) {
                    #all share a root_id
                    return(parse_id(make_spread_id(ss)))
                } else {
                    #all parts have an underscore, but they don't share a root.
                    #e.g. "CL_N1.HO_M1"
                    root <- x
                    suffix <- ""
                    type <- 'synthetic'
                    sufftype <- FALSE
                }  
            } else if (has.und(ss[1]) && !has.und(ss[2])) {
                #First part has underscore, but second doesn't. e.g. CL_H1.M1 --> "CL_H1" "M1"
                spl.und <- strsplit(x,"_")[[1]]
                root <- spl.und[1]
                suffix <- paste(spl.und[2:length(spl.und)],collapse=".")
                    
            } else { 
                suffix <- ""
                root <- x
                type <- 'synthetic'
                sufftype <- FALSE
            } 
        } else { #ES_U1, ES_U1M1, 
            root <- strsplit(x,"_")[[1]][1]
            suffix <- strsplit(x,"_")[[1]][2]
        }        
    }
    root <- gsub(" ","",root) 
    if (!is.na(suffix)) {
        suff <- parse_suffix(suffix, silent=silent)
        if (sufftype) type <- suff$type
        if (suffformat) format <- suff$format
    } else {
        suff <- NULL
        type <- "Asian"
        format <- "numeric"
    }
    structure(list(root=root, suffix=suffix, type=type, month=suff$month, 
                    year=suff$year, strike=suff$strike, right=suff$right, 
                    cm=suff$cm, cc=suff$cc, format=format),class='id.list')
}

#' parse a suffix_id
#' 
#' extract information from the suffix_id of an instrument
#'
#' These would be recognized as a Sep 2011 outright futures contract: 
#' U1, U11, SEP1, SEP11, U2011, Sep2011, SEP2011
#' 
#' These would be recognized as a call with a strike of 122.5 that expires Sep 17, 2011:
#' 110917C122.5, 20110917C122.5, 110917C00122500, 20110917C00122500
#'
#' These would be recognized as Sep 2011 single stock futures:
#' 1CU1, 1CU11, 1CSEP11, 1DU1 (dividend protected)
#'
#' These would be recognized as Adjusted futures:
#' cm.30 (30 day constant maturity future),
#' cc.OI (continuous contract rolled when Open Interest rolls),
#' cc.Vol (continuous contract roll when Volumn rolls),
#' cc.Exp.1 (continuous contract rolled 1 day before Expiration)
#'
#' Synthetics and spreads:
#'
#' SPY.DIA --> type == synthetic; 
#'
#' U1.Z1 or U11.Z11 --> type == "calendar", "spread"; month == 'SEP', year == 2011
#'
#' U1.0302 --> type == "ICS", "spread"; month == 'SEP', year == 2011
#'
#'
#' 110917C125.110917P125 --> type == option_spread, spread
#' @param x the suffix_id to be parsed
#' @param silent silence warnings? (warning will usually be about inferring a 4 digit year from a 1 or 2 digit year)
#' @return an object of class \sQuote{suffix.list} which is a list containing \sQuote{type} of instrument, 
#' \sQuote{month} of expiration, \sQuote{year} of expiration, \sQuote{strike} price of option, 
#' \sQuote{right} of option (\dQuote{C} or \dQuote{P}), \sQuote{cm} (maturity in days of a constant maturity contract),
#' \sQuote{cc} (method for calculating a continuous contract), \sQuote{format} (string that indicates the format of the unparsed id).
#' @author Garrett See
#' @seealso \code{\link{parse_id}}, \code{\link{format_id}}
#' @examples
#' parse_suffix("U11")
#' parse_suffix("110917C125")
#' @export
parse_suffix <- function(x, silent=TRUE) {
#TODO better support for spreads and flies; inter and intra. Maybe add a 'members' to suffix.list
    type <- 'outright'
    cm <- FALSE
    cc <- FALSE
    month <- 0
    year <- 0
    strike <- NA
    right <- NA
    format <- NA
    if (x == "") {
        type <- "root"
    } else if (!identical(gsub("cm.","",x), x)) {
    #A 30 day constant maturity synthetic futures contract 
    #on the vix would look like VX_cm.30 
        type <- c('outright', 'cm')
        cm <- as.numeric(strsplit(x,"\\.")[[1]][2])
        format <- 'cm'
    } else if (!identical(gsub("cc.","",x),x)) {
    #cc.OI for rolling on Open Interest, 
    #cc.Vol for rolling on Volume, 
    #cc.Exp.1 for rolling 1 day before Expiration date.  (Exp.0 would be rolling on expiration)    
        type <- c('outright', 'cc')
        cc <- gsub('cc.','',x)
        format <- 'cc'
    } else if (nchar(x) > 7 && (any(substr(x,7,7) == c("C","P")) || any(substr(x,9,9) == c("C","P"))) ) {
        # if the 7th or 9th char is a "C" or "P", it's an option
        # 110917C125 or 20110917C125 or 110917C00125000 or 20110917C00125000
        hasdot <- !identical(integer(0),grep("\\.",x))
        if (!hasdot
            || (hasdot 
                && !is.na(as.numeric(strsplit(x,"\\.")[[1]][2])))) {
                #&& nchar(strsplit(x,"\\.")[[1]][2]) <= 2)) {
            #if it doesn't have a dot, or it does have dot, but what follows 
            #the dot is numeric, then it's an option outright
            if (any(substr(x,7,7) == c("C","P"))) {
                type <- c("outright","option")
                month <- toupper(month.abb[as.numeric(substr(x,3,4))])
                year <- 2000 + as.numeric(substr(x,1,2))
                strike <- as.numeric(substr(x,8,nchar(x)))
                if (nchar(x) >= 15) strike <- strike/1000                 
                right <- substr(x,7,7)
                format <- 'opt2'    
            } else if (any(substr(x,9,9) == c("C","P"))) {
                type <- c("outright","option")
                month <- toupper(month.abb[as.numeric(substr(x,5,6))])
                year <- as.numeric(substr(x,1,4))
                strike <- as.numeric(substr(x,10,nchar(x)))
                if (nchar(x) >= 15) strike <- strike/1000
                right <- substr(x,9,9)
                format <- 'opt4'
            } else stop("how did you get here?")
        } else type <- c("option_spread","spread")
    } else if (!identical(gsub("\\.","",x),x)) { #has a dot. U1.Z1, U11.Z11, SPY.DIA, 
        if (identical(all.equal(nchar(x) - nchar( gsub("\\.","",x)),1), TRUE)) { #only 1 dot, so it's not a fly
            #U1.Z1, U11.Z11, SPY.DIA, EUR.USD, H2.0302
            s <- strsplit(x,"\\.")[[1]]
            s1 <- try(parse_suffix(s[1],silent=TRUE),silent=TRUE)
            s2 <- try(parse_suffix(s[2],silent=TRUE),silent=TRUE)
            if (inherits(s1,'try-error')) {
                s1 <- parse_id(s[1],silent=TRUE)
            }
            if (inherits(s2,'try-error')) {
                s2 <- parse_id(s[2],silent=TRUE)
            }
            
            if (all(c(s1$type,s2$type) == 'root')) {
                type='synthetic'
            } else {
                type <- if (!is.na(s1$format) 
                         && !is.na(s2$format)
                         && !s1$format %in% c("opt2", "opt4", "NNNN")
                         && (s2$format == "NNNN")) {
                         #H2.0302 (e.g. suffix of March 12 FYT ICS)
                            c("ICS", "spread")
                        } else c("calendar", "spread")
                month <- s1$month
                year <- s1$year                
                format <- paste(s1$format, s2$format, sep=".")
            }
        } else if (identical(all.equal(nchar(x) - nchar(gsub("\\.","",x)),2), TRUE)) { #2 dots; it's a fly
            #U1.Z1.H2, U11.Z11.H12, SPY.DIA.QQQ
            s <- strsplit(x, "\\.")[[1]]
            s1 <- try(parse_suffix(s[1],silent=TRUE),silent=TRUE)
            s2 <- try(parse_suffix(s[2],silent=TRUE),silent=TRUE)
            s3 <- try(parse_suffix(s[3],silent=TRUE),silent=TRUE)
            if (inherits(s1,'try-error')) s1 <- parse_id(s[1],silent=TRUE)
            if (inherits(s2,'try-error')) s2 <- parse_id(s[2],silent=TRUE)
            if (inherits(s3,'try-error')) s3 <- parse_id(s[3],silent=TRUE)          
            if (all(c(s1$type,s2$type,s3$type) == 'root')) {
                type='synthetic' #don't really know if it's fly-like or a basket
            } else { 
                type=c('butterfly','spread')
            }
            format <- unique(c(s1$format, s2$format, s3$format))
        } else {
            if (!silent) warning("limited functionality for suffix that implies more than 3 instruments")
            type='synthetic' #condors, baskets, strips, packs/bundles,     
        }
    ##End check for suffixes with a dot
    } else if (any(substr(x,1,2) == c("1C","1D"))) { # Single-stock future (SSF) 
        #1CU1, 1CU11, 1CSEP11 -- 1DU1 for dividend protected
        if (substr(x,1,2) == "1C") {
            type <- c('outright', 'SSF')
        } else if (substr(x,1,2) == "1D") {
            type <- c('outright', 'SSF', 'NoDivRisk')
        }
        suff <- parse_suffix(substr(x,3,nchar(x)),silent=silent)
        month <- suff$month
        year <- suff$year
        format <- paste('1x', suff$format, sep="")
    } else if (nchar(x) == 2) { #U1
        if (substr(x,1,1) %in% M2C() && !is.na(suppressWarnings(as.numeric(substr(x,2,2))))) {
            type <- c("outright","future")
            month <- toupper(C2M(substr(x,1,1)))
            year <- as.numeric(substr(x,2,2)) + 2010
            if (!silent)
                warning("Converting 1 digit year to 4 digit year assumes there are no futures before 2010")
            format <- 'CY'
        } else if (suppressWarnings(is.na(as.numeric(x)))) type <- 'root'
    } else if (nchar(x) == 3) { #U11
        if (substr(x,1,1) %in% M2C() && !is.na(suppressWarnings(as.numeric(substr(x,2,3))))) {
            type <- c("outright","future")
            month <- toupper(C2M(substr(x,1,1)))
            year <- as.numeric(substr(x,2,3)) + 2000
            if (year > 2040) year <- year - 100
            if (!silent)
                warning('Converting 2 digit year to 4 digit year will result in a year between 1941 and 2040')
            format <- 'CYY'
        } else if (suppressWarnings(is.na(as.numeric(x)))) {
            type <- 'root'
        } else if (!silent) warning("Could not parse 3 character suffix")
    } else if (nchar(x) == 4) { #SEP1, VXU1, 0911
        if (toupper(substr(x, 1, 3)) %in% toupper(C2M()) && !is.na(suppressWarnings(as.numeric(substr(x,4,4))))) { 
            #sep1, Sep1, SEP1
            suff <- paste(M2C(tolower(substr(x,1,3))), substr(x,4,4), sep="") #convert from Sep1 to U1
            out <- parse_suffix(suff,silent=silent) #call recursively with 2 character suffix
            out$format <- 'MMMY'
            return(out) 
#        } else if (substr(x,3,3) %in% M2C() && !is.na(as.numeric(substr(x,4,4)))) {
#            #xxU1, VXU1 #ignore the 1st 2 characters, and call recursively with 2 character suffix            
#            suff <- parse_suffix(substr(x,3,4),silent=silent) 
#            month <- suff$month
#            year <- suff$year
#            format <- 'xxCY'
        } else if (!is.na(as.numeric(x))) { 
            #0911
            #convert to U11 and call recursively
            suff <- paste(M2C()[as.numeric(substr(x,1,2))], substr(x, 3,4), sep="")
            out <- parse_suffix(suff,silent=silent)
            out$format <- "NNNN"
            return(out) 
        } else {
            if (!silent)
                warning("Could not parse 4 character suffix")
            #return(NULL)
        }

    } else if (nchar(x) == 5) { #SEP11, U2011
        type <- c("outright","future")
        if (toupper(substr(x,1,3)) %in% toupper(C2M()) && !is.na(as.numeric(substr(x,4,5)))) {
            #SEP11        
            month <- toupper(substr(x,1,3))
            year <- as.numeric(substr(x,4,5)) + 2000
            if (!silent) warning('Converting 2 digit year to 4 digit year assumes there are no futures before 2000')
            format <- 'MMMYY'
       } else if (!is.na(as.numeric(substr(x,2,5))) && (substr(x,1,1) %in% M2C()) ) {
            #U2011
            month <- toupper(C2M(substr(x,1,1)))
            year <- as.numeric(substr(x,2,5))
            format <- 'CYYYY'
       }
    } else if (nchar(x) == 6) {
        #201109, 092011, 091611
        if (!silent)
            warning("I'm not programmed to handle a 6 character suffix")
        return(NULL)
            
    } else if (nchar(x) == 7) {
        #Sep2011
        if (toupper(substr(x, 1, 3)) %in% toupper(C2M()) && !is.na(as.numeric(substr(x,4,7))) ) {
            type <- c("outright","future")
            month <- toupper(substr(x,1,3))
            year <- as.numeric(substr(x, 4,7))
            format <- 'MMMYYYY'
        }    
    } 
    structure(list(type=type, month=month,year=year, strike=strike, right=right, cm=cm, cc=cc, format=format), class='suffix.list')
}


#' id.list class print method
#' 
#' @keywords internal
#' @export
print.id.list <- function(x, ...) {
  str(x, comp.str="", give.length=FALSE, give.attr=FALSE)
  invisible(x)
}

#' suffix.list class print method
#' 
#' @keywords internal
#' @export
print.suffix.list <- function(x, ...) {
  str(x, comp.str="", give.length=FALSE, give.attr=FALSE)
  invisible(x)
}


