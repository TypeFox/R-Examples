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
# $Id: load.instruments.R 1656 2014-11-28 01:07:57Z gsee $
#
###############################################################################


#' load instrument metadata into the .instrument environment
#' 
#' This function will load instrument metadata (data about the data)
#' either from a file specified by the \code{file} argument or
#' from a \code{data.frame} specified by the \code{metadata} argument.
#' 
#' The function will attempt to make reasonable assumptions about what you're trying to do, but this isn't magic.
#' 
#' You will typically need to specify the \code{type} of instrument to be loaded, failure to do so will generate a Warning and \code{default_type} will be used.
#' 
#' You will need to specify a \code{primary_id}, or define a \code{id_col} that contains the data to be used as the primary_id of the instrument.
#' 
#' You will need to specify a \code{currency}, unless the instrument \code{type} is 'currency'
#' 
#' Use the \code{identifier_cols} argument to specify which fields (if any) in the CSV are to be passed to \code{\link{instrument}} as the \code{identifiers} argument
#'
#' Typically, columns will exist for \code{multiplier} and \code{tick_size}.
#' 
#' Any other columns necessary to define the specified instrument type will also be required to avoid fatal Errors.  
#' 
#' Additional columns will be processed, either as additional identifiers for recognized identifier names, or as custom fields.  See \code{\link{instrument}} for more information on custom fields.
#' 
#' @param file string identifying file to load, default NULL, see Details
#' @param ... any other passthru parameters
#' @param metadata optional, data.frame containing metadata, default NULL, see Details
#' @param id_col numeric column containing id if primary_id isn't defined, default 1
#' @param default_type character string to use as instrument type fallback, see Details
#' @param identifier_cols character vector of field names to be passed as identifiers, see Details
#' @param overwrite TRUE/FALSE. See \code{\link{instrument}}.
#' @seealso 
#' \code{\link{loadInstruments}},
#' \code{\link{instrument}}, 
#' \code{\link{setSymbolLookup.FI}}, 
#' \code{\link[quantmod]{getSymbols}}, 
#' \code{\link{getSymbols.FI}}
#' @examples
#' \dontrun{
#' load.instruments(system.file('data/currencies.csv',package='FinancialInstrument'))
#' load.instruments(system.file('data/root_contracts.csv',package='FinancialInstrument'))
#' load.instruments(system.file('data/future_series.csv',package='FinancialInstrument'))
#'
#' }
#' @export
load.instruments <- function (file=NULL, ..., metadata=NULL, id_col=1, default_type='stock', identifier_cols=NULL, overwrite=TRUE) {

    if(is.null(file) && is.null(metadata)) stop("You must pass either a file identifier string or a metadata object to be converted.")
    if(is.null(metadata)){
        if (file.exists(file)){
            filedata<-read.csv(file,stringsAsFactors=FALSE, ...=...)
        } else {
            stop("The specified file ",file," does not seem to exist, maybe specify the full path?")
        }        
    } else {
        filedata<-metadata
        rm(metadata)
    }
    
    # check required column headers
    if(!any(grepl('primary_id',colnames(filedata)))) {
        #no primary_id column name, use id_col as the name
        set_primary<-TRUE
    } else {
        set_primary<-FALSE
    }
    dotargs<-list(...)
    if(!any(grepl('type',colnames(filedata)))) {
        if (is.null(dotargs$type)) {
            warning("metadata does not appear to contain instrument type, using ",default_type,". This may produce incorrect valuations.")
            filedata$type<-rep(default_type,nrow(filedata))
        } else {
            filedata$type <- rep(dotargs$type, nrow(filedata))
            dotargs$type <- NULL        
        }
    }
    if (!is.null(dotargs$currency) && !is.currency.name(dotargs$currency)) currency(dotargs$currency)
    
    #now process the data
    for(rn in 1:nrow(filedata)){
        type <- as.character(filedata[rn,'type'])
        arg <- as.list(filedata[rn,])
        if(type=='spread' || type=='guaranteed_spread'){
            if(!is.null(arg$members)){
                arg$members<-unlist(strsplit(arg$members,','))
            }
            if(!is.null(arg$memberratio)){
                arg$memberratio<-unlist(strsplit(arg$memberratio,','))
            }
            if(!is.null(arg$ratio)){
                arg$memberratio<-unlist(strsplit(arg$ratio,','))
            }
        }
        arg$type <- NULL
        arg <- arg[!is.na(arg)]
        arg <- arg[!arg==""]
        if (set_primary) {
            arg$primary_id<-filedata[rn,id_col]
        }
        
        #do some name cleanup to make up for Reuters silliness
        if(substr(arg$primary_id,1,1)==1) arg$primary_id <- substr(arg$primary_id,2,nchar(arg$primary_id))
        arg$primary_id<-make.names(arg$primary_id)
        if(!is.null(arg$X.RIC)){
            if(substr(arg$X.RIC,1,1)==1) arg$X.RIC <- substr(arg$X.RIC,2,nchar(arg$X.RIC))
        }            
        if(!is.null(arg$RIC)){
            if(substr(arg$RIC,1,1)==1) arg$RIC <- substr(arg$RIC,2,nchar(arg$RIC))
        }            
        if(length(dotargs)) arg<-c(arg,dotargs)
        
        if(!is.null(identifier_cols) && any(identifier_cols %in% names(arg))){
            arg$identifiers <- arg[names(arg) %in% identifier_cols]
            arg[identifier_cols] <- NULL
        }
        
        arg$overwrite <- overwrite
        if(is.function(try(match.fun(type),silent=TRUE))){
            out <- try(do.call(type,arg))
            
            
            #TODO recover gracefully?
        } else {
            # the call for a function named for type didn't work, so we'll try calling instrument as a generic
            type=c(type,"instrument")
            arg$type<-type # set the type
            arg$assign_i<-TRUE # assign to the environment
            try(do.call("instrument",arg))
        }
    } # end loop on rows
}

#' set quantmod-style SymbolLookup for instruments
#' 
#' This function exists to tell \code{\link[quantmod]{getSymbols}} where to look for your repository of market data.
#' 
#' The \code{base_dir} parameter \emph{must} be set or the function will fail.  
#' This will vary by your local environment and operating system.  For mixed-OS environments,
#' we recommend doing some OS-detection and setting the network share to your data to a common 
#' location by operating system.  For example, all Windows machines may use \dQuote{M:/} 
#' and all *nix-style (linux, Mac) machines may use \dQuote{/mnt/mktdata/}. 
#' 
#' The \code{split_method} currently allows either \sQuote{days} or \sQuote{common}, and expects the 
#' file or files to be in sub-directories named for the symbol.  In high frequency data, it is standard practice to split
#' the data by days, which is why that option is the default.
#'     
#' @param base_dir string specifying the base directory where data is stored, see Details 
#' @param Symbols character vector of names of instruments for which to \code{setSymbolLookup}
#' @param \dots any other passthru parameters
#' @param storage_method currently only \sQuote{rda}, but we will eventually support \sQuote{indexing} at least, and maybe others
#' @param split_method string specifying the method files are split, currently \sQuote{days} or \sQuote{common}, see Details
#' @param use_identifier string identifying which column should be use to construct the \code{primary_id} of the instrument, default 'primary_id'
#' @param extension file extension, default "rda"
#' @param src which \code{\link[quantmod]{getSymbols}} sub-type to use, default \code{\link{getSymbols.FI}} by setting 'FI'
#' @seealso 
#' \code{\link{getSymbols.FI}},
#' \code{\link{instrument_attr}},
#' \code{\link{load.instruments}}, \code{\link{loadInstruments}},
#' \code{\link[quantmod]{setSymbolLookup}}
#' @importFrom zoo as.Date
#' @export
setSymbolLookup.FI<-function(base_dir, Symbols, ..., split_method=c("days","common"), storage_method='rda', use_identifier='primary_id', extension='rda', src='FI'){
    # check that base_dir exists
    if(!file.exists(base_dir)) stop('base_dir ',base_dir,' does not seem to specify a valid path' )
    
    # take split
    split_method<-split_method[1] # only use the first value

    #load all instrument names
    instr_names <- if(missing(Symbols)) {
        ls_non_currencies(ls(pos=.instrument)) #if roots begin with a dot, this will filter out roots and currencies
    } else Symbols
    
    #TODO add check to make sure that src is actually the name of a getSymbols function
    
    #initialize list
    params<-list()
    params$storage_method<-storage_method
    params$extension<-extension
    params$split_method<-split_method
    params$src<-src
    if(length(list(...))>=1){
        dlist<-list(...)
        params<-c(params,dlist)
    }
    new.symbols<-list()
    ndc<-nchar(base_dir)
    if(substr(base_dir,ndc,ndc)=='/') sepch='' else sepch='/'
    for (instr in instr_names){
        tmp_instr<-getInstrument(instr)
        if(!use_identifier=='primary_id'){
            instr_str<-make.names(tmp_instr$identifiers[[use_identifier]])
        } else {
            instr_str<-make.names(tmp_instr[[use_identifier]])
        } 
        if(!is.null(instr_str)) instr<-instr_str
        symbol<-list()
        symbol[[1]]<-params
        # construct $dir
        symbol[[1]]$dir<-base_dir
        names(symbol)[1]<-instr
        new.symbols<-c(new.symbols,symbol)
    }
    setSymbolLookup(new.symbols)
}


#' getSymbols method for loading data from split files
#' 
#' This function should probably get folded back into getSymbols.rda in 
#' quantmod.
#' 
#' Meant to be called internally by \code{\link[quantmod]{getSymbols}} .
#' 
#' The symbol lookup table will most likely be loaded by 
#' \code{\link{setSymbolLookup.FI}}
#' 
#' If date_format is NULL (the Default), we will assume an ISO date as changed 
#' by \code{\link{make.names}}, for example, 2010-12-01 would be assumed to be a 
#' file containing 2010.12.01
#' 
#' If \code{indexTZ} is provided, the data will be converted to that timezone
#'
#' If auto.assign is FALSE, \code{Symbols} should be of length 1.  Otherwise, 
#' \code{\link[quantmod]{getSymbols}} will give you an error that says 
#' \dQuote{must use auto.assign=TRUE for multiple Symbols requests}
#' However, if you were to call \code{getSymbols.FI} directly (which is 
#' \emph{NOT} recommended) with \code{auto.assign=FALSE} and more than one 
#' Symbol, a list would be returned.
#' 
#' Argument matching for this function is as follows.  If the user provides a 
#' value for an argument, that value will be used.  If the user did not provide
#' a value for an argument, but there is a value for that argument for the 
#' given \code{Symbol} in the Symbol Lookup Table (see 
#' \code{\link{setSymbolLookup.FI}}), that value will be used.  Otherwise, 
#' the formal defaults will be used.
#'
#' @param Symbols a character vector specifying the names of each symbol to be 
#'   loaded
#' @param from Retrieve data no earlier than this date. Default '2010-01-01'.
#' @param to Retrieve data through this date. Default Sys.Date().
#' @param ... any other passthru parameters
#' @param dir if not specified in getSymbolLookup, directory string to use.  
#'   default ""
#' @param return.class only "xts" is currently supported
#' @param extension file extension, default "rda"
#' @param split_method string specifying the method used to split the files, 
#'   currently \sQuote{days} or \sQuote{common}, see 
#'   \code{\link{setSymbolLookup.FI}}
#' @param use_identifier optional. identifier used to construct the 
#'   \code{primary_id} of the instrument. If you use this, you must have 
#'   previously defined the \code{\link{instrument}} 
#' @param date_format format as per the \code{\link{strptime}}, see Details
#' @param verbose TRUE/FALSE
#' @param days_to_omit character vector of names of weekdays that should not be 
#'   loaded.  Default is \code{c("Saturday", "Sunday")}.  Use \code{NULL} to 
#'   attempt to load data for all days of the week.
#' @param indexTZ valid TZ string. (e.g. \dQuote{America/Chicago} or 
#'   \dQuote{America/New_York}) See \code{\link[xts]{indexTZ}}.
#' @seealso 
#' \code{\link{saveSymbols.days}}
#' \code{\link{instrument}}
#' \code{\link{setSymbolLookup.FI}}
#' \code{\link{loadInstruments}}
#' \code{\link[quantmod]{getSymbols}}
#' @examples
#' \dontrun{
#' getSymbols("SPY", src='yahoo')
#' dir.create("tmpdata")
#' saveSymbols.common("SPY", base_dir="tmpdata")
#' rm("SPY")
#' getSymbols("SPY", src='FI', dir="tmpdata", split_method='common')
#' unlink("tmpdata/SPY", recursive=TRUE)
#' }
#' @export
getSymbols.FI <- function(Symbols,
                            from=getOption("getSymbols.FI.from", "2010-01-01"),
                            to=getOption("getSymbols.FI.to", Sys.Date()),
                            ..., 
                            dir=getOption("getSymbols.FI.dir", ""),
                            return.class=getOption("getSymbols.FI.return.class", 
                                                   "xts"),
                            extension=getOption("getSymbols.FI.extension", "rda"),
                            split_method=getOption("getSymbols.FI.split_method",
                                                   c("days", "common")),
                            use_identifier=getOption("getSymbols.FI.use_identifier",
                                                     NA),
                            date_format=getOption("getSymbols.FI.date_format"),
                            verbose=getOption("getSymbols.FI.verbose", TRUE),
                            days_to_omit=getOption("getSymbols.FI.days_to_omit",
                                                   c("Saturday", "Sunday")),
                            indexTZ=getOption("getSymbols.FI.indexTZ", NA)
                         )
{
    if (is.null(date_format)) date_format <- "%Y.%m.%d"
    if (is.null(days_to_omit)) days_to_omit <- 'NULL'
    this.env <- environment()
    for(var in names(list(...))) {
        assign(var,list(...)[[var]], this.env)
    }

    #The body of the following function comes from Dominik's answer here: 
    #browseURL("http://stackoverflow.com/questions/7224938/can-i-rbind-be-parallelized-in-r")
    #it does what do.call(rbind, lst) would do, but faster and with less memory usage
    do.call.rbind <- function(lst) {
        while(length(lst) > 1) {
            idxlst <- seq(from=1, to=length(lst), by=2)

            lst <- lapply(idxlst, function(i) {
                if(i==length(lst)) { return(lst[[i]]) }

                return(rbind(lst[[i]], lst[[i+1]]))
            })
        }
        lst[[1]]
    }

    # Find out if user provided a value for each formal
    if (hasArg.from <- hasArg(from)) .from <- from
    if (hasArg.to <- hasArg(to)) .to <- to
    if (hasArg.dir <- hasArg(dir)) .dir <- dir
    if (hasArg.return.class <- hasArg(return.class)) 
        .return.class <- return.class
    if (hasArg.extension <- hasArg(extension)) .extension <- extension
    if (hasArg.split_method <- hasArg(split_method)) 
        .split_method <- split_method
    if (hasArg.use_identifier <- hasArg(use_identifier)) 
        .use_identifier <- use_identifier
    if (hasArg.date_format <- hasArg(date_format)) .date_format <- date_format
    if (hasArg.verbose <- hasArg(verbose)) .verbose <- verbose
    if (hasArg.days_to_omit <- hasArg(days_to_omit)) 
        .days_to_omit <- days_to_omit
    if (hasArg.indexTZ <- hasArg(indexTZ)) .indexTZ <- indexTZ

    importDefaults("getSymbols.FI")

    # Now get the values for each formal that we'll use if not provided
    # by the user and not found in the SymbolLookup table
    default.from <- from
    default.to <- to
    default.dir <- dir
    default.return.class <- return.class
    default.extension <- extension
    default.split_method <- split_method[1]
    default.use_identifier <- use_identifier
    default.date_format <- date_format
    default.verbose <- verbose
    default.days_to_omit <- days_to_omit
    default.indexTZ <- indexTZ
    
    # quantmod:::getSymbols will provide auto.assign and env
    # so the next 2 if statements should always be TRUE
    auto.assign <- if(hasArg(auto.assign)) {auto.assign} else TRUE
    env <- if(hasArg(env)) {env} else .GlobalEnv 

    # make an argument matching function to sort out which values to use for each arg
    pickArg <- function(x, Symbol) {
        if(get(paste('hasArg', x, sep="."))) {
            get(paste(".", x, sep=""))
        } else if(!is.null(SymbolLookup[[Symbol]][[x]])) {
            SymbolLookup[[Symbol]][[x]]
        } else get(paste("default", x, sep="."))
    }

    SymbolLookup <- getSymbolLookup()
    fr <- NULL
    datl <- lapply(1:length(Symbols), function(i) {
        #FIXME? Should nothing be saved if there are errors with any of 
        # the Symbols (current behavior)?  Or, if auto.assign == TRUE, should
        # we assign the data as we get it instead of making a list of data and 
        # assigning at the end.
        from <- pickArg("from", Symbols[[i]])
        to <- pickArg("to", Symbols[[i]])
        dir <- pickArg("dir", Symbols[[i]])
        return.class <- pickArg("return.class", Symbols[[i]])
        extension <- pickArg('extension', Symbols[[i]])
        split_method <- pickArg('split_method', Symbols[[i]])
        use_identifier <- pickArg('use_identifier', Symbols[[i]])
        date_format <- pickArg('date_format', Symbols[[i]])
        verbose <- pickArg('verbose', Symbols[[i]])
        days_to_omit <- pickArg('days_to_omit', Symbols[[i]])
        indexTZ <- pickArg('indexTZ', Symbols[[i]])
        # if 'dir' is actually the 'base_dir' then we'll paste the instrument name (Symbol) to the end of it.
        # First, find out what the instrument name is
        instr_str <- NA
        if(!is.na(use_identifier)) { 
            tmp_instr <- try(getInstrument(Symbols[[i]], silent=FALSE))
            if (inherits(tmp_instr,'try-error') || !is.instrument(tmp_instr)) 
                stop("must define instrument first to call with 'use_identifier'")
            if (!use_identifier=='primary_id') {
                instr_str<-make.names(tmp_instr$identifiers[[use_identifier]])
            } else  instr_str <- make.names(tmp_instr[[use_identifier]])
            if (length(instr_str) == 0L) stop("Could not find instrument. Try with use_identifier=NA")
        }
        Symbol <- ifelse(is.na(instr_str), make.names(Symbols[[i]]), instr_str)
        ndc<-nchar(dir)
        if(substr(dir,ndc,ndc)=='/') dir <- substr(dir,1,ndc-1) #remove trailing forward slash
        dir <- paste(dir, Symbol, sep="/")
        
        if(!dir=="" && !file.exists(dir)) {
            if (verbose) cat("\ndirectory ",dir," does not exist, skipping\n")
        } else {
            if(verbose) cat("loading ",Symbols[[i]],".....\n")
            switch(split_method[1],
                    days={
                        StartDate <- as.Date(from) 
                        EndDate <- as.Date(to) 
                        date.vec <- as.Date(StartDate:EndDate)
                        date.vec <- date.vec[!weekdays(date.vec) %in% days_to_omit]  
                        date.vec <- format(date.vec, format=date_format)
                        sym.files <- paste(date.vec, Symbol, extension, sep=".")
                        if (dir != "") sym.files <- file.path(dir, sym.files)
                        dl <- lapply(sym.files, function(fp) {
                            sf <- strsplit(fp, "/")[[1]]
                            sf <- sf[length(sf)]
                            if (verbose) cat("Reading ", sf, "...")                            
                            if(!file.exists(fp)) {       
                                if (verbose) cat(" failed. File not found in ", dir, " ... skipping\n")
                            } else {
                                if (verbose) cat(' done.\n')
                                local.name <- load(fp)
                                dat <- get(local.name)
                                if (!is.na(indexTZ) && !is.null(dat)) indexTZ(dat) <- indexTZ
                                dat
                            }
                        })
                        if (verbose) cat('rbinding data ... ')
                        fr <- do.call.rbind(dl)
                    },
                    common = , {
                        sym.file <- paste(Symbol,extension,sep=".")
                        if(dir != "") sym.file <- file.path(dir, sym.file)
                        if(!file.exists(sym.file)) {
                            if (verbose) cat("file ",paste(Symbol,extension,sep='.')," does not exist in ",dir,"....skipping\n")
                        } else {
                            #fr <- read.csv(sym.file)
                            local.name <- load(sym.file)
                            dat <- get(local.name)
                            if (!is.na(indexTZ) && !is.null(dat)) indexTZ(dat) <- indexTZ
                            assign('fr', dat)
                            if(verbose) cat("done.\n")
                            #if(!is.xts(fr)) fr <- xts(fr[,-1],as.Date(fr[,1],origin='1970-01-01'),src='rda',updated=Sys.time())
                        }
                    } # end 'common'/default method (same as getSymbols.rda)    
                ) # end split_method switch
            fr <- convert.time.series(fr=fr,return.class=return.class)
            Symbols[[i]] <-make.names(Symbols[[i]]) 
            tmp <- list()
            tmp[[Symbols[[i]]]] <- fr
            if(verbose) cat("done.\n")
            tmp     
        }
    }) #end loop over Symbols

    if (length(Filter("+", lapply(datl, length))) == 0) {
        warning("No data found.")
        return(NULL) 
    }

    datl.names <- do.call(c, lapply(datl, names))
    missing <- Symbols[!Symbols %in% datl.names]
    if (length(missing) > 0) warning('No data found for ', paste(missing, collapse=" "))
    if(auto.assign) {
        #invisible(lapply(datl, function(x) if (length(x) > 0) assign(names(x), x[[1]], pos=env)))
        out <- Filter(function(x) length(x) > 0, datl)
        invisible(lapply(out, function(x) assign(names(x), x[[1]], pos=env)))
        return(datl.names)
    } else {
        #NOTE: Currently, NULLs aren't filtered out.  If there are data for any Symbol,
        # the returned list will have an element for each symbol requested even if some don't contain data.
        out <- lapply(datl, function(x) {
            if (length(x) > 0) x[[1]]
        })
        if (length(out) == 1)
            return(out[[1]])
        else {
            names(out) <- Symbols
            return(out)
        }
    }
}


#' currency metadata to be used by \code{\link{load.instruments}}
#'
#' @name currencies
#' @docType data
#' @keywords data
NULL

#' future metadata to be used by \code{\link{load.instruments}}
#'
#' @name root_contracts
#' @docType data
#' @keywords data
NULL

