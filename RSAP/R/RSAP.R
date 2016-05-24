# file RSAP/R/RSAP.R
# copyright (C) 2012 Piers Harding
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#
#  Low level wrappers for odbc driver
#
#
#
.onLoad <- function(libname, pkgname)
{
    if(is.null(getOption("dec")))
        options(dec = Sys.localeconv()["decimal_point"])
}

.onUnload <- function(libpath)
{
    .Call(C_RSAPTerm)
    library.dynam.unload("RSAP", libpath)
}


RSAPConnect <- function (...)
{
    args <- list(...)
    if (length(args) == 0) {
        stop("No arguments supplied")
    }
    if (typeof(args[[1]]) == "list") {
        args = args[[1]]
    }
    # did we get passed a config file?
    if (typeof(args[[1]]) == "character" && file.exists(args[[1]])) {
        library(yaml)
        config <- yaml.load_file(args[[1]])
        newargs <- list()
        for (x in names(config)) { newargs[[x]] <- as.character(config[[x]]); }
        return(RSAPConnect(newargs))
    }
    # format client
    if (exists("client", where=args)) {
        args[['client']] <- sprintf("%03d", as.integer(args[['client']]))
    }
    # format sysnr
    if (exists("sysnr", where=args)) {
        args[['sysnr']] <- sprintf("%02d", as.integer(args[['sysnr']]))
    }
    res <- .Call(C_RSAPRFCConnect, args)
    return(res)
}

RSAPshowArgs <-
    function(...)
{
    res <- .Call(C_RSAPshowArgs, list(...))
    #warning(res)
    return(res)
}

RSAPValidHandle <-  function(con)
{
    if (!is.integer(con)) {
        print("con is not an integer")
        return(FALSE)
    }
    if (!is.element("handle_ptr", names(attributes(con)))) {
        print("handle_ptr does not exist")
        return(FALSE)
    }
    res <- .Call(C_RSAPValidHandle, con)
    return(res)
}


RSAPGetInfo <- function(con)
{
    if(!RSAPValidHandle(con))
       stop("argument is not a valid RSAP con")
    res <- .Call(C_RSAPGetInfo, con)
    return(res)
}


RSAPInvoke <- function(con, func, parms)
{
    if(!RSAPValidHandle(con))
       stop("argument is not a valid RSAP con")
    res <- .Call(C_RSAPInvoke, con, func, parms)
    return(res)
}


RSAPReadTable <- function(con, saptable, options=list(), fields=list())
{
    if(!RSAPValidHandle(con))
       stop("argument is not a valid RSAP con")
    library(reshape)
    parms <- list('DELIMITER' = ';',
              'QUERY_TABLE' = saptable,
              'OPTIONS' = list('TEXT' = options),
              'FIELDS' = list('FIELDNAME' = fields)
              )
	res <- RSAPInvoke(con, "RFC_READ_TABLE", parms)
	flds <- sub("\\s+$", "", res$FIELDS$FIELDNAME)
    data <- NULL
    if (length(res$DATA$WA) == 0) {
        data <- data.frame()
    }
    else {
	    data <- data.frame(res$DATA, colsplit(res$DATA$WA, split = ";", names = flds))
    }
	
	for (i in 1:length(flds)) {
		f <- flds[i]
		typ <- res$FIELDS$TYPE[i]
		if (typ == 'N' || typ == 'I' || typ == 'P') {
			data[[f]] <- as.numeric(unlist(lapply(data[[f]], FUN=function (x) {sub("[^\\d\\.\\-\\,]", "", x, perl=TRUE)})));
			data[[f]][is.na(data[[f]])] <- 0 
		} else {
		    data[[f]] <- sub("\\s+$", "", data[[f]]);
		}
	}
	data$WA <- NULL
    return(data)
}


RSAPReadCube <- function(con, cube, ref_date=NULL, chars=list(), kfigures=list(), options=list())
{
	if(!RSAPValidHandle(con))
		stop("argument is not a valid RSAP con")
	
	if (is.null(ref_date))
		ref_date <- format(Sys.time(), "%Y%m%d")
	
	# RSSEM_CHA_VALUES_GET - I_CHANM, I_DATETO=ref_date, I_READ_TEXT='X'
	
	caliases <- as.list(sub("^\\d+", "", chars))
	kaliases <- as.list(sub("^\\d+", "", kfigures))
	parms <- list(
			'I_INFOCUBE' = cube,
			'I_REFERENCE_DATE' = ref_date,
			'I_T_CHA' = list('CHANM' = chars, 'CHAALIAS' = caliases),
			'I_T_KYF' = list('KYFNM' = kfigures, 'KYFALIAS' = kaliases),
			'I_T_RANGE' = options
	)
	res <- RSAPInvoke(con, "RSDPL_CUBE_DATA_READ", parms)
	cube <- RSAPReadTable(con, res$E_TABLENAME)
	# add the text label names as another column
	for (el in chars[chars != '0CALMONTH']) {
		v <- sub("^\\d+", "", el) # the data.frame vector name
		vt <- paste(v, "_TEXT", sep="")
		cube[[vt]] <- cube[[v]]
		ref_date <- format(Sys.time(), "%Y%m%d")
		parms <- list('I_CHANM' = el,
				'I_READ_TEXT' = 'X',
				'I_DATETO' = ref_date)
		res <- RSAPInvoke(con, "RSSEM_CHA_VALUES_GET", parms)
		res$ETH_CHAVL$IOBJNM   <- sub("\\s+$", "", res$ETH_CHAVL$IOBJNM);
		res$ETH_CHAVL$VALUE   <- sub("\\s+$", "", res$ETH_CHAVL$VALUE);
		elements <- res$ETH_CHAVL$VALUE[res$ETH_CHAVL$IOBJNM == el]
		nms <- res$ETH_CHAVL$VALUE[res$ETH_CHAVL$IOBJNM == '0TXTSH']
		for (k in elements) {
			if (length(cube[[v]][cube[[v]] == k]) > 0) {
				cube[[vt]][cube[[v]] == k] <- nms[elements == k]
			}
		}
	}
	
	return(cube)
}


RSAPGetCube <- function(con, cube)
{
    if(!RSAPValidHandle(con))
       stop("argument is not a valid RSAP con")

    parms <- list(
              'INFOCUBE' = cube,
              'OBJVERS' = as.list('A')
              )
	res <- RSAPInvoke(con, "BAPI_CUBE_GETDETAIL", parms)
    res$INFOOBJECTS$INFOCUBE   <- sub("\\s+$", "", res$INFOOBJECTS$INFOCUBE);
    res$INFOOBJECTS$INFOOBJECT <- sub("\\s+$", "", res$INFOOBJECTS$INFOOBJECT);
    res$RETURN <- NULL
    res$OBJVERS <- NULL
    res$INFOCUBE <- NULL
    return(res)
}


RSAPListCubes <- function(con)
{
    if(!RSAPValidHandle(con))
       stop("argument is not a valid RSAP con")

    parms <- list(
              'CUBETYPE' = as.list('%'),
              'OBJVERS' = as.list('A')
              )
	res <- RSAPInvoke(con, "BAPI_CUBE_GETLIST", parms)
    res$INFOCUBELIST$INFOCUBE <- sub("\\s+$", "", res$INFOCUBELIST$INFOCUBE);
    res$INFOCUBELIST$TEXTLONG <- sub("\\s+$", "", res$INFOCUBELIST$TEXTLONG);
    res$INFOCUBELIST$INFOAREA <- sub("\\s+$", "", res$INFOCUBELIST$INFOAREA);
    return(as.data.frame(res$INFOCUBELIST))
}

     
RSAPExecInfoQuery <- function(con, infoprovider, infoquery)
{
    if(!RSAPValidHandle(con))
       stop("argument is not a valid RSAP con")

    parms <- list(
              'I_INFOPROVIDER' = infoprovider,
              'I_QUERY' = infoquery
              )
	res <- RSAPInvoke(con, "RRW3_GET_QUERY_VIEW_DATA", parms)
    return(res)
}

print.RSAP_Connector <- function(x, ...) RSAPGetInfo(x)

close.RSAP_Connector <- function(con, ...) RSAPClose(con)

RSAPClose <- function(con)
{
    if(!RSAPValidHandle(con))
       stop("argument is not a valid RSAP con")
    res <- .Call(C_RSAPClose, con)
    return(res)
}

readTable <- function(con, ...){
                        res <- RSAPReadTable(con, ...)
                        return(res)
                }

readCube <- function(con, ...){
                        res <- RSAPReadCube(con, ...)
                        return(res)
                }

