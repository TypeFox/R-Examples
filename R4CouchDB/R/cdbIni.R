#' Ini function
#'
#' Function returns a list with some default settings and often used functions
#' such as \code{cdb$baseUrl}.
#'
#' The list: \code{ cdb <- list(serverName = "localhost", ... )} is returned if
#' the packages \code{library(RCurl)} and \code{library(RJSONIO)} are
#' successfully loaded.
#'
#' @author wactbprot, parisni 
#' @export
#' @usage cdbIni(serverName="localhost",
#' port="5984",
#' prot = "http",
#' DBName="",
#' uname = "",
#' pwd = "",
#' newDBName = "",
#' removeDBName = "",
#' id  = "",
#' fileName = "",
#' design = "",
#' view = "",
#' list = "",
#' queryParam = "",
#' encSub = "?",
#' count = 10,
#' dataList = list(),
#' attachmentsWithPath=TRUE,
#' digits = 7)
#' @param serverName server name
#' @param port port
#' @param prot name of the protocol default is http
#' @param DBName name of database
#' @param uname name of the user
#' @param pwd password
#' @param newDBName name of the database for cdbMakeDB()
#' @param removeDBName name of the database to remove with cdbRemoveDB()
#' @param id the document id to get, put, post or delete
#' @param fileName for use in cdbAddAttachment
#' @param design the name of the design used when asking a view or list
#' @param view the name of the view to query
#' @param list the name of the list to query
#' @param queryParam additional query params
#' @param encSub a character which is used as a replacement for chars who can not be converted by iconv
#' @param count how many uuids should be returned by cdbGetUuidS()
#' @param dataList a list containing data to post or update
#' @param attachmentsWithPath effects the result of the function cdbAddAttachment in the way the variable is named
#' @param digits digits kept at toJSON conversion
#' @import RCurl RJSONIO bitops
#' @export
#' @examples
#'\dontrun{
#' ccc <- cdbIni(digits=13,
#'               DBName="r4couch_db",
#'               attachmentsWithPath=FALSE,
#'               dataList=list(normalDistRand =  rnorm(20)))
#'}
#' @return \item{cdb}{The R4CouchDB (method) chain(ing) list }
#' @keywords misc
#'

cdbIni <- function(serverName   = "localhost",
                   port         = "5984",
                   prot         = "http",
                   DBName       = "",
                   uname        = "",
                   pwd          = "",
                   newDBName    = "",
                   removeDBName = "",
                   id           = "",
                   fileName     = "",
                   design       = "",
                   view         = "",
                   list         = "",
                   queryParam   = "",
                   encSub       = "?",
                   count        = 10,
                   dataList     = list(),
                   attachmentsWithPath = TRUE,
                   digits       = 7){


    cdb <- list(
        DBName       = DBName,
        serverName   = serverName,
        prot         = prot,
        port         = port,
        uname        = uname,
        pwd          = pwd,
        newDBName    = newDBName,
        removeDBName = removeDBName,
        id           = id,
        dataList     = dataList,
        fileName     = fileName,
        design       = design,
        view         = view,
        list         = list,
        queryParam   = queryParam,
        count        = count,
        encSub       = encSub,
        error        = "",
        res          = "",
        date         = toString(Sys.Date()),
        curl         = getCurlHandle(),
        localEnc     = "UTF-8",
        serverEnc    = "UTF-8",
        attachmentsWithPath = TRUE,
        digits       = digits)

    cdb$opts <- function(cdb){
        if(cdb$uname != "" & cdb$pwd != ""){
            opts <- curlOptions(header   = FALSE,
                                httpauth = 1L,
                                userpwd  = paste(cdb$uname,
                                    ":",
                                    cdb$pwd,
                                    sep="")
                                )
        }else{
            opts <- curlOptions(header = FALSE)
        }
        
        return(opts)
    }

    cdb$baseUrl <- function(cdb){
        return(paste(cdb$prot,
                     "://",
                     cdb$serverName,
                     ":",
                     cdb$port,
                     "/",
                     sep="")
               )
    }

    cdb$fromJSON <- function(jsn){
        jsn <- iconv(jsn,
                     cdb$serverEnc,
                     cdb$localEnc,
                     sub=cdb$encSub)

        return(fromJSON(jsn,
                        nullValue         = NA,
                        simplify          = FALSE,
                        simplifyWithNames = FALSE))
    }

    cdb$toJSON <- function(lst){
        jsn <- toJSON(lst,
                      collapse = "",
                      digits   = digits)
        jsn <- iconv(jsn,
                     cdb$localEnc,
                     cdb$serverEnc,
                     sub=cdb$encSub)

        ## one can {"a":"\r"} have in the
        ## database but one can not send it
        ## this way. A \r is here replaced by \\r
        ## resulting in \r in the database
        jsn <- gsub("\\r","\\\\r",jsn)

        return(jsn)
    }

    cdb$getDocRev <- function(cdb){
        adrString <- paste(cdb$baseUrl(cdb),
                           cdb$DBName,"/",
                           cdb$id,
                           sep="")
        res <- url.exists(adrString, .header=TRUE)["ETag"]
        if(is.na(res)){
            return(NA)
        }else{
            return(paste("", gsub("\\\"", "", res), sep = ""))
        }
    }

    cdb$checkRes <- function(cdb,res){
        if(!(cdb$error == "")){
            stop( paste("local error:", cdb$error))
        }

        res <- cdb$fromJSON(res)

        if(length(res$error) > 0){
            stop(paste("server error:", res$error,
                       "server reason:", res$reason))
        }else{
            cdb$res <- res
            return( cdb )
        }
    }

    cdb$checkCdb <- function(cdb,fname){
	    switch(fname,
		   cdbGetDoc = {
			   cdb <- chk.server.name(cdb)
			   cdb <- chk.id(cdb)
			   cdb <- chk.db.name(cdb)
		   },
		   cdbAddAttachment = {
			   cdb <- chk.server.name(cdb)
			   cdb <- chk.id(cdb)
			   cdb <- chk.db.name(cdb)
			   cdb <- chk.doc.exists(cdb)
			   cdb <- chk.file.name(cdb)
		   },
		   cdbAddDoc = {
			   cdb <- chk.server.name(cdb)
			   cdb <- chk.db.name(cdb)
			   cdb <- chk.data.list(cdb)
		   },
		   cdbDeleteDoc = {
			   cdb <- chk.server.name(cdb)
			   cdb <- chk.db.name(cdb)
			   cdb <- chk.id(cdb)
		   },
		   cdbGetConfig = {
			   cdb <- chk.server.name(cdb)
		   },
		   cdbGetList = {
			   cdb <- chk.server.name(cdb)
			   cdb <- chk.db.name(cdb)
			   cdb <- chk.design.name(cdb)
			   cdb <- chk.list.name(cdb)
			   cdb <- chk.view.name(cdb)
		   },
		   cdbGetUuid = {
			   cdb <- chk.server.name(cdb)
		   },
		   cdbGetUuidS = {
			   cdb <- chk.server.name(cdb)
			   cdb <- chk.count(cdb)
		   },
		   cdbGetView = {
			   cdb <- chk.server.name(cdb)
			   cdb <- chk.db.name(cdb)
			   cdb <- chk.design.name(cdb)
			   cdb <- chk.view.name(cdb)
		   },
		   cdbListDB = {
			   cdb <- chk.server.name(cdb)
		   },
		   cdbMakeDB = {
			   cdb <- chk.server.name(cdb)
			   cdb <- chk.newdb.name(cdb)
		   },
		   cdbRemoveDB = {
			   cdb <- chk.server.name(cdb)
			   cdb <- chk.rmdb.name(cdb)
		   },
		   cdbUpdateDoc = {
			   cdb <- chk.server.name(cdb)
			   cdb <- chk.id(cdb)
			   cdb <- chk.db.name(cdb)
			   cdb <- chk.data.list(cdb)
		   },
		   cdbAddDocS = {
			   cdb <- chk.server.name(cdb)
			   cdb <- chk.db.name(cdb)
			   cdb <- chk.data.lists(cdb)
		   },
		   stop(paste0(fname," is not defined"))
		   )

        return(cdb)
    }

    ## ----------------------chk.fns-----------------v
    chk.count <- function(cdb){
        if(!is.numeric(cdb$count) | (cdb$count  < 1)){
            
            cdb$error <- paste(cdb$error,
                               ";cdb$count is not numeric or smaller than 1")
        }
        return(cdb)

    }
    chk.doc.exists <- function(cdb){
        res <- cdb$getDocRev(cdb)

        if(is.na(res)){
            cdb$error <- paste(cdb$error,
                               ";document cdb$id does not exist")
        }
        return(cdb)
     
    }
    chk.newdb.name <- function(cdb){
        if(cdb$newDBName == ""){
            cdb$error <- paste(cdb$error,
                               ";no cdb$newDBName given")
        }else{
            DBNames  <- cdbListDB(cdb)$res
            DBexists <- which(DBNames == cdb$newDBName)

            if(length(DBexists) > 0){
                cdb$error <- paste(cdb$error,
                                   ";cdb$newDBName already exists ")
            }
        }
        return(cdb)
    }

    chk.rmdb.name <- function(cdb){
        if(cdb$removeDBName == ""){
            cdb$error <- paste(cdb$error, ";no cdb$removeDBName given")
        }else{
            DBNames <- cdbListDB(cdb)$res
            DBexists <- which(DBNames == cdb$removeDBName)

            if(length(DBexists) == 0){
                cdb$error <- paste(cdb$error,
                                   ";cdb$removeDBName does not exist")
            }
        }
        return(cdb)
    }

    chk.design.name <- function(cdb){
        if(cdb$design == "") {
            cdb$error <- paste(cdb$error,
                               ";no cdb$design given")
        }
        return(cdb)
    }

    chk.list.name <- function(cdb){
        if(cdb$list == "") {
            cdb$error <- paste(cdb$error,
                               ";no cdb$list given")
        }
        return(cdb)
    }

    chk.view.name <- function(cdb){
        if(cdb$view == "") {
            cdb$error <- paste(cdb$error,
                               ";no cdb$view given")
        }
        return(cdb)
    }


    chk.data.list <- function(cdb){
        if( (length(cdb$dataList) < 1)){
            cdb$error <- paste(cdb$error,
                               ";no cdb$dataList given")
        }
        return(cdb)
    }

    chk.data.lists <- function(cdb){

        if(!is.list(cdb$dataList)){
            cdb$error <- paste(cdb$error,
                               ";cdb$dataList is not a list")
        }else{
            if((length(cdb$dataList) < 1)){
                cdb$error <- paste(cdb$error,
                                   ";cdb$dataList has length zero")
            }else{
                if(!is.list(cdb$dataList[[1]])){
                    cdb$error <- paste(cdb$error,
                                       ";cdb$dataList is not a list of lists")
                    
                }
            }
        }
        
        return(cdb)
    }

    chk.file.name <- function(cdb){
        if( !(file.exists(cdb$fileName))){
            cdb$error <- paste(cdb$error,
                               ";no cdb$fileName given or does not exist")
        }
        return(cdb)
    }
    chk.server.name <- function(cdb){
        if(cdb$serverName == ""){
            cdb$error <- paste(cdb$error,
                               ";no cdb$serverName given")
        }
        return(cdb)
    }

    chk.db.name <- function(cdb){
        if(cdb$DBName == ""){
            cdb$error <- paste(cdb$error,
                               ";no cdb$DBName given")
        }
        return(cdb)
    }

    chk.id <- function(cdb){
        if( cdb$id == ""){
            cdb$error <- paste(cdb$error,
                               ";no cdb$id given ")
        }
        return(cdb)
    }
    ## ----------------------chk.fns-----------------^
    return( cdb )
}
