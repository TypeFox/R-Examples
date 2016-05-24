#' This function adds multiple database documents with one request
#'
#' This is done via the _bulk_docs API provided by
#' an already existing database.
#'
#' The _bulk_docs endpoint requires that \code{cdb$dataList} resolves
#' to an json array. This is reached with e.g.
#' \code{cdb$dataList <- list(list(...),list(...),...)}.
#' Furthermore, _bulk_docs requires the documents to be wrapped in a key
#' named \code{docs:[...]}; this is done by \code{cdbAddDocS()} if
#' \code{cdb$dataList} is a list of lists. The user dont need to care.
#'
#' At the moment the resulting \code{_rev} and  \code{_id} will be not
#' written back to the \code{cdb$dataList}. This means that a second
#' call of  \code{cdbAddDocS()} generates new Documents. 
#'
#' @usage cdbAddDocS(cdb)
#' @param cdb \code{cdb$dataList} has to be a list of lists,
#' \code{cdb$DBName}, \code{cdb$serverName} is needed.
#' @return \item{cdb}{The couchdb response is stored in \code{cdb$res} }
#' @author parisni, wactbprot
#' @export
#' @examples
#' \dontrun{
#' ccc               <- cdbIni()
#' # I assume a database at localhost:5984 already exists
#' ccc$DBName        <- "r4couchdb_db"
#' docs <- list()
#' for(i in 1:10){
#'  docs[[i]] <- list(normalDistRand =  rnorm(20))
#' }
#' # docs is noe a list of 10 lists
#' ccc$dataList <- docs
#' # generating 10 database documents
#' cccAddDocS(ccc)$res
#' }
#' @seealso \code{cdbAddDoc()}
#' @keywords misc
#'

cdbAddDocS <- function( cdb){
    
    fname <- deparse(match.call()[[1]])
    cdb   <- cdb$checkCdb(cdb,fname)
    
    if(cdb$error == ""){
        adrString <- paste(cdb$baseUrl(cdb),
                           cdb$DBName,"/",
                           "_bulk_docs",
                           sep="")
        docs  <- list(docs = cdb$dataList)
        
        res   <- getURL(utils::URLencode(adrString),
                        customrequest   = 'POST',
                        postfields      = cdb$toJSON(docs),
                        httpheader      = c('Content-Type: application/json;charset=utf-8'),
                        .opts           = cdb$opts(cdb),
                        curl            = cdb$curl)
        res <- cdb$fromJSON( res )
        if(length(res$error) == 0){
            cdb$res <- res
        return( cdb )
        }else{
            cdb$error <- paste(cdb$error, res$error)
    }
}
    
    if(!(cdb$error == "")){
        stop( cdb$error )
    }
}
