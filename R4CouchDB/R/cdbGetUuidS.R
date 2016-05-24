#' Function for request some ids
#' 
#' Function returns a 128bit uuid requested from CouchDB
#' 
#' CouchDB API provides the url http://serverName:port/_uuids for those clients
#' who aren't able to create those ids. The number N
#' of ids received from a CouchDB can be set by \code{cdb$count <- N}
#' since version 0.6.  The function writes to cdb$res (in opposite to
#' \code{cdbGetUuid()} whitch writes to \code{cdb$id})
#' 
#' @param cdb Only the connection settings \code{cdb$port},
#' \code{cdb$serverName} and \code{cdb$count} is needed.
#' @return
#' 
#' \item{cdb }{The result of the request is stored in \code{cdb$res} after
#' converting the answer into a list using \code{fromJSON()}.}
#' @author wactbprot
#' @export
#' @examples
#'\dontrun{
#' ccc            <- cdbIni()
#' ccc$count      <- 100
#' cdbGetUuidS(ccc)$res
#'}
#' @seealso \code{cdbMakeDB}
#' @keywords misc
#'

cdbGetUuidS <- function(cdb){

    fname <- deparse(match.call()[[1]])
    cdb   <- cdb$checkCdb(cdb,fname)
    
    if (cdb$error == ""){
        
        adrString     <- paste(cdb$baseUrl(cdb),
                               "_uuids?count=",cdb$count,
                               sep="")
        
        res <- getURL(utils::URLencode(adrString),
                      customrequest = "GET",
                      curl          = cdb$curl,
                      .opts         = cdb$opts(cdb))
        
        return(cdb$checkRes(cdb,res))
        
    }else{
        stop(cdb$error)
  }
}
