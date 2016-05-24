#' Function for request one id
#'
#' Function returns a 128bit uuid requested from CouchDB
#'
#' Simple CouchDB API end point to http://serverName:port/_uuids.
#'
#' @param cdb Only the connection settings \code{cdb$port} and
#' \code{cdb$serverName} is needed.
#' @return
#'
#' \item{cdb }{The result of the request is stored in \code{cdb$id} after
#' converting the answer into a list using \code{fromJSON()}.}
#' @author wactbprot
#' @export
#' @examples
#' \dontrun{
#' cdbGetUuid(cdbIni())$res
#'}
#'
#' @seealso \code{cdbMakeDB}
#' @keywords misc

cdbGetUuid <- function(cdb){

    fname <- deparse(match.call()[[1]])
    cdb   <- cdb$checkCdb(cdb,fname)

    if(cdb$error == ""){

        adrString <- paste(cdb$baseUrl(cdb),
                           "_uuids?count=1",
                           sep = "")
        
        res           <- getURL(utils::URLencode(adrString),
                                customrequest = "GET",
                                curl          = cdb$curl,
                                .opts         = cdb$opts(cdb))

        cdb <-  cdb$checkRes(cdb,res)

        if(cdb$error == ""){
            cdb$id <- unlist(cdb$res$uuids)
            return(cdb)
        }else{
            stop(cdb$error)
        }
    }else{
        stop(cdb$error)
    }
}
