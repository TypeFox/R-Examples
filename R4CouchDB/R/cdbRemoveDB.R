#' Function to remove a database
#'
#' Removing a database means sending a http- "DELETE"- request to
#' \code{http://cdb$serverName:cdb$port/ ...}
#'
#' In \code{cdb} a entry \code{cdb$delDBName} should be provided for more
#' explicit deleting respectively more secure removing.
#'
#' @usage cdbRemoveDB(cdb)
#' @param cdb The \code{cdb} has to provide \code{cdb$serverName},
#' \code{cdb$port} and \code{cdb$DBName}
#' @return \item{cdb}{The CouchDB answer is stored in \code{cdb$res}. Any
#' problems on the R side are reportet in \code{cdb$error} }
#' @author wactbprot
#' @export
#'@examples
#'\dontrun{
#' ccc               <- cdbIni()
#' ccc$newDBName     <- "r4couchdb_db"
#' ccc               <- cdbMakeDB(ccc)
#' ccc$res
#' ccc$removeDBName  <- ccc$DBName
#' cdbRemoveDB(ccc)$res
#'}
#'
#' @seealso \code{cdbMakeDB}
#' @keywords misc

cdbRemoveDB <- function(cdb){

  fname <- deparse(match.call()[[1]])
  cdb   <- cdb$checkCdb(cdb,fname)

  if( cdb$error == ""){
    adrString <- paste(cdb$baseUrl(cdb),
                       cdb$removeDBName,
                       sep="")

    res       <- getURL(utils::URLencode(adrString),
                        customrequest = "DELETE",
                        curl          = cdb$curl,
                        .opts         = cdb$opts(cdb))

    cdb$removeDBName <- ""

    return(cdb$checkRes(cdb,res))

  }else{
    stop(cdb$error)
  }
}
