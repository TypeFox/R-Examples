#' Get a doc from CouchDB
#' 
#' With a given \code{cdb$id} this function requests the document.
#' 
#' 
#' @usage cdbGetDoc(cdb)
#' @param cdb Beside \code{cdb$serverName}, \code{cdb$port} and
#' \code{cdb$DBName} the \code{cdb$id} must be given. R errors are reported
#' 
#' in cdb$errors
#' @return
#' \item{cdb }{The result of the request is stored in \code{cdb$res} after
#' converting the answer into a list using \code{fromJSON()}. If a list entry
#' needed in \code{cdb} is not provided \code{cdb$error} gives some
#' information.
#' 
#' }
#' @author wactbprot
#' @export
#' @examples
#' \dontrun{
#' ccc               <- cdbIni()
#' ccc$newDBName     <- "r4couchdb_db"
#' ccc$dataList      <- list(normalDistRand =  rnorm(20))
#' ccc               <- cdbAddDoc(ccc)
#' cdbGetDoc(ccc)$res
#' }
#' @seealso \code{cdbAddDoc()}
#' @keywords misc
#'

cdbGetDoc <- function(cdb){

  fname <- deparse(match.call()[[1]])
  cdb   <- cdb$checkCdb(cdb,fname)
    
  if(cdb$error == ""){
    adrString <- paste(cdb$baseUrl(cdb),
                       cdb$DBName,"/",
                       cdb$id,
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
