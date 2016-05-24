#' Returns all databases on the server
#' 
#' Gives a list of all databases available at \code{cdb$serverName}.
#' 
#' The function uses the \code{_all_dbs}  API end point . 
#' 
#' @usage cdbListDB(cdb)
#' @param cdb Only the connection settings \code{cdb$port} and
#' \code{cdb$serverName} is needed.
#' @return
#' 
#' \item{cdb }{The result of the request is stored in cdb$res after converting
#' the json answer into a list using \code{cdb$fromJSON()}.}
#' @author wactbprot
#' @export
#' @examples
#' \dontrun{
#' cdbListDB(cdbIni())$res
#'}
#' @seealso \code{cdbMakeDB}
#' @keywords misc

cdbListDB <- function(cdb){
  
  fname <- deparse(match.call()[[1]])
  cdb   <- cdb$checkCdb(cdb,fname)
  
  if(cdb$error == ""){
    adrString <- paste(cdb$baseUrl(cdb), "_all_dbs",
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
