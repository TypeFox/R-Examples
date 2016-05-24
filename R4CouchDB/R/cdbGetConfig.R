#' Request couchdb config
#' 
#' Function provides access to the \code{_config} api end point.
#' 
#' 
#' @param cdb Only the connection settings \code{cdb$port} and
#' \code{cdb$serverName} is needed.
#' @return
#' 
#' \item{cdb }{The result of the request is stored in \code{cdb$res} after
#' converting the answer into a list using \code{fromJSON()}.}
#' @author wactbprot
#' @export
#' @examples
#'\dontrun{
#' cdbGetConfig(cdbIni())$res
#'}
#' @seealso \code{cdbMakeDB}
#' @keywords misc

cdbGetConfig <- function(cdb){
  
  fname <- deparse(match.call()[[1]])
  cdb   <- cdb$checkCdb(cdb,fname)
  
  if (cdb$error ==""){

    adrString <- paste(cdb$baseUrl(cdb),
                       "_config",
                       sep="")

    res       <- getURL(utils::URLencode(adrString),
                        customrequest = "GET",
                        .opts         = cdb$opts(cdb),
                        curl          = cdb$curl)
    
    return(cdb$checkRes(cdb,res))
  }else{
    stop( cdb$error )
  }
}
