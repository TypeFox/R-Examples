#' Receive view results from CouchDB
#' 
#' The function provides accesses to CouchDB views.
#' 
#' Query params e.g. \code{"reduce=false"} or \code{"group_level=1"} can be
#' provided in \code{cdb$queryParam}
#' 
#' @usage cdbGetView(cdb)
#' @param cdb Beside the connection details (\code{cdb$port},\code{cdb$DAName}
#' ...) the \code{cdb$design} and \code{cdb$view} is needed.
#' @return
#' 
#' \item{cdb }{The result of the request is stored in cdb$res after converting
#' the json answer into a list using fromJSON(). If a needed cdb list entry was
#' not provided cdb$error says something about the R side
#' 
#' }
#' @note For the moment only one \code{cdb$queryParam} is possible. In the
#' future maybe Duncans \code{RJavaScript} package can be used to generate
#' views without leaving R.
#' @author wactbprot
#' @export
#' @keywords misc

cdbGetView <- function( cdb ){
  
  fname <- deparse(match.call()[[1]])
  cdb   <- cdb$checkCdb(cdb,fname)
  
  if(cdb$error ==""){
    
    if(cdb$queryParam == ""){
      queryString <- ""
    }else{
      queryString <- paste("?",cdb$queryParam, sep="")
    }
    
    adrString <- paste(cdb$baseUrl(cdb),
                       cdb$DBName,
                       "/_design/",
                       cdb$design,
                       "/_view/",
                       cdb$view,
                       queryString,
                       sep="")
    
    res <- getURL(utils::URLencode(adrString),
                  customrequest = "GET",
                  curl          = cdb$curl,
                  .opts         = cdb$opts(cdb))
    
    return(cdb$checkRes(cdb,res))
    
  }else{
    stop( cdb$error )
  }
}
