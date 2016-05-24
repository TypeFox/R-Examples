#' Receive list results from CouchDB
#' 
#' The function provides accesses to CouchDB lists.
#' 
#' Query params e.g. \code{"reduce=false"} or \code{"group_level=1"} can be
#' provided in \code{cdb$queryParam} By now multible params must be given in
#' one string e.g.  \code{"a=b&c=d&e=f"}.
#' 
#' @usage cdbGetList(cdb)
#' @param cdb Beside the connection details (\code{cdb$port},\code{cdb$DAName}
#' ...) the \code{cdb$design} \code{cdb$view} and \code{cdb$list} is needed.
#' @return
#' 
#' \item{cdb }{The result of the request is stored in cdb$res after converting
#' the json answer into a list using fromJSON(). If a needed cdb (design, list,
#' view, ...) entry was not provided cdb$error says something about the R side.
#' 
#' }
#' @author wactbprot
#' @export
#' @keywords misc
#'

cdbGetList <- function( cdb ){
  
  fname <- deparse(match.call()[[1]])
  cdb   <- cdb$checkCdb(cdb,fname)
  
  if(cdb$error ==""){

    if(cdb$queryParam == ""){
      queryString <- ""
    }else{
      ## todo: let usr specify list(a=b, c=d) ...
      queryString <- paste("?",cdb$queryParam, sep="")
    }

    adrString <- paste(cdb$baseUrl(cdb),
                       cdb$DBName,
                       "/_design/",
                       cdb$design,
                       "/_list/",
                       cdb$list,
                       "/",cdb$view,
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
