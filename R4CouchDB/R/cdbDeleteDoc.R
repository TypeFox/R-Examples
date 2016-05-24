#' Deletes a document from a database
#' 
#' With a given \code{cdb$id} this function sends a http \code{"DELETE"}
#' request to the url \code{.../cdb$id?rev=cdb$rev}.
#' 
#' 
#' @usage cdbDeleteDoc(cdb)
#' @param cdb Beside \code{cdb$serverName}, \code{cdb$port} and
#' \code{cdb$DBName} the \code{cdb$id} must be given. R errors are reported in
#' \code{cdb$errors}
#' @return
#' 
#' \item{cdb }{The result of the delete request is stored in
#' \code{cdb$res}(whatever this means).  }
#' @author wactbprot
#' @export
#' @seealso \code{cdbAddDoc()}
#' @keywords misc
#'

cdbDeleteDoc <- function( cdb ){
  
  fname <- deparse(match.call()[[1]])
  cdb   <- cdb$checkCdb(cdb,fname)

  if(cdb$error == ""){
    cdb     <- cdbGetDoc(cdb)
    cdb$rev <- cdb$res$'_rev'
    
    adrString <- paste(cdb$baseUrl(cdb),
                       cdb$DBName,"/",
                       cdb$id,
                       "?rev=",
                       cdb$rev,
                       sep="")
    
    res <- getURL(utils::URLencode(adrString),
                  customrequest = "DELETE",
                  curl          = cdb$curl,
                  .opts         = cdb$opts(cdb))
    
    return(cdb$checkRes(cdb,res))
    
  }else{
    stop(cdb$error)
  }
}
