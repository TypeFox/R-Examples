#' This function updates an existing doc
#'
#' This essentially means that a
#' revision, corresponding to the '_id' has to be provided. If no '_rev' is
#' given in the \code{cdb} list the function gets the doc from the db
#' and takes the rev number for the update
#'
#' Updating a doc at couchdb means executing a http "PUT" request. The
#' \code{cdb} list must contain the \code{cdb$serverName}, \code{cdb$port},
#' \code{cdb$DBName}, \code{cdb$id}. Since v0.6 the revision of the document
#' should exist at the intended place: \code{cdb$dataList$'_rev'}.
#'
#' \code{getURL()} with \code{customrequest = "PUT"} does the work.  If a
#' needed \code{cdb$} list entry is not provided \code{cdb$error} maybe says
#' something about the R side.
#'
#' @usage cdbUpdateDoc(cdb)
#' @param cdb the cdb connection configuration list must contain the
#' \code{cdb$serverName}, \code{cdb$port}, \code{cdb$DBName} and \code{cdb$id}.
#' The data which updates the data stored in the doc is provided in
#' \code{cdb$dataList}
#' @return \item{cdb }{The response of the request is stored in \code{cdb$res}
#' after converting the answer by means of \code{fromJSON()}. The revision
#' provided by the respons is used for updating the \code{cdb$dataList$'_rev'}.
#' }
#' @author wactbprot
#' @export
#' @examples
#' \dontrun{
#' ccc               <- cdbIni()
#' # I assume a database at localhost:5984 already exists
#' ccc$DBName        <- "r4couchdb_db"
#' ccc$dataList      <- list(normalDistRand =  rnorm(20))
#' ccc               <- cdbAddDoc(ccc)
#'
#' ccc$dataList$Date <- date()
#' ccc               <- cdbUpdateDoc(ccc)
#'}
#'
#' @seealso \code{cdbInit()}
#' @keywords misc
#'

cdbUpdateDoc <- function( cdb){

    fname <- deparse(match.call()[[1]])
    cdb   <- cdb$checkCdb(cdb,fname)


    rev <- cdb$getDocRev(cdb)
    if(!is.na(rev)){
        cdb$dataList[["_rev"]] <- rev
    }

    if(cdb$error == ""){

        adrString   <- paste(cdb$baseUrl(cdb),
                             cdb$DBName,"/",
                             cdb$id,
                             sep="")

        res <- getURL(utils::URLencode(adrString),
                      customrequest = "PUT",
                      postfields    = cdb$toJSON(cdb$dataList),
                      httpheader    = c('Content-Type: application/json;charset=utf-8'),
                      curl          = cdb$curl,
                      .opts         = cdb$opts(cdb))

        cdb <-  cdb$checkRes(cdb,res)

        if((length(cdb$res$ok)) > 0 ){
            cdb$dataList[["_rev"]] <- cdb$res$rev
            cdb$rev <- cdb$res$rev
        }
        return(cdb)
    }else{
        stop(cdb$error)
    }
}
