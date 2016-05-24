#' List of history in the Elliptic DB
#' 
#' List of unique history reflecting the bootstrap activities.
#'
#' @method history ecdb
#'
#' @param object an object of ecdb class.
#' 
#' @return list of history
#'
#' @keywords ecdb
#'
#' @author Stephen H-T. Lihn
#'
#' @export
#' 
#' @import RSQLite
#'
### <======================================================================>
"history.ecdb" <- function(object)
{
    conn <- object@conn
    sql <- "SELECT DISTINCT activity FROM HISTORY h ORDER BY activity"
    RSQLite::dbGetQuery(conn, sql)$activity
}
### <---------------------------------------------------------------------->
#' @rdname history.ecdb
setGeneric("history", function(object) standardGeneric("history"))
#' @rdname history.ecdb
setMethod("history", "ecdb", history.ecdb)
### <---------------------------------------------------------------------->
# This is an internal helper library for ecdb, primarily for bootstrap
#
### <---------------------------------------------------------------------->
# insert activity to history
setGeneric("history<-", function(object, value) standardGeneric("history<-"))

setMethod("history<-", "ecdb", function(object, value){
    conn <- object@conn
        
    sql <- "INSERT INTO HISTORY VALUES (:iso_dttm, :activity, :time_stamp)"
    history <- data.frame(iso_dttm=c(as.character(Sys.time())),
                          activity=c(value),
                          time_stamp=c(as.integer(Sys.time()))
                         )
    rs <- RSQLite::dbGetPreparedQuery(conn, sql, bind.data=history)
    ecdb.protectiveCommit(object)
    Sys.sleep(2) # to ensure uniqueness of iso_sttm
    invisible(object)
})
### <---------------------------------------------------------------------->
