#' Summary for the Elliptic DB (ECDB)
#' 
#' Summary for the Elliptic DB (ECDB)
#'
#' @method summary ecdb
#'
#' @param object an object of ecdb class.
#' @param ... more arguments for summary. Currently not used.
#'
#' @keywords ecdb
#'
#' @author Stephen H-T. Lihn
#'
#' @export
#' 
#' @import RSQLite
#'
#' @examples
#' summary(ecdb())
### <======================================================================>
"summary.ecdb" <- function(object, ...)
{
    conn <- object@conn
    # row count
    cnt <- .ecdb.row_count(object)
    
    print(paste("ecdb file:", object@file, 
                ";", ceiling(file.info(object@file)$size/1000),"KB",
                ";", ceiling(cnt/1000),"K rows",
                ";", ifelse(object@is.internal,"Internal", "External")
                ))
    # version
    sql <- "SELECT version, note, iso_dttm 
            FROM VERSION v 
            ORDER BY time_stamp DESC"
    print(RSQLite::dbGetQuery(conn, sql))
    # history
    sql <- "SELECT iso_dttm, activity FROM HISTORY h ORDER BY time_stamp DESC"
    print(RSQLite::dbGetQuery(conn, sql))

    invisible(object)
}
### <---------------------------------------------------------------------->
#' @rdname summary.ecdb
setGeneric("summary", function(object, ...) standardGeneric("summary"))
#' @rdname summary.ecdb
setMethod("summary", signature("ecdb"), summary.ecdb)
### <---------------------------------------------------------------------->