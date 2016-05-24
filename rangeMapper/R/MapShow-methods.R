
#' Method tables
#'
#' @param object an sqlite connection object.
#' @rdname tables-methods
#' @exportMethod tables

setGeneric("tables", function(object)       standardGeneric("tables") )

#' Tables and column names  of an sqlite db
#'
#' @param object an sqlite connection object.
#' @return data.frame
#' @export
#'
setMethod("tables",
    signature  = "SQLiteConnection",
        definition = function(object) {

        tbs = dbGetQuery(object, "select name from sqlite_master where type = 'table' ")$name

        X = lapply(tbs, function(x) dbGetQuery(object, paste0("PRAGMA table_info(", shQuote(x), ")") )  )

        mapply(FUN = function(X,tbs) { X$tabname =tbs; X } , X = X, tbs = tbs, SIMPLIFY = FALSE) %>%
        do.call(rbind, .)

        }
    )


#' rangeMap file info
#'
#' @param con a connection to a rangeMapper project.
#'
#' @return data.table
#' @export
#'
rangeMapProjInfo <- function(con) {

    tabs = tables(con)
    tabs = tabs[!tabs$tabnam %in% c('canvas', 'ranges') & tabs$name != 'id', c('name', 'tabname', 'type')]

    tabs$sql = paste('select', tabs$name, 'from', tabs$tabnam, 'limit 1')

    tabs$firstValue = sapply(tabs$sql, FUN = dbGetQuery, conn = con)

    tabs$sql = paste('select count(', tabs$name, ') from', tabs$tabnam)

    tabs$nrows = sapply(tabs$sql, FUN = dbGetQuery, conn = con)

    tabs$sql  = NULL

     tabs

}



setMethod("show", signature(object = "rangeMap"), function(object){
    out = list()
    tbs = dbGetQuery(object@CON, "select name from sqlite_master where type = 'table' ")$name

    out[["class"]] = class(object)
    dbinfo = dbGetInfo(object@CON)
    out[["Project_location"]] = object@CON@dbname
    out[["SQLite_version"]]   = dbinfo$serverVersion

    if( is.empty(object@CON, object@CANVAS) )
    out[["empty_project"]]  = "Empty rangeMapper project." else {
        out[["Proj4"]]      = dbReadTable(object@CON, object@PROJ4STRING)[1,1]
        out[["CellSize"]]   = dbReadTable(object@CON, object@GRIDSIZE)[1,1] %>% prettyNum
        Extent              = as.list(dbReadTable(object@CON, object@BBOX))
        out[["Extent"]]     = Extent  %>% prettyNum %>%
                              paste(names(Extent), ., sep  = '=') %>%
                              paste(.,collapse =',')
        out[["BIO_tables"]] = paste(  gsub(object@BIO, "", tbs[grep(object@BIO, tbs)]), collapse = ";" )
        out[["MAP_tables"]] = paste(  gsub(object@MAP, "", tbs[grep(object@MAP, tbs)]), collapse = ";" )
        mtd = is.empty(object@CON, object@METADATA_RANGES)

        out[["METADATA_RANGES"]] = dbGetQuery(object@CON, paste("select * from", object@METADATA_RANGES, "limit 1") ) %>%
                                    names  %>% setdiff(., 'bioid') %>%
                                    paste(.,collapse =',')

      }

     names(out)  = sprintf(paste0('%-', max(nchar(names(out))), 's'), names(out))


     delim = paste0('+', rep('_', nchar(names(out)[1])-2) %>%
              paste0(., collapse = ''), '+',collapse = '')
     message(delim)
     paste(names(out), out) %>%
     paste(., collapse = "\n")  %>%
     message
     message(delim)

     })



