#' Functions to wrap SQLite calls
#' 
#' sqliter helps users, mainly data munging practioneers, to organize
#' their sql calls in a clean structure. It simplifies the process of
#' extracting and transforming data into useful formats.
#'
#' @name sqliter-package
#' @docType package
#' @import stringr
#' @import DBI
#' @import RSQLite
#' @import functional
NULL

#' Creates the sqliter a kinf of SQLite database manager, but not that far.
#' 
#' \code{sqliter} object works pretty much like a database manager helping users to execute queries and transform data through a clean interface.
#' 
#' @export
#' @param ... arguments such as \code{path} must be provided during object instantiation.
#' @examples
#' \dontrun{DBM <- sqliter(path=c("data", "another/project/data"))}
#' 
sqliter <- function(...) {
    defaults <- list(...)
    
    get <- function(name, drop=TRUE) {
        if (missing(name))
            defaults
        else {
            if (drop && length(name) == 1)
                defaults[[name]]
            else
                defaults[name]
        }
    }
    set <- function(...) {
        dots <- list(...)
        if (length(dots) == 0) return()
        if (is.null(names(dots)) && length(dots) == 1 && is.list(dots[[1]]))
            if (length(dots <- dots[[1]]) == 0) return()
        defaults <<- merge(dots)
        invisible(NULL)
    }
    
    this <- list(get=get, set=set)
    class(this) <- 'sqliter'
    this
}

#' returns the paths of the given database
#' 
#' @param object \code{sqliter} object
#' @param database the SQLite database filename without extension
#' @export
#' @examples
#' \dontrun{
#' DBM <- sqliter(path=c("data", "another/project/data"))
#' find_database(DBM, "dummydatabase")
#' # "data/dummydatabase.db"
#' }
find_database <- function(object, database) UseMethod('find_database', object)

#' @rdname find_database
#' @method find_database sqliter
#' @S3method find_database sqliter
find_database.sqliter <- function(object, database) {
    path <- paste0(object$get('path'), '/', database, '.db')
    res <- sapply(path, file.exists)
    path[res]
}

#' execute query into a given database
#' 
#' Once you have a \code{sqliter} database properly set you can execute queries into that database and get your data transformed.
#' By default this function returns a data.frame object, but if you transform your data you can get whatever you need.
#' 
#' @param object \code{sqliter} object
#' @param database the SQLite database filename without extension
#' @param query the query string
#' @param post_proc a function to transform data, it receives a database and returns whatever you need.
#' @param ... additional arguments used by prepared queries
#' @export
#' @examples
#' \dontrun{
#' DBM <- sqliter(path=c("data", "another/project/data"))
#' ds <- execute(DBM, "dummydatabase", "select count(*) from dummytable")
#' ds <- execute(DBM, "dummydatabase", "select * from dummytable where 
#'       name = :name", name=c("Macunamima", "Borba Gato"))
#' ds <- execute(DBM, "dummydatabase", "select * from dummytable where
#'       name = :name", name=c("Macunamima", "Borba Gato"),
#'         post_proc=function(ds) {
#' ds <- transform(ds, birthday=as.Date(birthday))
#' ds
#' })
#' }
execute <- function(object, ...) UseMethod('execute', object)

#' @rdname execute
#' @method execute sqliter
#' @S3method execute sqliter
execute.sqliter <- function(object, database, query, post_proc=identity, ...) {
    path <- find_database(object, database)
    stopifnot(length(path) == 1)
    conn <- dbConnect('SQLite', path)
    if (length(list(...)) != 0) {
        ds <- dbGetPreparedQuery(conn, query, data.frame(...))
    } else {
        ds <- dbGetQuery(conn, query)
    }
    dbDisconnect(conn)
    post_proc(ds)
}

#' @method $ sqliter
#' @S3method $ sqliter
'$.sqliter' <- function(object, name) {
    if (str_detect(name, "^query_(.*)$")) {
        database <- unlist(str_split_fixed(name, "_", 2))[2]
        Curry(execute, object, database)
    } else {
        object[[name]]
    }
}

#' query functions
#' 
#' **query functions** are dynamic functions which connect to a database, execute queries in it and transform data.
#' Actually it is a decorator for \code{execute} function.
#' \code{execute} has 5 arguments.
#' The first argument is an instance of the \code{sqliter} class and the second is the database name.
#' The call to a query function is executed like a method call to the \code{sqliter} object through the \code{$} operator.
#' The function name must have the following pattern: \code{query_<database name without extension>}.
#' This call returns an \code{execute} function with the first 2 argument already set.
#' The first parameter is the \code{sqliter} object on which the \code{$} operator have been called and the second argument is extracted from the query function name, the name after the preffix \code{query_}.
#' 
#' @name query-functions
#' @examples
#' \dontrun{
#' DBM <- sqliter(path=c("data", "another/project/data"))
#' DBM$query_dummydatabase("select count(*) from dummytable")
#' }
NULL