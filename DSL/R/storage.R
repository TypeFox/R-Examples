################################################################################
## DStorage: the distributed storage class
## All functions dealing with the storage implementation
## of the DList class are provided here.
## Each storage is of general class 'DStorage' and of a subclass defining
## the actual storage, e.g. 'LFS', 'HDFS'
################################################################################


################################################################################
## DStorage (DS) high level constructor
################################################################################

DStorage <- function( type = c("LFS", "HDFS"),
                      base_dir,
                      chunksize = 1024^2 ){
    type <- match.arg( type )
    .DS_init( switch(type,
                     "LFS"  = .make_LFS_storage(base_dir, chunksize),
                     "HDFS" = .make_HDFS_storage(base_dir, chunksize)) )
}

################################################################################
## Default storage (local disk)
################################################################################

DS_default <- function(){
    .DS_init( .make_LFS_storage(base_dir  = tempdir(),
                                chunksize = 10 * 1024^2) )
}

################################################################################
## DS_storage low level constructor
################################################################################

.DStorage <- function( type,
                       base_directory,
                       chunksize,
                       dir_create,
                       fetch_last_line,
                       get,
                       list_directory,
                       put,
                       read_lines,
                       unlink,
                       write_lines ){
    structure( list(base_directory = base_directory,
                    chunksize = chunksize,
                    dir_create = dir_create,
                    fetch_last_line = fetch_last_line,
                    get = get,
                    list_directory = list_directory,
                    put = put,
                    read_lines = read_lines,
                    unlink = unlink,
                    write_lines = write_lines),
              class = c(type, "DStorage"))
}

.make_LFS_storage <- function( base_dir, chunksize, ... ){
    .DStorage( type = "LFS",
               base_directory  = base_dir,
               chunksize       = chunksize,
               dir_create      = function(x) base::dir.create(x, showWarnings = FALSE),
               fetch_last_line = function(x) utils::tail(base::readLines(as.character(x)), n = 1L),
               get             = function(x) {con <- file(x, open = "r")
                                              obj <- unserialize(con)
                                              close.connection(con)
                                              obj},
               list_directory  = base::dir,
               put             = function( obj, x ) {con <- file( x, open = "w" )
                                                     serialize( obj, con )
                                                     close.connection(con)
                                                 },
               read_lines      = function(x) base::readLines(as.character(x)),
               unlink          = function(x) LFS_remove(x), ## removed DSL:::
               write_lines     = function(text, fil) base::writeLines(text, con = as.character(fil)), ...
              )
}

.make_HDFS_storage <- function( base_dir, chunksize, ... ){
    .DStorage( type = "HDFS",
               base_directory  = base_dir,
               chunksize       = chunksize,
               dir_create      = function(x) hive::DFS_dir_create(x, henv = hive::hive()),
               fetch_last_line = function(x) hive::DFS_tail(n = 1L, as.character(x), henv = hive::hive()),
               get             = function(x) hive::DFS_get_object(as.character(x), henv = hive::hive()),
               list_directory  = function(x) hive::DFS_list(x, henv = hive::hive()),
               put             = function(obj, x) hive::DFS_put_object(obj, as.character(x), henv = hive::hive()),
               read_lines      = function(x) hive::DFS_read_lines(as.character(x), henv = hive::hive()),
               unlink          = function(x) hive::DFS_delete(x, recursive = TRUE, henv = hive::hive()),
               write_lines     = function(text, fil) hive::DFS_write_lines(text, as.character(fil), henv = hive::hive()), ...
              )
}

################################################################################
## DStorage S3 methods
################################################################################

is.DStorage <- function( ds )
    inherits( ds, "DStorage" )

as.DStorage <- function( x )
    UseMethod("as.DStorage")

as.DStorage.DStorage <- identity

print.DStorage <- function( x, ... ){
    writeLines( "DStorage." )
    writeLines( sprintf("- Type: %s", class(x)[1] ) )
    writeLines( sprintf("- Base directory on storage: %s", DS_base_dir(x)) )
    writeLines( sprintf("- Current chunk size [bytes]: %s", DS_chunksize(x)) )
}

summary.DStorage <- function( object, ... ){
    print( object )
    meths <- names(object)[!(names(object) %in% c("description", "chunksize", "base_directory"))]
    lenpermeth <- ceiling(length(meths)/2)
    writeLines( sprintf("- Registered methods:\n  %s\n  %s",
                        paste(meths[1:lenpermeth], collapse = ", "),
                        paste(meths[(lenpermeth+1):length(meths)], collapse = ", ")) )
}

################################################################################
## DStorage helper functions
################################################################################

## helper function for DList -> list coercion
## checks if all chunks are already in place
## FIXME: we need to provide some checksums here
.check_contents_of_storage <- function(x, value){
    all( basename(.get_chunks_from_current_revision(x))	 %in%
        DS_list_directory(value, .revisions(x)[1]) )
}



.DS_empty <- function(){
    NA
}

## .storage_init() initializes the storage to be used for the distributed corpus.
.DS_init <- function( x ) {
    if( missing(x) )
        x <- DS_default()
    ## Create storage base directory
    if( is.DStorage(x) )
       x$dir_create( DS_base_dir(x) )

    out <- if( inherits(tryCatch(check_storage_for_sanity(x),
                                 error = identity), "error") )
        .DS_empty()
    else
        x
    out
}

DS_base_dir <- function( storage )
    storage$base_directory

DS_chunksize <-function( storage )
    storage$chunksize

DS_dir_create <- function( storage, dir ){
    stopifnot( is.DStorage(storage) )
    dir <- sub("/./", "", file.path(DS_base_dir(storage), dir ), fixed = TRUE )
    storage$dir_create( dir )
}


DS_fetch_last_line <- function( storage, file ){
    stopifnot( is.DStorage(storage) )
    file <- sub("/./", "", file.path(DS_base_dir(storage), file ), fixed = TRUE)
    storage$fetch_last_line( file )
}

DS_get <-function( storage, x ){
    stopifnot( is.DStorage(storage) )
    storage$get( file.path(DS_base_dir(storage), x) )
}

DS_list_directory <- function( storage, dir = "." ){
    stopifnot( is.DStorage(storage) )
    if( dir == "." )
        dir <- ""
    storage$list_directory( file.path(DS_base_dir(storage), dir ) )
}

DS_put <-function( storage, obj, x ){
    stopifnot( is.DStorage(storage) )
    storage$put( obj, file.path(DS_base_dir(storage), x) )
}

DS_read_lines <- function( storage, file ){
    stopifnot( is.DStorage(storage) )
    #file <- sub("/./", "", file.path(storage$base_directory, file ))
    file <- file.path(DS_base_dir(storage), file )
    storage$read_lines( file )
}

## NOTE: due to security reasons no recursive unlinking is permitted!
## deletes a link on the corresponding storage (file, directory)
DS_unlink <- function( storage, link ){
    stopifnot( is.DStorage(storage) )
    link <- sub("/./", "", file.path(DS_base_dir(storage), link ), fixed = TRUE)
    storage$unlink( link )
}

DS_write_lines <- function( storage, text, file ){
    stopifnot( is.DStorage(storage) )
    file <- sub("/./", "", file.path(DS_base_dir(storage), file ), fixed = TRUE)
    storage$write_lines( text, file )
}

check_storage_for_sanity <- function( storage ){
    stopifnot( is.DStorage(storage) )

    dir <- basename( tempfile() )
    DS_dir_create( storage, dir )

    contents <- DS_list_directory( storage )
    stopifnot(any(contents == dir))

    stopifnot(DS_unlink( storage, dir ))

    text <- c("some text", "other text", "more text")
    file <- "some_file"
    DS_write_lines( storage, text, file )
    text_read <- DS_read_lines( storage, file )
    stopifnot( all(text == text_read) )

    line <- DS_fetch_last_line( storage, file )
    stopifnot( line == text[length(text)] )

    stopifnot( DS_unlink( storage, file ) )
    invisible( TRUE )
}

################################################################################
## local file system manipulators
################################################################################

LFS_remove <- function( x ){
    foo <- file.remove
    ## on Windows file.remove does not remove empty directories
    ## so we need to check for it
    if( .Platform$OS.type == "windows" )
        if( length(dir(x, recursive = TRUE, all.files = TRUE)) )
            warning( sprintf("Directory %s not empty, so not removed.", x) )
        else
            foo <- function( x ) {
                unlink( x, recursive = TRUE )
                TRUE
            }
    foo( x )
}
