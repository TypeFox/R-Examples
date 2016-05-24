################################################################################
## This file defines the map/lapply function on DList objects
## DMap is conceptually based on Google's MapReduce paradigm
## DLapply() is based on DMap() and behaves like R's lapply
################################################################################

################################################################################
## DLapply function
################################################################################

DLapply <- function( x, FUN, parallel, ..., keep = FALSE ){
    args <- list(...)
    foo <- match.fun(FUN)
    MAP <- function( keypair )
        list( key = keypair$key, value = do.call(foo, c(list(keypair$value), args)) )
    DMap( x = x, MAP = MAP, parallel = parallel, keep = keep )

}

## (
## idea: DSL constructor might include a modifyable plain text
## reader to read a bunch of documents
## )

################################################################################
## DMap function
################################################################################

DMap <- function( x, MAP, parallel, keep = FALSE ){
    stopifnot( inherits(x, "DList") )
    if( missing(parallel) )
        parallel <- FALSE
    ## HDFS is always parallel since we cannot easily control parallel
    ## execution on Hadoop clusters
    if( inherits(DL_storage(x), "HDFS") )
        parallel <- TRUE
    new_rev <- .DMap( DL_storage(x), x = x, MAP = MAP, parallel = parallel )

    if( keep )
        #### FIXME: o currently this has some possibly unforeseen side effects
        ####        o garbage collection
        out <- .DList_add_revision( x, new_rev )
    else
        out <- .DList( list(),
                       .make_chunk_handler(file.path(new_rev, basename(.get_chunks( x ))),
                                           new_rev,
                                           DL_storage(x)),
                       attr( x, "Keys" ),
                       attr( x, "Mapping" ),
                       DL_storage( x )
                      )
    chunks <- DGather( out, keys = TRUE, names = FALSE )
    keys <- unlist( chunks )
    ## reconstruct mapping if we have more keys after the map step
    if( length(keys) > dim(attr(out, "Mapping" ))[1] ){
        len <- unlist( lapply( chunks, length ) )
        new_mapping <- cbind( rep.int(seq_along(len), len), unlist(lapply(len, seq_len)))
        colnames( new_mapping ) <- c( "Chunk", "Position" )
        rownames( new_mapping ) <- keys
        attr( out, "Mapping" ) <- new_mapping
       }
    attr( out, "Keys" ) <- keys

    out
}

################################################################################
## .DMap methods (depend on storage type)
################################################################################

.DMap <- function( storage, x, MAP, parallel )
    UseMethod( ".DMap" )

## not yet used:
DPair <- function( key, value )
    list( key = make.names(key), value = value )

.get_keys_from_current_revision <- function( x ){
    structure( unlist(DGather( x, keys = TRUE)), names = NULL )
}

## updates given list with new revision
.DList_add_revision <- function( x, rev ){
    ## add new revision
    chunks <- .get_chunks_from_current_revision( x )
    assign(rev, file.path(rev, basename(chunks)), envir = attr(x, "Chunks"))
    .revisions( x ) <- c( rev, .revisions(x) )
    ## update to active revision
    x
}
