################################################################################
## This file defines the reduce function on DList objects
## DReduce() is conceptually based on Google's MapReduce paradigm
################################################################################

################################################################################
## DReduce function
################################################################################

DReduce <- function( x, REDUCE = identity, parallel, ... ){
    x <- as.DList( x )
    if( missing(parallel) )
        parallel <- FALSE
    ## HDFS is always parallel since we cannot easily control parallel
    ## execution on Hadoop clusters
    if( inherits(DL_storage(x), "HDFS") )
        parallel <- TRUE

    new_rev <- .DReduce( DL_storage(x), x = x, REDUCE = REDUCE, parallel = parallel, ... )

    ## there are possibly fewer chunks after reduce (check for them)
    cn <- basename(.get_chunks(x))
    chunks <- cn[cn %in% DS_list_directory(DL_storage(x), new_rev)]
    out <- .DList( list(),
                  .make_chunk_handler(file.path(new_rev, chunks),
                                      new_rev,
                                      DL_storage(x)),
                  attr( x, "Keys" ),
                  attr( x, "Mapping" ),
                  DL_storage( x )
                  )
    gathered <- DGather( out, keys = TRUE, names = FALSE )
    keys <- unlist( gathered )
    ## reconstruct mapping if we have more keys after the map step
    if( length(keys) != dim(attr(out, "Mapping" ))[1] ){
        len <- unlist( lapply( gathered, length ) )
        new_mapping <- cbind( rep.int(seq_along(len), len), unlist(lapply(len, seq_len)))
        colnames( new_mapping ) <- c( "Chunk", "Position" )
        rownames( new_mapping ) <- keys
        attr( out, "Mapping" ) <- new_mapping
       }
    attr( out, "Keys" ) <- keys
    out
}

################################################################################
## .DReduce methods (depend on storage type)
################################################################################

.DReduce <- function( storage, x, REDUCE, parallel, ... )
    UseMethod( ".DReduce" )
