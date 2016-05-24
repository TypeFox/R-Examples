
## collects all elements in a chunk file and returns a named list of those elements
DGather <- function( x, keys = FALSE, n = -1L, names = TRUE ){
    chunks <- .get_chunks_from_current_revision(x)
    if( n > 0L ){
        len <- length(chunks)
        n <- ifelse( n > len, len, n )
        chunks <- chunks[1L:n] ## utils.R
    }

    out <- lapply(chunks,
                  function(f){
                      lines <- DS_read_lines( DL_storage(x),
                                     f )
                      ## note, the last line just contains information about the keys
                      neof <- grep("^<<EOF-", lines, invert = TRUE )
                      lapply(lines[ neof ], function(line)
                             if( !keys )
                             DSL_unserialize_object( strsplit( line, "\t" )[[ 1L ]][ 2L ] )
                             else
                             strsplit( line, "\t" )[[ 1L ]][ 1L ]
                             )
                  })
    if( names )
        names(out) <- chunks
    out
}

as.list.DList <- function(x, ...)
    structure( do.call(c, DGather(x)), names = names(x) )
