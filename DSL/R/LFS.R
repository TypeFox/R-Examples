################################################################################
## DStorage: LFS
## The functions below implement DList operations for the "Local File System"
################################################################################

.DMap.LFS <- function( storage, x, MAP, parallel ){

    new_rev <- .make_DSL_revision()
    DS_dir_create(storage, new_rev)

    LAPPLY <- if( parallel ){
        parallel::mclapply
    } else {
        lapply
    }

    LAPPLY( basename(.get_chunks(x)), function( chunk )
            .LFS_mapper( MAP = MAP,
                         input  = file.path(storage$base_directory, .revisions( x )[1], chunk),
                         output = file.path(storage$base_directory, new_rev, chunk)) )

    new_rev
}

.DReduce.LFS <- function( storage, x, REDUCE, parallel, ... ){

    intermed_rev <- .make_DSL_revision()
    DS_dir_create(storage, intermed_rev)

    LAPPLY <- if( parallel ){
        parallel::mclapply
    } else {
        lapply
    }

    outchunk <- basename(.get_chunks(x)[1])
    LAPPLY( basename(.get_chunks(x)), function( chunk )
            .LFS_reducer( REDUCE = REDUCE,
                         input  = file.path(storage$base_directory, .revisions( x )[1], chunk),
                         output = file.path(storage$base_directory, intermed_rev, outchunk)),
            ...
         )

    if(length(.get_chunks(x)) > 1L){
        new_rev <- .make_DSL_revision()
        DS_dir_create(storage, new_rev)
        .LFS_reducer( REDUCE = REDUCE,
                      input  = file.path(storage$base_directory, intermed_rev, outchunk),
                      output = file.path(storage$base_directory, new_rev, outchunk) )
    } else {
        new_rev <- intermed_rev
    }
    new_rev
}

.LFS_mapper <- function( MAP, input, output ){

    con <- file( input, open = "r" )
    con2 <- file( output, open = "w" )

    while (length(line <- readLines(con, n = 1L, warn = FALSE)) > 0) {
        input <- DSL_split_line( line )
        if( length(grep("^<<EOF-", input$key)) ){
            chunk <- as.character(input$value["Chunk"])
            break
        }

        result <- MAP( input )

        ## FIXME: should be an object oriented approach here
        ##        there should be a test if given structure is really a (list of) key value pair(s).
        if( (length(result) == 2) && identical(names(result), c("key", "value")))
            writeLines( sprintf("%s\t%s", result$key, DSL_serialize_object(result$value)), con2 )
        else
            for( i in seq_along(result) )
                writeLines( sprintf("%s\t%s", result[[i]]$key, DSL_serialize_object(result[[i]]$value)), con2 )
    }
    ## In the last step we need to add a stamp to this chunk
    ## <key:randomstring, value_serialized:c(firstdocumentkey,lastdocumentkey)>
    writeLines( .make_chunk_signature( chunk ),
                con2 )
    close(con)
    close(con2)

    invisible( TRUE )
}

.LFS_reducer <- function( REDUCE, input, output){

    ## INIT

    ## use efficient collector for integer pairlists
    CONCATENATE <- function( collector = FALSE )
        if( collector ){
            .collector2 ## FIXME: untested remove of DSL:::
        } else {
            base::c
        }

    chunk <- NA
    INTPAIRLIST <- NULL

    ## initialize hash table holding reduce results
    env <- new.env( hash = TRUE, size = 10240 )

    ## CON

    con  <- file( input, open = "r" )
    con2 <- file( output, open = "at" )

    ## STREAM

    while (length(line <- readLines(con, n = 1L, warn = FALSE)) > 0) {
        input <- DSL_split_line( line ) ## FIXME: untested remove of DSL:::
        ## Skip end of line
        if( length(grep("^<<EOF-", input$key)) ){
            chunk <- as.character(input$value["Chunk"])
            next
        }

        ## we have an efficient collector for integer pair lists (based on linked lists)
        if( is.null(INTPAIRLIST) )
            INTPAIRLIST <- is.list(input$value) && all(unlist(lapply(input$value, is.integer)))

        tryCatch( assign(input$key,
                         CONCATENATE(INTPAIRLIST)(if(tryCatch(exists(input$key, envir = env, inherits = FALSE),
                                                              error = function(x) FALSE))
                                                  get(input$key, envir = env, inherits = FALSE)
                         else
                                                  NULL,
                                                  input$value
                                                  ),
                         envir = env
                         ), error = function(x) FALSE )
    }

    ## CLOSE
    close(con)

    ## OUTPUT
    env <- as.list(env)
    if( INTPAIRLIST ){
        env <- lapply(env, .collector2, NULL) ## FIXME: untested remove of DSL:::
    }
    keys <- names(env)
    for( i in seq_along(keys) )
        writeLines( sprintf("%s\t%s", keys[i],
                     DSL_serialize_object(REDUCE(env[[ i ]]))), con = con2 ) ## FIXME: untested remove of DSL:::
    writeLines( .make_chunk_signature( chunk ),
               con2 )
    close(con2)
}
