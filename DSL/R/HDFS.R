################################################################################
## DStorage: HDFS
## The functions below implement DList operations for the
## "Hadoop Distributed File System"
################################################################################

.DMap.HDFS <- function( storage, x, MAP, parallel ){

    ## SET a user specific command environment variable here
    cmdenv_arg <- NULL
    cmdenv_arg <- c( cmdenv_arg,
                     sprintf("_HIVE_FUNCTION_TO_APPLY_=%s",
                             DSL_serialize_object(MAP)) )

    rev <- .MapReduce( x, .HDFS_mapper(), cmdenv_arg = cmdenv_arg  )

    ## now handle part-xxxx stuff
    ## read all chunk signatures
    .fix_chunknames( DL_storage(x), rev )

    rev
}

.DReduce.HDFS <- function( storage, x, REDUCE, parallel, ... ){

    ## SET a user specific command environment variable here
    cmdenv_arg <- NULL
    cmdenv_arg <- c( cmdenv_arg,
                     sprintf("_HIVE_REDUCER_TO_APPLY_=%s",
                             DSL_serialize_object(REDUCE)) )

    rev <- .MapReduce( x, .HDFS_identity_mapper(), .HDFS_reducer(), cmdenv_arg = cmdenv_arg, ... )

    ## now handle part-xxxx stuff
    ## read all chunk signatures
    .fix_chunknames( DL_storage(x), rev )

    rev
}

## Hadoop Streaming mapper
.HDFS_mapper <- function() {
    function(){
        ## we need this in order to let only the actual output be
        ## written to stdout, does not work with current stemDocument
        ## so not using until needed
        #sink(stderr(), type = "output")
        hive:::redirect_java_output( NULL )

        serialized <- Sys.getenv( "_HIVE_FUNCTION_TO_APPLY_" )
        MAP <- if( serialized == "" ){
            identity
        } else {
            unserialize( charToRaw(gsub("\\n", "\n", serialized, fixed = TRUE)) )
        }

        con <- file("stdin", open = "r")

        chunk <- NA
        while (length(line <- readLines(con, n = 1L, warn = FALSE)) > 0) {
            input <- getFunction("DSL_split_line", where = getNamespace("DSL"))( line )
            if( length(grep("^<<EOF-", input$key)) ){
                chunk <- as.character(input$value["Chunk"])
                break
            }

            result <- MAP( input )

            ## FIXME: should be an object oriented approach here (associative array vs dictionary)
            if( length(result) > 2 )
                for( i in seq_along(result) )
                    writeLines( sprintf("%s\t%s", result[[i]]$key, getFunction("DSL_serialize_object", where = getNamespace("DSL"))(result[[i]]$value)) )
            else
                writeLines( sprintf("%s\t%s", result$key, getFunction("DSL_serialize_object", where = getNamespace("DSL"))(result$value)) )
        }

        ## In the last step we need to add a stamp to this chunk
        ## <key:randomstring, value_serialized:c(firstdocumentkey,lastdocumentkey)>
        if( !is.na(chunk) )
            writeLines( getFunction(".make_chunk_signature", where = getNamespace("DSL"))( chunk ) )

        close(con)

        invisible( TRUE )
    }
}

## For reduce-only jobs we need sort of an identity mapper since there
## must always be a map task before a reduce task according to the
## MapReduce paradigm.
.HDFS_identity_mapper <- function(){
    function(){

        hive:::redirect_java_output( NULL )

        con <- file("stdin", open = "r")

        while (length(line <- readLines(con, n = 1L, warn = FALSE)) > 0) {
                writeLines( line )
        }

        close(con)

        invisible( TRUE )
    }
}

## Hadoop Streaming reducer
.HDFS_reducer <- function() {
    function(){

        ## INIT

        hive:::redirect_java_output( NULL )

        serialized <- Sys.getenv( "_HIVE_REDUCER_TO_APPLY_" )
        REDUCE <- if( serialized == "" ){
            identity
        } else {
            unserialize( charToRaw(gsub("\\n", "\n", serialized, fixed = TRUE)) )
        }

        ## use efficient collector for integer pairlists
        CONCATENATE <- function( collector = FALSE )
            if( collector ){
                getFunction(".collector2", where = getNamespace("DSL"))
            } else {
                base::c
            }

        chunk <- NA
        INTPAIRLIST <- NULL

        ## initialize hash table holding reduce results
        env <- new.env( hash = TRUE, size = 10240 )

        ## CON
        con <- file("stdin", open = "r")

        ## STREAM

        while( length(line <- readLines(con, n = 1L, warn = FALSE)) > 0 ) {
            input <- getFunction("DSL_split_line", where = getNamespace("DSL"))( line )
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
        if(!is.null(INTPAIRLIST))
            if( INTPAIRLIST ){
                env <- lapply(env, getFunction(".collector2", where = getNamespace("DSL")), NULL)
            }
        keys <- names(env)
        for( i in seq_along(keys) )
            writeLines( sprintf("%s\t%s", keys[i],
                                getFunction("DSL_serialize_object", where = getNamespace("DSL"))(REDUCE(env[[ i ]]))) )
        writeLines( getFunction(".make_chunk_signature", where = getNamespace("DSL"))( chunk ) )
    }
}

.MapReduce <- function( x, MAP, REDUCE = NULL, cmdenv_arg = NULL ) {

    x <- as.DList( x )
    rev <- .make_DSL_revision()

    #nmaps <- as.integer(hive::hive_get_parameter("mapred.tasktracker.map.tasks.maximum"))

    ## very important: we don't want to split existing chunks into
    ## separate parts for MapReduce jobs. Hadoop is doing this
    ## automatically unless we set the "mapred.min.split.size"
    ## properly. We set it to the overall block size
    blocksize <- hive::hive_get_parameter("dfs.block.size")

    streaming_args <- sprintf( "-Dmapred.min.split.size=%s", blocksize )
    if( is.null(REDUCE) )
        streaming_args <- sprintf("-Dmapred.reduce.tasks=0 %s", streaming_args)

    ## MAP/REDUCE are functions e.g., provided by R/packages or any user defined
    ## function. It is supplied to the Rscript via an object file written to
    ## disk and exported as environment variable

    ## start the streaming job
    hive::hive_stream( MAP, REDUCE,
                       input  = file.path(DL_storage(x)$base_directory, .revisions(x)[1]),
                       output = file.path(DL_storage(x)$base_directory, rev),
                       cmdenv_arg = cmdenv_arg,
                       streaming_args = streaming_args)

    ## in case the streaming job failed to create output directory return error
    stopifnot( hive::DFS_dir_exists(file.path(DL_storage(x)$base_directory, rev)) )

    rev
}

.get_chunks_after_MapReduce <- function(storage, rev)
    file.path( rev, grep("part-", DS_list_directory(storage, rev), ## removed DSL:::
                              value =TRUE) )


.fix_chunknames <- function( storage, rev ){
    parts <- .get_chunks_after_MapReduce( storage, rev )
    chunks <- lapply( parts, function(part)
                     .read_chunk_signature( storage, part)$value["Chunk"] )
    names(chunks) <- parts
    ind <- !unlist( lapply(chunks, is.null) )
    ## FIXME: need to provide a DS_<> method here?
    for( part in parts[ind] )
        hive::DFS_rename( from = file.path(DS_base_dir(storage), part),
                           to   = file.path(DS_base_dir(storage), rev, basename(chunks[[part]])) )
    invisible( TRUE )
}


