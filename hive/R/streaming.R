## Function to facilitate the usage of Hadoop Streaming

## High-level wrapper for hadoop-streaming (MapReduce)
## TODO: ncpu argument -> probably via -D mapred.jobtracker.maxtasks.per.job=..
## TODO: command env arg in henv ?
## TODO: what to do with mapper_args, reducer_args?
## TODO: what if you want to supply more than 1 mapper/reducer function?
hive_stream <- function( mapper, reducer, input, output, henv = hive(),
                         mapper_args = NULL, reducer_args = NULL, cmdenv_arg=NULL, streaming_args=NULL) {
  ## check directories in DFS
  stopifnot( DFS_dir_exists(input, henv) )

  if( missing(reducer) ){
    reducer <- NULL
  }

  ## TODO: encapsulate in separate function something like "prepare..."
  ## are the mapper and reducer really functions?
  .hadoop_check_function_sanity( mapper )
  ## generate mapper and reducer executables
  mapper_exec <- .generate_executable( mapper, .get_hadoop_executable(type = "mapper") )
  ## check if mapper scripts exists
  stopifnot( file.exists(mapper_exec) )
  ## now the reducer (if available)
  reducer_exec <- NULL
  if( is.null(streaming_args) )
      streaming_args <- "-D mapred.reduce.tasks=0"
  if( !is.null(reducer) ){
    .hadoop_check_function_sanity( reducer )
    ## set additional args to Hadoop Streaming
    streaming_args <- paste( sprintf("-D mapred.reduce.tasks=%d", hive_get_nreducer(henv)), streaming_args, collapse = " " )
    reducer_exec <- .generate_executable( reducer, .get_hadoop_executable(type = "reducer") )
    stopifnot( file.exists(reducer_exec) )
  }
  ## check args
  if( is.null(mapper_args) )
    mapper_args <- ""
  if( is.null(reducer_args) )
    reducer_args <- ""

  ## start hadoop streaming
  msg <- .hadoop_streaming( mapper_exec, reducer_exec, input, output, mapper_args, reducer_args, streaming_args, cmdenv_arg, henv )

  ## delete temporary created mapper/reducer scripts
  .hadoop_cleanup( files = c(mapper_exec, reducer_exec) )
  invisible( msg )
}

.hadoop_streaming <- function( mapper, reducer, input, output, mapper_args, reducer_args,
                               streaming_args, cmdenv_arg, henv ){
  files <- paste( unique(c("", mapper, reducer)), collapse = " -file " )
  mapper_arg <- sprintf( "-mapper '%s %s'", mapper, as.character(mapper_args) )
  reducer_arg <- ""
  #writeLines(sprintf("DEBUG: reducer: %s", as.character(!is.null(reducer))))
  if(!is.null(reducer))
    reducer_arg <- sprintf("-reducer '%s %s'", reducer, as.character(reducer_args))
  if(is.null(cmdenv_arg))
    cmdenv_arg <- ""
  else
    cmdenv_arg <- sprintf( "-cmdenv '%s'", paste(as.character(cmdenv_arg), collapse = " ") )
  out <- system(sprintf('%s jar %s %s -input %s -output %s %s %s %s %s',
                 hadoop(henv), hadoop_streaming(henv), streaming_args, input, output,
                 mapper_arg, reducer_arg, files, cmdenv_arg),
         intern = TRUE)
  out
}

## creates the mapper and/or reducer script
## see TODOs above: should be stored in DFS
.generate_executable <- function(x, script){
  foo <- paste(deparse(functionBody(x)), collapse = "\n")
  ## we use 'Rscript'
  prefix <- "#!/usr/bin/env Rscript\n"
  ## write file
  cat(sprintf("%s%s", prefix, foo), file = script, sep = "\n")
  ## make file executable
  status <- system(sprintf("chmod 775 %s", script), ignore.stderr = TRUE)
  if(status){
    warning("no executable found!")
    invisible(FALSE)
  }
  invisible(script)
}

.hadoop_check_function_sanity <- function(x){
  stopifnot(is.function(x))
  ## TODO: how many arguments? typically, x has no arguments
}

.hadoop_cleanup <- function(files){
  if(all(file.exists(files)))
    unlink(files)
}

## temp file
## TODO: store that one in DFS
.get_hadoop_executable <- function(type = "mapper"){
  sprintf("_hadoop_%s_", type)
}

