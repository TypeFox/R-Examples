## Functions related to the Hadoop Distributed File System (HDFS)
## Author: Stefan Theussl

## out of simplicity queries status of / in DFS (Java)
DFS_is_available <- function( henv = hive() ) {
    ## first check if DFS is registered
    if( DFS_is_registered(henv) ){
        ## then check if we can access DFS
        stat <- .DFS_stat( "/", henv )
        if( !(is.null(stat) || is.na(stat)) )
            return( TRUE )
    }
    FALSE
}

DFS_is_registered <- function(henv = hive()){
    if( is.null(HDFS(henv)) || is.null(IOUTILS(henv)) ){
        return(FALSE)
    }
    TRUE
}

## does file exist in DFS? (Java)
DFS_file_exists <- function( file, henv = hive() ) {
    stopifnot(DFS_is_registered(henv))
    hdfs <- HDFS(henv)
    hdfs$exists(HDFS_path(file))
}

## does dir exist in DFS? (Java)
DFS_dir_exists <- function( path, henv = hive() ) {
  path <- file.path(path)
  status <- tryCatch(.DFS_getFileStatus(path, henv ), error = identity)
  if(inherits(status, "error"))
    return(FALSE)
  status$isDir()
}

## create dir in DFS (Java)
## TODO: throws getClass error although function calls work in global env
## interestingly, when debugging and calling DFS_mkdir() twice it works ...
DFS_dir_create <- function( path, henv = hive() ) {
  if( DFS_dir_exists(path, henv) ) {
    warning( sprintf("directory '%s' already exists.", path) )
    return( invisible(FALSE) )
  }
  if( DFS_file_exists(path, henv) ) {
    warning( sprintf("'%s' already exists but is not a directory", path) )
    return( invisible(FALSE) )
  }
  status <- .DFS_mkdir( path, henv )
  if( is.null(status) ) {
    warning( sprintf("cannot create dir '%s'.", path) )
    return( invisible(FALSE) )
  }
  invisible( TRUE )
}

## Delete files in DFS (Java)
DFS_delete <- function( file, recursive = FALSE, henv = hive() ) {
  if( DFS_dir_exists(file, henv) && !recursive){
    warning(sprintf("cannot remove directory '%s'. Use 'recursive = TRUE' instead.", file))
    return(FALSE)
  }

  status <- .DFS_delete( file, henv )
  if(!status){
    warning(sprintf("cannot remove file '%s'.", file))
    return(FALSE)
  }
  TRUE
}

DFS_dir_remove <- function(path, recursive = TRUE, henv = hive()){
  if( DFS_dir_exists(path, henv) ){
    DFS_delete(path, recursive, henv)
    TRUE
  } else {
    warning(sprintf("'%s' is not a directory.", path))
    FALSE
  }
}


## private int ls(String srcf, boolean recursive) throws IOException {
##    Path srcPath = new Path(srcf);
##    FileSystem srcFs = srcPath.getFileSystem(this.getConf());
##    FileStatus[] srcs = srcFs.globStatus(srcPath);
##    if (srcs==null || srcs.length==0) {
##      throw new FileNotFoundException("Cannot access " + srcf +
##          ": No such file or directory.");
##    }
##
##    boolean printHeader = (srcs.length == 1) ? true: false;
##    int numOfErrors = 0;
##    for(int i=0; i<srcs.length; i++) {
##      numOfErrors += ls(srcs[i], srcFs, recursive, printHeader);
##    }
##    return numOfErrors == 0 ? 0 : -1;
##  }


DFS_list <- function( path = ".", henv = hive() ) {
  globstat <- .DFS_stat(path, henv)
  if( is.null(globstat) ){
    warning(sprintf("'%s' is not a readable directory", path))
    return(character(0))
  }

  splitted <- strsplit(grep(path, .DFS_intern("-ls", path, henv), value = TRUE), path)
  sapply(splitted, function(x) basename(x[2]))
}

DFS_cat <- function( file, con = stdout(), henv = hive() ){
  stopifnot( DFS_file_exists(file, henv) )
  .DFS("-cat", file, henv)
}

DFS_rename <- function( from, to, henv = hive() ){
    stopifnot( DFS_file_exists(from, henv) )
    .DFS_rename( from, to, henv )
}


DFS_tail <- function(file, n = 6L, size = 1024, henv = hive() ){
  stopifnot( as.integer(n) > 0L )
  stopifnot( DFS_file_exists(file, henv) )
  out <- .DFS_tail(file, size, henv = henv)
  len <- length(out)
  out[(len - (n - 1)) : len]
}

.DFS_tail <- function(file, size = 1024, henv = hive()){
    stopifnot(DFS_is_registered(henv))
    hdfs <- HDFS(henv)
    ioutils <- IOUTILS(henv)

    hdfs_file <- HDFS_path(file)
    len <- hdfs$getFileStatus(hdfs_file)$getLen()
    offset <- ifelse(len > size, len - size, 0)

    inputstream <- hdfs$open(hdfs_file)
    inputstream$seek(.jlong(offset))

    ## we need to copy the contents of the file to an output stream
    ## Thus, for the time being we use the JRI class to RConsoleOutputStream divert
    ## the outputstream to the R console
    routput <- .jnew("org/rosuda/JRI/RConsoleOutputStream", .jengine(TRUE), as.integer(0))
    ## now we need to capture the contents from the console usingg a text connection
    ## we save the results in the object out
    out <- character(0)
    con <- textConnection("out", open = "w", local = TRUE)
    sink(file = con)
    ioutils$copyBytes(inputstream, routput, as.integer(1024), TRUE)
    sink()
    close(con)
    out
}

# Load local files into hadoop and distribute them along its nodes
DFS_put <- function( files, path = ".", henv = hive() ) {
  if(length(files) == 1)
    status <- .DFS("-put", paste(files, path), henv )
  else {
    if( !DFS_dir_exists(path, henv) )
      DFS_dir_create( path, henv )
    status <- .DFS("-put", paste(paste(files, collapse = " "), path), henv )
  }
  if( status ){
    warning( sprintf("Cannot put file(s) to '%s'.", path) )
    return( invisible(FALSE) )
  }
  invisible( TRUE )
}

## serialize R object to DFS
DFS_put_object <- function( obj, file, henv = hive() ) {
  con <- .DFS_pipe( "-put", file, open = "w", henv = henv )
  status <- tryCatch(serialize( obj, con ), error = identity)
  close.connection(con)
  if(inherits(status, "error"))
    stop("Serialization failed.")
  invisible(file)
}

## worse performance than read_lines2, reason: paste
DFS_write_lines <- function( text, file, henv = hive() ) {
    stopifnot(DFS_is_registered(henv = henv))
    if(DFS_file_exists(file)){
        warning(sprintf("file '%s' already exists.", file))
        return(NA)
    }

    if(!length(text))
        stop("text length of zero not supported.")

    hdfs <- HDFS(henv)

    outputstream <- hdfs$create(HDFS_path(file))
    for( i in seq_along(text) ){
        outputstream$writeBytes(text[i])
        outputstream$writeBytes("\n")
    }
    outputstream$close()

    invisible(file)
}

## TODO: there is a line reader class provided by Hadoop, see also
## http://hadoop.apache.org/common/docs/r0.20.1/api/org/apache/hadoop/util/LineReader.html
##       cannot instantiate using LineReader <- .jnew("org/apache/hadoop/util/LineReader")
DFS_read_lines <- function( file, n = -1L, henv = hive() ) {
    if(!DFS_file_exists(file)){
        warning(sprintf("file '%s' does not exists.", file))
        return(NA)
    }
    hdfs <- HDFS(henv)
    ioutils <- IOUTILS(henv)
    offset <- 0
    inputstream <- hdfs$open(HDFS_path(file))
    if( n <= 0 ){
        inputstream$seek(.jlong(offset))

        ## we need to copy the contents of the file to an output stream
        ## Thus, for the time being we use the JRI class to RConsoleOutputStream divert
        ## the outputstream to the R console
        routput <- .jnew("org/rosuda/JRI/RConsoleOutputStream", .jengine(TRUE), as.integer(0))
        ## now we need to capture the contents from the console usingg a text connection
        ## we save the results in the object out
        con <- textConnection("out", open = "w", local = TRUE)
        sink(file = con)
        ioutils$copyBytes(inputstream, routput, as.integer(1024), TRUE)
        sink()
        close(con)
        #inputstream$close()
    }
    else {
        out <- character(n)
        for(i in 1:n)
            out[i] <- inputstream$readLine()
        inputstream$close()
    }
    out
}

## serialize R object from DFS
DFS_get_object <- function( file, henv = hive() ) {
  con <- .DFS_pipe( "-cat", file, open = "r", henv = henv )
  obj <- tryCatch( unserialize(con), error = identity)
  close.connection(con)
  if(inherits(obj, "error"))
     return(NA)
  obj
}

## deletes a file or empty directory
## returns TRUE if successful and FALSE otherwise
## caution: always deletes recursively!
.DFS_delete <- function(x, henv){
    stopifnot( DFS_is_registered(henv) )
    hdfs <- HDFS(henv)
    hdfs$delete(HDFS_path(x))
}

## creates directory on DFS
## returns TRUE if successful and NULL otherwise
.DFS_mkdir <- function(x, henv){
    stopifnot( DFS_is_registered(henv) )
    hdfs <- HDFS(henv)
    hdfs$mkdirs(HDFS_path(x))
}

.DFS_stat <- function(x, henv){
    stopifnot( DFS_is_registered(henv) )
    hdfs <- HDFS(henv)
    stat <- hdfs$globStatus(HDFS_path(x))
    if(is.null(stat)){
        warning(sprintf("cannot stat '%s': No such file or directory", x))
        return(NULL)
    }
    ## for the time being return TRUE
    ## TODO: this should return an R object containing the stat information
    TRUE
}

.DFS_rename <- function( from, to, henv ){
    stopifnot( DFS_is_registered(henv) )
    hdfs <- HDFS(henv)
    hdfs$rename(HDFS_path(from), HDFS_path(to))
}

.DFS_getFileStatus <- function(x, henv){
    stopifnot( DFS_is_registered(henv) )
    hdfs <- HDFS(henv)
    hdfs$getFileStatus(HDFS_path(x))
}

################################################################################
## Old command line wrappers
## FIXME: all of them have to be replaced by Java routines

.DFS <- function( cmd, args, henv )
  system( .DFS_create_command(cmd, args, henv), ignore.stderr = TRUE )

.DFS_pipe <- function( cmd, args, open = "w", henv ){
  if(open == "w")
    pipe(.DFS_create_command(cmd, sprintf("- %s", args), henv), open = open)
  else
    pipe(.DFS_create_command(cmd, args, henv), open = open)
}

.DFS_intern <- function( cmd, args, henv )
  system( .DFS_create_command(cmd, args, henv), intern = TRUE, ignore.stderr = TRUE )

.DFS_create_command <- function( cmd, args, henv )
  sprintf("%s fs %s %s", hadoop(henv), cmd, args)

################################################################################
## use with caution
## FIXME: not working yet, too dangerous
.DFS_format <- function(henv){
  ##machines, DFS_root= "/var/tmp/hadoop"
  stopifnot(hive_stop(henv))
  machines <- unique(c(hive_get_slaves(henv), hive_get_masters(henv)))
  DFS_root <- gsub("\\$\\{user.name\\}", system("whoami", intern = TRUE),
                   hive_get_parameter("hadoop.tmp.dir", henv))
  for(machine in machines){
    ## delete possibly corrupted file system
    command <- sprintf("ssh %s 'rm -rf %s/*
rm -rf %s-*' ", machine, DFS_root, DFS_root)
    system(command)
  }
  # reformat DFS
  system(sprintf("%s namenode -format", hadoop(henv)))
}
