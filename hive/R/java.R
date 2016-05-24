############################################################
## Java DFS support functions
############################################################

add_java_DFS_support <- function(henv){
    if( !all(hadoop_get_jars(henv) %in% .jclassPath()) ){
        .jaddClassPath( hadoop_get_jars(henv) )
    }
    if( hadoop_version(henv) >= "2.6.0" ){
        hadoop_src <- file.path(hadoop_home(henv), "share", "doc", "hadoop")
        hadoop_lib <- file.path(hadoop_home(henv), "share", "hadoop", "common", "lib")
        hadoop_mapreduce <- file.path(hadoop_home(henv), "share", "hadoop", "mapreduce")
        core_default <- .jnew( "org/apache/hadoop/fs/Path", file.path(hadoop_src, "hadoop-project-dist", "hadoop-common", "core-default.xml") )
        core_site    <- .jnew( "org/apache/hadoop/fs/Path", file.path(hadoop_home(henv), "etc", "hadoop", "core-site.xml") )
        hdfs_default <- .jnew( "org/apache/hadoop/fs/Path", file.path(hadoop_src, "hadoop-project-dist", "hadoop-hdfs", "hdfs-default.xml"))
        hdfs_site    <- .jnew( "org/apache/hadoop/fs/Path", file.path(hadoop_home(henv), "etc", "hadoop", "hdfs-site.xml"))
        mapred_default <- .jnew( "org/apache/hadoop/fs/Path", file.path(hadoop_src, "hadoop-mapreduce-client", "hadoop-mapreduce-client-core", "mapred-default.xml") )
        mapred_site  <- .jnew( "org/apache/hadoop/fs/Path", file.path(hadoop_home(henv), "etc", "hadoop", "mapred-site.xml"))
    } else {
        ## add paths to Hadoop configuration files
        core_default <- .jnew("org/apache/hadoop/fs/Path", file.path(hadoop_home(henv), "src", "core", "core-default.xml"))
        core_site <- .jnew("org/apache/hadoop/fs/Path", file.path(hadoop_home(henv), "conf", "core-site.xml"))
        hdfs_default <- .jnew("org/apache/hadoop/fs/Path", file.path(hadoop_home(henv), "src", "hdfs", "hdfs-default.xml"))
        hdfs_site <- .jnew("org/apache/hadoop/fs/Path", file.path(hadoop_home(henv), "conf", "hdfs-site.xml"))
        mapred_default <- .jnew("org/apache/hadoop/fs/Path", file.path(hadoop_home(henv), "src", "mapred", "mapred-default.xml"))
        mapred_site <- .jnew("org/apache/hadoop/fs/Path", file.path(hadoop_home(henv), "conf", "mapred-site.xml"))
    }
    configuration <- .jnew("org/apache/hadoop/conf/Configuration")
    .jcall(configuration, "V", "addResource", core_default)
    .jcall(configuration, "V", "addResource", core_site)
    .jcall(configuration, "V", "addResource", hdfs_default)
    .jcall(configuration, "V", "addResource", hdfs_site)
    .jcall(configuration, "V", "addResource", mapred_default)
    .jcall(configuration, "V", "addResource", mapred_site)
    

  ## Update class loader
  jcl <- .jclassLoader()
  .jcall(configuration, "V", "setClassLoader", .jcast(jcl, "java/lang/ClassLoader"))

  ## create hdfs handler
  hdfs <- .jcall("org/apache/hadoop/fs/FileSystem", "Lorg/apache/hadoop/fs/FileSystem;", "get", configuration)
  ## ioutils class
  ioutils <- .jnew("org/apache/hadoop/io/IOUtils")
  ## store everything in hadoop environment
  assign("configuration", configuration, henv)
  assign("hdfs", hdfs, henv)
  assign("ioutils", ioutils, henv)

  invisible(TRUE)
}

remove_java_DFS_support <- function(henv){
  suppressWarnings( rm("configuration", envir = henv) )
  suppressWarnings( rm("hdfs", envir = henv) )
  suppressWarnings( rm("ioutils", envir = henv) )
  invisible(TRUE)
}


############################################################
## Extractors
############################################################

## Returns the Hadoop DFS configuration object (Java) or NULL
HDFS <- function(henv = hive()){
  hdfs <- tryCatch(get("hdfs", henv), error = identity)
  if(inherits(hdfs, "error"))
    hdfs <- NULL
  hdfs
}

## Returns path as Hadoop DFS path object (Java)
HDFS_path <- function(x)
  .jnew("org/apache/hadoop/fs/Path", x)

## returns the ioutils (Java) object
IOUTILS <- function(henv = hive()){
  ioutils <- tryCatch(get("ioutils", henv), error = identity)
  if(inherits(ioutils, "error"))
    ioutils <- NULL
  ioutils
}

## java output cannot be sinked with standard R tools. Thus, as RWeka
## functions send some strange messages to stdout and err which is
## suboptimal in an Hadoop session we need to redirect the output.
redirect_java_output <- function(x = NULL){
    ## need this when called in streaming jobs (no jvm started)
    rJava::.jinit()
    if(is.null(x)){
        bos <- rJava::.jnew( "java/io/ByteArrayOutputStream" )
        out <- rJava::.jfield( "java/lang/System", , "out" )
        rJava::.jcall( "java/lang/System", "V", "setOut",
               rJava::.jnew("java/io/PrintStream", rJava::.jcast(bos,"java/io/OutputStream")) )
        err <- rJava::.jfield( "java/lang/System", , "err" )
        rJava::.jcall( "java/lang/System", "V", "setErr",
               rJava::.jnew("java/io/PrintStream", rJava::.jcast(bos,"java/io/OutputStream")) )
        out <- list(out = out, err = err)
    }
    else {
        if(is.list(x)){
            rJava::.jcall("java/lang/System", "V", "setOut", out)
            rJava::.jcall("java/lang/System", "V", "setErr", err)
            ## Display them if wanted.
            message(rJava::.jcall(bos, "Ljava/lang/String;", "toString"))
        }
        out <- NULL
    }
    invisible(out)
}
