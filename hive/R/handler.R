## Handling the Hadoop framework

## .hinit() initializes the Hadoop framework and returns the corresponding environment.
.hinit <- function( hadoop_home ) {
  ## use installation in HADOOP_HOME if it exists, otherwise the Debian default
  if( missing(hadoop_home) )
    hadoop_home <- ifelse( utils::file_test("-d", Sys.getenv("HADOOP_HOME")), Sys.getenv("HADOOP_HOME"), "/etc/hadoop" )
  tmp <- tryCatch( hive_create(hadoop_home), error = identity )
  hive <- if( inherits(tmp, "error") )
    .hive_default_env()
  else
    tmp
  hive
}

## See also .create_hive_from_installation.
## Returns on object of class 'hive'.
hive_create <- function( hadoop_home ){
  if( missing(hadoop_home) ){
    hadoop_home_exists <- file.exists( .hadoop_home_dirs() )
    ifelse( any(hadoop_home_exists),
            hadoop_home <- .hadoop_home_dirs()[hadoop_home_exists],
            stop( "could not find Hadoop home directory!")
           )
  }
  hive <- .create_hive_from_installation( file_path_as_absolute(hadoop_home) )
  hive_set_nreducer( hive_default_nreducer(hive), hive )
  class( hive ) <- "hive"
  hive
}

.hadoop_home_dirs <- function(){
  c("/etc/hadoop")
}
## Given a pointer to a Hadoop installation directory, this function
## creates an environment containing all information about the Hadoop
## cluster.  We store the hadoop home directory, the hadoop version,
## and the parsed configuration files in a separate R environment.
.create_hive_from_installation <- function( hadoop_home ){
  if( !file.exists(hadoop_home) )
    stop( sprintf("There is no directory '%s'.", hadoop_home) )
  hive <- new.env()
  hvers <- hadoop_get_version( hadoop_home )
  streaming_home <- c( file.path(hadoop_home, "share/hadoop/tools/lib"),
                       file.path(hadoop_home, "contrib/streaming"),
                       file.path("/usr/lib/hadoop/contrib/streaming") )

  jar <- lapply(lapply(streaming_home, dir), function(x) grep("(hadoop.+streaming.+jar$)", x, value = TRUE))
  ind <- unlist(lapply(jar, function(x) length(x) > 0))
  hadoop_streaming <- file.path(streaming_home[ind][1], unlist(jar[ind])[1])
  ## OLD CONFIG: config files have been split up and located in different places since version 0.20.0
  if( hvers < "0.20.0" ){
      local( {
          hadoop <- file.path(hadoop_home, "bin", "hadoop")
          version <- hvers
          stopifnot(file.exists(hadoop))
          config_files <- list(hadoop_default = get_hadoop_config("hadoop-default.xml", file.path(hadoop_home,"conf")),
                               hadoop_site = get_hadoop_config("hadoop-site.xml", file.path(hadoop_home,"conf")),
                                slaves = readLines(file.path(hadoop_home, "conf", "slaves")),
                                masters = readLines(file.path(hadoop_home, "conf", "masters")))
      }, hive )
      
  } else {
  ## CURRENT CONFIGS
      if( hvers >= "2.6.0" ) {
          hadoop_src       <- file.path( hadoop_home, "share", "doc", "hadoop" )
          hadoop_common    <- file.path( hadoop_home, "share", "hadoop", "common" )
          hadoop_lib       <- file.path( hadoop_common, "lib")
          hadoop_mapreduce <- file.path( hadoop_home, "share", "hadoop", "mapreduce" )
          hadoop_hdfs      <- file.path( hadoop_home, "share", "hadoop","hdfs" )
          hadoop_tools     <- file.path( hadoop_home, "share", "hadoop","tools", "lib" )
          local( {
              hadoop <- if( utils::file_test("-x", file.path(hadoop_home, "bin", "hadoop")) )
                            file.path(hadoop_home, "bin", "hadoop")
                        else
                            Sys.which("hadoop")
              version <- hvers
              stopifnot(file.exists(hadoop))
              
              config_files <- list(core_default = get_hadoop_config("core-default.xml", file.path(hadoop_src, "hadoop-project-dist", "hadoop-common")),                          
                                   core_site = get_hadoop_config("core-site.xml", file.path(hadoop_home, "etc", "hadoop")),
                                   hdfs_default = get_hadoop_config("hdfs-default.xml", file.path(hadoop_src, "hadoop-project-dist", "hadoop-hdfs")),
                                   hdfs_site = get_hadoop_config("hdfs-site.xml", file.path(hadoop_home, "etc", "hadoop")),
                                   mapred_default = get_hadoop_config("mapred-default.xml", file.path(hadoop_src, "hadoop-mapreduce-client", "hadoop-mapreduce-client-core")),
                                   mapred_site = get_hadoop_config("mapred-site.xml", file.path(hadoop_home, "etc", "hadoop")),
                                   slaves = suppressWarnings(tryCatch(readLines(file.path(hadoop_home, "etc", "hadoop", "slaves")), error = function(x) NA)),
                                   masters = suppressWarnings(tryCatch(readLines(file.path(hadoop_home, "etc", "hadoop", "masters")), error = function(x) NA)))
              
              commonsloggingjar       <- grep( "commons-logging-[0-9].*[.]jar", dir(hadoop_lib), value = TRUE )
              commonsconfigurationjar <- grep( "commons-configuration-[0-9].*[.]jar", dir(file.path(hadoop_lib)), value = TRUE )
              commonslangjar <- grep( "commons-lang-[0-9].*[.]jar", dir(hadoop_lib), value = TRUE )
              commonscli <- grep( "commons-cli-[0-9].*[.]jar", dir(hadoop_lib), value = TRUE )
              commonsiojar <- grep( "commons-io-[0-9].*[.]jar", dir(hadoop_lib), value = TRUE )
              guavajar <- grep( "guava-[0-9].*[.]jar", dir(hadoop_tools), value = TRUE )
              commonscollectionsjar <- grep( "commons-collections-[0-9].*[.]jar", dir(hadoop_tools), value = TRUE )
              hadoop_auth <- sprintf( "hadoop-auth-%s.jar", hvers )
              slf4jjar <- grep( "slf4j-api-[0-9].*[.]jar", dir(hadoop_lib), value = TRUE )
              slf4jlogjar <- grep( "slf4j-log4j12-[0-9].*[.]jar", dir(hadoop_lib), value = TRUE )
              servletjar <- grep( "servlet-api-[0-9].*[.]jar", dir(hadoop_lib), value = TRUE )
              log4jjar <- grep( "log4j-[0-9].*[.]jar", dir(hadoop_tools), value = TRUE )
              protobufjavajar <- grep( "protobuf-java-[0-9].*[.]jar", dir(hadoop_lib), value = TRUE )
              htracejar <- grep("htrace-core-[0-9].*[.]jar", dir(hadoop_lib), value = TRUE )
              hadoop_jars <- c( file.path(hadoop_common, c(sprintf("hadoop-common-%s.jar", hvers))),
                                file.path(hadoop_hdfs, c(sprintf("hadoop-hdfs-%s.jar", hvers))),
                                file.path(hadoop_mapreduce, c(sprintf("hadoop-mapreduce-client-common-%s.jar", hvers))),
                                file.path(hadoop_lib, c(commonscli, commonsloggingjar, commonsconfigurationjar, commonslangjar, commonsiojar)),
                                file.path(hadoop_tools, c(guavajar, commonscollectionsjar, hadoop_auth, servletjar)),
                                file.path(hadoop_lib, c(slf4jjar, slf4jlogjar, log4jjar, protobufjavajar, htracejar)) )
              stopifnot( all(file.exists(hadoop_jars)) )
              
      }, hive )
  } else {
      ## VERSION 1.0 CONFIGS
      if( file.exists( file.path(hadoop_home, "src/core" )) ){
          hadoop_src <- file.path(hadoop_home, "src" )
      } else if( file.exists( file.path("/usr/src/hadoop-0.20", "core" )) ){
          hadoop_src <- "/usr/src/hadoop-0.20"
      } else {
          hadoop_src <- system.file("defaults", package = "hive")
      }
      local( {
          hadoop <- if( utils::file_test("-x", file.path(hadoop_home, "bin", "hadoop")) )
              file.path(hadoop_home, "bin", "hadoop")
          else
              Sys.which("hadoop")
          version <- hvers
          stopifnot(file.exists(hadoop))

          config_files <- list(core_default = get_hadoop_config("core-default.xml", file.path(hadoop_src, "core")),
                                core_site = get_hadoop_config("core-site.xml", file.path(hadoop_home, "conf")),
                                hdfs_default = get_hadoop_config("hdfs-default.xml", file.path(hadoop_src, "hdfs")),
                                hdfs_site = get_hadoop_config("hdfs-site.xml", file.path(hadoop_home, "conf")),
                                mapred_default = get_hadoop_config("mapred-default.xml", file.path(hadoop_src, "mapred")),
                                mapred_site = get_hadoop_config("mapred-site.xml", file.path(hadoop_home, "conf")),
                                slaves = suppressWarnings(tryCatch(readLines(file.path(hadoop_home, "conf", "slaves")), error = function(x) NA)),
                                masters = suppressWarnings(tryCatch(readLines(file.path(hadoop_home, "conf", "masters")), error = function(x) NA)))

      
          ## FIXME: Debian packages not available anymore, thus using CDH paths
          ## hadoop_jars <- file.path("/usr/share/java/", c("hadoop-core.jar", "commons-logging.jar"))
          commonsloggingjar <- grep("commons-logging-[0-9].*[.]jar", dir(file.path("/usr/lib/hadoop", "lib")), value = TRUE)
          hadoop_jars <- file.path("/usr/lib/hadoop", c("hadoop-core.jar", file.path("lib",commonsloggingjar)))

          if( !all(file.exists(hadoop_jars)) ){
              commonsloggingjar <- grep("commons-logging-[0-9].*[.]jar", dir(file.path(hadoop_home, "lib")), value = TRUE)
              commonsconfigurationjar <- grep("commons-configuration-[0-9].*[.]jar", dir(file.path(hadoop_home, "lib")), value = TRUE)
              commonslangjar <- grep("commons-lang-[0-9].*[.]jar", dir(file.path(hadoop_home, "lib")), value = TRUE)
              hadoop_jars <- if(version >= "0.20.203")
                  file.path(hadoop_home, c(sprintf("hadoop-core-%s.jar", version), file.path("lib", commonsconfigurationjar), file.path("lib", commonslangjar), file.path("lib", commonsloggingjar)))
              else
                  file.path(hadoop_home, c(sprintf("hadoop-%s-core.jar", version),  file.path("lib", commonsloggingjar)))
        }
      }, hive )
  }}
  hive
}

## Default environment: NA
.hive_default_env <- function(){
  NA
}

## Checks if henv inherits from class 'hive'
hive_is_valid <- function( henv ){
  inherits(henv, "hive")
}

## Provides information about the "hive"
print.hive <- function( x, ... ){
  writeLines( "HIVE: Hadoop Cluster" )
  writeLines( sprintf("- Avail. datanodes: %d", length(hive_get_slaves(x))) )
  writeLines( sprintf("'- Max. number Map tasks per datanode: %s",
                      ifelse( hadoop_version(x) >= "2.6.0",
                             hive_get_parameter("mapreduce.tasktracker.map.tasks.maximum", x),
                             hive_get_parameter("mapred.tasktracker.map.tasks.maximum", x))) )
  writeLines( sprintf("'- Configured Reducer tasks: %d",
                      hive_get_nreducer(x)) )
}

summary.hive <- function( object, ... ){
    print(object)
    writeLines( "---" )
    writeLines( sprintf("- Hadoop version: %s", hadoop_version(hive(object))) )
    writeLines( sprintf("- Hadoop home/conf directory: %s", hadoop_home(object)) )
    writeLines( sprintf("- Namenode: %s", hive_get_masters(object)) )
    writeLines( "- Datanodes:")
    writeLines( sprintf("'- %s\n", hive_get_slaves(object)) )
}

## Start and stop a Hadoop cluster.
## NOTE: Java DFS support is only available for the current cluster.
##       Thus, add/remove DFS support in each call to hive_start/stop
hive_start <- function( henv = hive() ){
    ## does nothing if hadoop daemons are already running
    if( DFS_is_registered(henv) )
        return( invisible(TRUE) )
    ## otherwise start the Hadoop framework ...
    ## FIXME: on Debian systems Hadoop daemons are controlled via
    ##        /etc/init.d scripts (started automatically)
    start_dfs_sh <- ifelse( hadoop_version(henv) >= "2.6.0",
                            file.path(hadoop_home(henv), "sbin", "start-dfs.sh"),
                            file.path(hadoop_home(henv), "bin", sprintf("%s-all.sh", "start")) )
    if( file.exists(start_dfs_sh) )
        hadoop_framework_control( "start", henv )
    ## ... and add Java DFS support
    ## NOTE: really important that Java DFS support is added AFTER
    ## framework has been started
    status <- add_java_DFS_support( henv = hive() )
    ## if there are problems starting hive, close it
    if( !status ){
        # writeLines(msg)
        suppressWarnings( hive_stop(henv) )
        return( invisible(FALSE) )
    }
    invisible( TRUE )
}

hive_stop <- function( henv = hive() ){
    if( DFS_is_registered(henv) ){
        remove_java_DFS_support( henv )
        ## FIXME: on Debian systems Hadoop daemons are controlled via
        ##        /etc/init.d scripts (stopped automatically)
        stop_dfs_sh <- ifelse( hadoop_version(henv) >= "2.6.0",
                              file.path(hadoop_home(henv), "sbin", "stop-dfs.sh"),
                              file.path(hadoop_home(henv), "bin", sprintf("%s-all.sh", "stop")) )
        if( file.exists(stop_dfs_sh) )
            hadoop_framework_control( "stop", henv )
    }
    else
        warning( "No Hadoop cluster running. Nothing to stop." )
    invisible( TRUE )
}

## FIXME: Very simple query of hadoop status. Probably better to use pid?
hive_is_available <- function( henv = hive() ){
    stopifnot( hive_is_valid(henv) )
    suppressWarnings( DFS_is_available(henv) )
}


hive_get_nreducer <- function( henv = hive() ){
  get( "nreducer", henv )
}

hive_set_nreducer <- function( n, henv = hive() ) {
  assign( "nreducer", as.integer(n), henv )
}

hive_default_nreducer <- function( henv = hive() ){
    as.integer(round( length(hive_get_slaves(henv)) / 1.5 ))
}

## Internal extractor functions
hadoop <- function( henv )
  get( "hadoop", henv )

hadoop_home <- function( henv )
  get( "hadoop_home", henv )

hadoop_streaming <- function( henv )
  get( "hadoop_streaming", henv )

hadoop_version <- function( henv )
  get( "version", henv )

hadoop_get_jars <- function( henv )
  get( "hadoop_jars", henv )

hadoop_get_version <- function( hadoop_home ){
  version_file <- file.path(hadoop_home, "contrib", "hod", "bin", "VERSION")
  if( file.exists( version_file ) )
    readLines( version_file )
  else
    gsub("Hadoop ", "", system("hadoop version", intern = TRUE, ignore.stderr = TRUE)[1])
}

## Controlling the Hadoop framework: currently using the start/stop_all.sh scripts
## in the $HADOOP_HOME/bin directory
## FIXME: This may not be platform independent
hadoop_framework_control <- function( action = c("start", "stop"), henv ){
    action <- match.arg(action)
    dfs_action_sh <- ifelse( hadoop_version(henv) >= "2.6.0",
                             file.path(hadoop_home(henv), "sbin", sprintf("%s-dfs.sh", action)),
                             file.path(hadoop_home(henv), "bin", sprintf("%s-all.sh", action)) )
    system( dfs_action_sh, intern = TRUE, ignore.stderr = TRUE )
}
