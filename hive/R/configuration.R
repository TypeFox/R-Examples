## Getters
hive_get_parameter <- function(x, henv = hive()){
    if( hadoop_version(henv) < "0.20.0" ) {
        ## first search in hadoop-site configuration (overrules defaults)
        out <- .hadoop_configuration("hadoop_site", henv)[x]
        ## if not found then return value from default configuration
        if(is.na(out))
            out <- .hadoop_configuration("hadoop_default", henv)[x]
    } else {
        site <- list()
        for( i in grep("_site", names(get("config_files", henv))) )
            site <- unlist( c( site, get("config_files", henv)[[ i ]]) )
        default <- list()
        for( i in grep("_default", names(get("config_files", henv))) )
            default <- unlist( c( default, get("config_files", henv)[[ i ]]) )
        out <- site[x]
        ## if not found then return value from default configuration
        if(is.na(out))
            out <- default[x]
    }
    out
}

hive_get_slaves <- function(henv = hive()){
 .hadoop_configuration("slaves", henv)
}

hive_get_masters <- function(henv = hive()){
  .hadoop_configuration("masters", henv)
}

.hadoop_configuration <- function(x, henv){
  get("config_files", henv)[[x]]
}

## Setters

## FIXME: not updated yet
hive_set_slaves <- function(slaves, henv){
  hive_stop(henv)
  slave_conf <- file.path(hadoop_home(henv), "conf", "slaves")
  writeLines(slaves, con = slave_conf)
  hive_start(henv)
}

## Hadoop config XML parser

## FIXME: con argument
get_hadoop_config <- function(x, dir){
    if( !file.exists(file.path(dir, x)) )
        out <- NA
    else {
        infile <- xmlRoot(xmlTreeParse(file.path(dir, x)))
        ind_prop <- unlist(xmlApply(infile, xmlName) == "property")
        prop_names <- hadoop_parse_xml(infile[ind_prop], "name")
        out <- hadoop_parse_xml(infile[ind_prop], "value")
        names(out) <- prop_names
    }
    out
}

## returns the right function for coercion
## FIXME: 'value' can be of different type (integer, character, logical, NA)
hadoop_xml_return_type <- function(x){
  switch(x,
         "name"  = return(as.character),
         "value" = return(as.character),
         "description" = return(as.character)
         )
  stop("'x' can only be 'name', 'value' or 'description'")
}

## parse values from xml tree
hadoop_parse_xml <- function(x, what){
    as_type <- hadoop_xml_return_type(what)
    as_type( unlist(lapply(x, function(x) {out <- xmlValue(x[[what]])
                                           if(!length(out))
                                               out<-NA
                                           out})) )
}
