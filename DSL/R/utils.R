################################################################################
## Utility functions
################################################################################

## key value pairs
## DKey_value_pair<- function( key, value )
##     structure( list(key = key, value = value), class = c("DKey_value_pair", "list") )

## print.DKey_value_pair <- function( x, ... )
##     writeLines( make_string_from_Dkey_value_pair(x) )

## make_string_from_Dkey_value_pair <- function( x )
##     sprintf( "< %s , %s >", as.character(x$key), make_printable(x$value) )

## make_printable <- function(x)
##     UseMethod("make_printable")
## make_printable.character <- identity
## make_printable.numeric <- identity
## make_printable.DKey_value_pair <- function(x)
##     make_string_from_Dkey_value_pair( x )
## make_printable.default <- function( x )
##     sprintf( "( %s )", class(x)[1] )


## Operations on DList objects (getters)

.get_chunks_from_current_revision <- function(x){
    .get_chunks( x )
}

.get_chunks <- function( x, rev ){
    if( missing(rev) )
        rev <- .revisions( x )[1]
    rev <- as.character( rev )
    get(rev, attr(x, "Chunks"))
}

DSL_get_text_mapping_from_revision <- function( x, rev = .revisions(x)[1] )
  attr( x, "Mapping" )[[ rev ]]

.revisions <- function( x )
    get("Revisions", envir = attr( as.DList(x), "Chunks"))

`.revisions<-` <- function( x, value ){
    revs <- .revisions( x )
    to_delete <- revs[! revs %in% value]
    for(rev in to_delete){
        lapply(.get_chunks(x, rev), function(f) DS_unlink(DL_storage(x), f))
        DS_unlink(DL_storage(x), rev)
    }
    assign("Revisions", value, envir = attr( as.DList(x), "Chunks"))
    x
}

## chunk signature (each chunk contains a signature determining the final line)

.make_chunk_signature <- function(chunk)
    sprintf("%s\t%s",
            .stamp(),
            DSL_serialize_object(c(Chunk = chunk)) )

.read_chunk_signature <- function( storage, chunk ){
    split <- strsplit( storage$fetch_last_line(file.path(DS_base_dir(storage), chunk)), "\t" )
    out <- NULL
    if( length(split) )
        out <- list( key = split[[1]][1], value = DSL_unserialize_object(split[[1]][2]) )
    out
}

## make key for final line in MapReduce operations
.stamp <- function(){
    paste("<<EOF", paste(sample(c(letters, 0:9), 10, replace = TRUE),
                collapse = ""),
          Sys.info()["nodename"], sep = "-")
}

DKeys <- function( x )
    attr(x, "Keys")

## thanks to ceeboo: improved collector version to make reduce step run more efficiently
.collector2 <- function(x = NULL, y, ...)
      .Call("_collector2", x, y, package = "DSL")
