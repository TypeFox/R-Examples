## Authors: Ingo Feinerer, Stefan Theussl

## "DCorpus" class
.DCorpus <- function( x, keep, cmeta, dmeta ) {
  ## use revisions? Default: TRUE. This can be turned off using
  ## setRevision() replacement function.
  if( missing(keep) )
      keep <- TRUE
  structure(list(content = x, meta = cmeta, dmeta = dmeta, keep = keep),
            class = c("DCorpus", "Corpus"))
}

DistributedCorpus <-
DCorpus <-
function(x,
         readerControl = list(reader = reader(x), language = "en"),
         storage = NULL, keep = TRUE, ...)
{
    ## For the moment we
    ##   - only support a directory as source (DirSource)
    ## TODO: add DList source
    ## FIXME: in earlier versions DCorpus had a keys argument for supplying user chosen keys
    if( !inherits(x, "DirSource") )
            stop("unsupported source type (use DirSource instead)")

    if (is.null(readerControl$reader))
        readerControl$reader <- reader(x)
    if (inherits(readerControl$reader, "FunctionGenerator"))
        readerControl$reader <- readerControl$reader(...)
    if (is.null(readerControl$language))
        readerControl$language <- "en"

    # Check for parallel element access
    if (is.function(getS3method("pGetElem", class(x), TRUE))) {
        elem <- pGetElem(x)
        # NOTE: DirSource guarantees !is.null(x$uri)
        names(elem) <- unlist(lapply(elem, function(x) basename(x$uri)))
        tdl <- DMap(as.DList(elem, DStorage = storage), function(keypair) list(key = keypair$key, value = readerControl$reader(keypair$value, readerControl$language, keypair$key)) )
    }
    else
        stop("Non-vectorized operation not yet implemented.")

    df <- data.frame(row.names = seq_along(tdl))
    cm <- structure(list(), class = "CorpusMeta")
    .DCorpus( tdl, keep, cm, df )
}

`[.DCorpus` <- getS3method("[", "VCorpus")

`[[.DCorpus` <- getS3method("[[", "VCorpus")

`[[<-.DCorpus` <- getS3method("[[<-", "VCorpus")

as.list.DCorpus <-
function(x, ...)
    as.list(x$content)

as.DistributedCorpus <-
as.DCorpus <- function(x, storage = NULL, ...){
  UseMethod("as.DCorpus")
}

as.DCorpus.DCorpus <- function(x, storage = NULL, ...)
    x

as.DCorpus.VCorpus <- function(x, storage = NULL, ...){
    dl <- as.DList(as.list(x), DStorage = storage, ...)
    .DCorpus( dl,
              keep = TRUE,
              cmeta = meta(x, type = "corpus"),
              dmeta = meta(x, type = "indexed") )
}

as.VCorpus.DCorpus <- function(x)
    structure(list(content = as.list(x),
                   meta = meta(x, type = "corpus"),
                   dmeta = meta(x, type = "indexed")),
              class = c("VCorpus", "Corpus"))

format.DCorpus <- getS3method("format", "VCorpus")

length.DCorpus <- getS3method("length", "VCorpus")

meta.DCorpus <- getS3method("meta", "VCorpus")

names.DCorpus <- getS3method("names", "VCorpus")

print.DCorpus <- getS3method("print", "VCorpus")

summary.DCorpus <- function( object, ... ) {
    print(object)
    cat( "\nDCorpus revisions:\n" )
    cat( strwrap(paste(unlist(getRevisions(object)), collapse = " "), indent = 2, exdent = 2), "\n" )
    cat( sprintf("DCorpus active revision: %s\n\n", DSL:::.revisions(object$content)[1]) )
    print( DL_storage(object$content) )
}

## Get all available revisions from the DC
getRevisions <- function( corpus ){
    DSL:::.revisions( corpus$content )
}

## Set active revision in the DC to the specified revision
setRevision <- function( corpus, revision ){
    pos <- as.character(revision) == getRevisions(corpus)
    if( !any(pos) )
        warning( "invalid revision" )
    DSL:::.revisions( corpus$content ) <- c( revision, getRevisions(corpus)[!pos] )
    invisible(corpus)
}

## the setRevision replacement function is used to turn revisions on and off
keepRevisions <- function( corpus )
    corpus$keep

`keepRevisions<-` <- function( corpus, value ){
    stopifnot( length(value) == 1L )
    stopifnot( is.logical(value) )
    corpus$keep <- value
    corpus
}

## remove a given revision
removeRevision <- function( corpus, revision ){
    pos <- revision == getRevisions(corpus)
    if( !any(pos) )
        stop( "Revision to remove does not exist." )
    DSL:::.revisions(corpus$content) <- getRevisions(corpus)[!pos]
    invisible( corpus )
}

## as.DList.DirSource <- function( x, DStorage = NULL, ... ){
##     if( is.null(DStorage) )
##         DStorage <- DS_default()
##     ## we like to have one file in one chunk so that
##     DStorage$chunksize = 1L
##     as.DList( as.list(x$FileList), DStorage = DStorage, ... )
## }
