initDict <-
function(pathData = "")
{
    validPath <- FALSE
    for(path in
        c(## Try user supplied path
          pathData,
          ## Try WNHOME (UNIX) environment variable
          file.path(Sys.getenv("WNHOME"), "dict"),
          ## Windows editions provide a registry key
          ## Try UNIX Wordnet 3.0 default path
          "/usr/local/WordNet-3.0/dict",
          ## Try UNIX Wordnet 2.1 default path
          "/usr/local/WordNet-2.1/dict",
          ## Try Debian WordNet default path
          "/usr/share/wordnet"
          )) {
        .jcall("com.nexagis.jawbone.Dictionary", "V", "initialize", path)
        validPath <-
            .jcall("com.nexagis.jawbone.Dictionary", "Z", "pathIsValid")
        if(validPath) break
    }

    if(!validPath)
        warning("cannot find WordNet 'dict' directory: please set the environment variable WNHOME to its parent")

    validPath
}

getDictInstance <-
function()
{
    .jnew("com.nexagis.jawbone.Dictionary")
}

setDict <-
function(pathData)
{
    if(initDict(pathData))
        dict(getDictInstance())
    else
        stop("could not find WordNet installation")
}

getDict <-
function() {
    if(!is.null(d <- dict()))
        d
    else
        stop("could not find Wordnet dictionary")
}

getIndexTerms <-
function(pos, maxLimit, filter)
{
    pos <- .expand_synset_type(pos[1L])
    iterator <-
        .jcall(getDict(),
               "Ljava/util/Iterator;",
               "getIndexTermIterator",
               .jfield("com.nexagis.jawbone.PartOfSpeech",
                       "Lcom/nexagis/jawbone/PartOfSpeech;",
                       pos),
               as.integer(maxLimit),
               .jcast(filter, "com.nexagis.jawbone.filter.TermFilter"))
    .jevalIterator(iterator)
}

WN_synset_types <-
    c("NOUN", "VERB", "ADJECTIVE", "ADJECTIVE_SATELLITE", "ADVERB")

.expand_synset_type <-
function(x)
{
    y <- charmatch(x, WN_synset_types)
    if(is.na(y))
        stop(sprintf("Unknown synset type '%s'", x))
    if(y == 0) {
        if(nchar(x) < 3L)
            stop(sprintf("Ambiguous synset type abbrev '%s'", x))
        if(substring(x, 3L, 3L) == "J")
            "ADJECTIVE"
        else
            "ADVERB"
    } else {
        WN_synset_types[y]
    }
}
