getRelatedSynsets <-
function(synset, pointerSymbol)
{
    l <- .jcall(synset, "Ljava/util/List;", "getRelatedSynsets",
                pointerSymbol)
    iterator <- .jcall(l, "Ljava/util/Iterator;", "iterator")
    i <- .jevalIterator(iterator)
    lapply(i, .jcast, "Lcom/nexagis/jawbone/Synset;")
}

getWord <-
function(synset)
{
    l <- .jcall(synset, "Ljava/util/List;", "getWord")
    iterator <- .jcall(l, "Ljava/util/Iterator;", "iterator")
    i <- .jevalIterator(iterator)
    i <- lapply(i, .jcast, "Lcom/nexagis/jawbone/WordData;")
    sapply(i, .jcall, "S", "getWord")
}
