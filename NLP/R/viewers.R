words <-
function(x, ...)
    UseMethod("words")

sents <-
function(x, ...)
    UseMethod("sents")

paras <-
function(x, ...)
    UseMethod("paras")

tagged_words <-
function(x, ...)
    UseMethod("tagged_words")

tagged_sents <-
function(x, ...)
    UseMethod("tagged_sents")

tagged_paras <-
function(x, ...)
    UseMethod("tagged_paras")

chunked_sents <-
function(x, ...)
    UseMethod("chunked_sents")

parsed_sents <-
function(x, ...)
    UseMethod("parsed_sents")

parsed_paras <-
function(x, ...)
    UseMethod("parsed_paras")


chunk_tree_from_chunk_info <-
function(words, ptags, ctags)
{
    ind <- grepl("^[BO]", ctags)
    ## <FIXME>
    ## Should this also use Tagged_Token()?
    chunks <- split(sprintf("%s/%s", words, ptags), cumsum(ind))
    ## </FIXME>
    nms <- sub(".*-", "", ctags[ind])
    ind <- nms != "O"
    chunks[ind] <- Map(Tree, nms[ind], chunks[ind])
    Tree("S", chunks)
}

POS_tag_mapper <-
function(map, set)
{
    if(is.function(map))
        return(map)
    if(is.list(map))
        map <- map[[set]]
    function(pos) map[pos]
}
