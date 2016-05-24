synonyms <-
function(word, pos)
{
    filter <- getTermFilter("ExactMatchFilter", word, TRUE)
    terms <- getIndexTerms(pos, 1L, filter)
    if (is.null(terms))
        character()
    else
        getSynonyms(terms[[1L]])
}
