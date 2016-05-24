getLemma <-
function(indexterm)
{
    .jcall(indexterm, "S", "getLemma")
}

getSynsets <-
function(indexterm)
{
    .jcall(indexterm, "[Lcom/nexagis/jawbone/Synset;", "getSynsets")
}

getSynonyms <-
function(indexterm)
{
    synsets <-
        .jcall(indexterm, "[Lcom/nexagis/jawbone/Synset;", "getSynsets")
    sort(unique(unlist(lapply(synsets, getWord))))
}
