nb.feuilles <-
function(noeud)
{
    nb.feuilles <- 0
    if (class(noeud) == "f.spodt")
    {
        nb.feuilles <- nb.feuilles + 1
    }
    else
    {
        nb.feuilles <- nb.feuilles(noeud@fg) + nb.feuilles(noeud@fd)
    }
    return(nb.feuilles)
}
