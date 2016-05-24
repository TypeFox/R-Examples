complexite <-
function(noeud, cplx=1)
{
    if (class(noeud) == "sp.spodt"  |  class(noeud) == "vql.spodt"  |  class(noeud) == "vqt.spodt")
    {
        cplx <- max(complexite(noeud@fg, cplx), complexite(noeud@fd, cplx)) + 1
    }
    return(cplx)
}
