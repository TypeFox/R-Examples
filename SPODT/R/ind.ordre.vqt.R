ind.ordre.vqt <-
function(ordre.vqt, rownams)
{
    if (is.vector(ordre.vqt))
    {
        return(ind.ordre.vqt.i(ordre.vqt, rownams=rownams))
    }
    else
    {
        return(as.vector(apply(ordre.vqt, MARGIN=2, ind.ordre.vqt.i, rownams=rownams)))
    }
}
