part.ordre.vqt <-
function(ordre.vqt=ordre.vqt, rownams)
{
    if (is.vector(ordre.vqt))
    {
        return(part.ordre.vqt.i(ordre.vqt, rownams=rownams))
    }
    else
    {
        return(apply(ordre.vqt, MARGIN=2, part.ordre.vqt.i, rownams=rownams))
    }
}
