part.ordre.vqt.i <-
function(ordre.vqt.i, rownams)
{
    return(ordre.vqt.i[which(ordre.vqt.i %in% rownams == TRUE)])
}
