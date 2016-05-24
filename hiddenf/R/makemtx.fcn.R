makemtx.fcn <-
function(tall)
{
tapply(tall$y,list(tall$rows,tall$cols),identity)
}
