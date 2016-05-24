coor <-
function (a, b) 
{
    out = list(net = a, coord = b)
    class(out) = c("GML")
    out
}
