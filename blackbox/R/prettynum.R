prettynum <-
function (x, extradigits = 0)
{
    sapply(x, myprettyW, extradigits = extradigits)
}
