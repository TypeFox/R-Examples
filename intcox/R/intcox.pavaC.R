intcox.pavaC <-
function (w, y)
{
    n <- length(y)
    index <- rep(0, n)
    weight <- rep(0, n)
    ghat <- rep(0, n)
    .C("pavaC", as.double(w), as.double(y), as.integer(n), as.integer(index),
        as.double(weight), as.double(ghat), PACKAGE = "intcox" )[[6]]
}
