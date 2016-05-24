redl <-
function (a, b) 
{
    red = list(full = a, reduc = b)
    class(red) = c("Galois", "reduced")
    red
}
