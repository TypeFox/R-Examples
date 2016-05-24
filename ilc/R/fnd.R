fnd <-
function (x, ...) 
{
    if (mode(x) == "character" && length(x) == 1) 
        find(x, simple.words = F)
    else find(deparse(substitute(x)), ...)
}
