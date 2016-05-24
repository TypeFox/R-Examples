'label<-' <-
function(x, value)
{
    if (!is.character(value))
        value <- deparse(value)
    attr(x, "label") <- value
    return(x)
}

label <-
function(x)
{
    attr(x, "label")
}
