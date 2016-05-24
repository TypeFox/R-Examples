"typeToPrototype" <-
function (type, prototype, classes) 
{
    x <- match(type, classes[, 1])
    if (is.null(x) || all(is.na(x))) 
        x <- match(type, classes[, 2])
    if (!is.null(x) && !all(is.na(x))) 
        prototype <- paste(classes[, 2][x])
    return(prototype)
}
