print.default <-
function (x, ...) 
{
    class.x <- class(x)
    methods.print <- row.names(attributes(utils::methods("print"))$info)
    classes <- substr(methods.print, 7, nchar(methods.print))
    if (any(class.x %in% classes)) 
        base::print(x, ...)
    else base::print.default(x, ...)
}
