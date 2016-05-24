"plot.car" <-
function (x, type=c("spec", "pred", "diag"), ...)  
{
    if(class(x) != "car")
    stop("Not a 'car' object\n")
    type <- match.arg(type)
    if(type=="spec")
    plot(x, ...)
    #plot.spec(x, ...)
    if(type=="pred")
    plot.predict.car(x, ...)
    if(type=="diag")
    diag(x, ...)
    invisible(x)
}
