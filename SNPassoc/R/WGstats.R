`WGstats` <-
function (object, ...) 
{
    if (!inherits(object, "WGassociation")) 
        stop("object must be an object of class 'WGassociation'")

    if (!is.null(attr(object,"fast")))
       stop("\n summary is implemented only for 'WGassociation' function")

    x <- attr(object,"tables")
    mostattributes(x)<-NULL
    
    print(x, na.print = "", ...)
    invisible(x)
}

