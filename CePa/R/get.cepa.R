# return cepa object
get.cepa = function(x, id = NULL, cen = 1) {
    
    if(class(x) != "cepa.all") {
        stop("x should be cepa.all object.\n")
    }
    
    if(is.null(id)) {
        stop("id cannot be null")
    }
    
    if(length(cen) > 1) {
        stop("Length of cen must be equal to 1.\n")
    }
    
    if(is.function(cen)) {
        cen = deparse(substitute(cen))
    }
    else if(mode(cen) == "name") {
        cen = deparse(cen)
    }
    
    return(x[[id]][[cen]])
}
    
