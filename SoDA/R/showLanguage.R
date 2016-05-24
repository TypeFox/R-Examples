showLanguage <- function(object, indent = "") {
    switch(mode(object),
           call = , "(" = 
           showCall(object, indent),
           name = 
           cat(indent, "Name object: ", dQuote(deparse(object)), "\n", sep=""),
       {##  non-language, presumably, since "language" is sealed.
           ## normally will be scalar 
           cat(indent)
           if(length(object)>1)
               showDefault(object)
           else {
               cat("An object of class ",dQuote(class(object)), "\n",
                   indent, deparse(object), "\n", sep="")
           }
       }
           )
    invisible(object)
}

showCall <- function(object, indent = "") {
    nargs <- length(object)-1
    if(nargs < 0)
        stop("Invalid call to showCall: object of class ", dQuote(class(object)),
             " and length 0")
    callTo <- object[[1]]
    cat(indent, "A call object of class ", dQuote(class(object)), " with ", nargs,
        ifelse(identical(nargs, 1), " argument\n"," arguments\n"), sep="")
    indent2 <- paste(indent, "    ", sep="")
    if(is(callTo, "name"))
        cat(indent, "Call to ", dQuote(deparse(callTo)), "\n", sep="")
    else {
        cat(indent, "Call supplied as object of class ", dQuote(class(callTo)),
            "\n", sep="")
        showLanguage(callTo, indent2)
    }
    for(i in seq(length = nargs)) {
        argi <- object[[i+1]]
        cat(indent, "Argument number ", i, "\n", sep="")
        showLanguage(argi, indent2)
    }
}

setMethod("show", "language", function(object)showLanguage(object))
