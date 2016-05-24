setMethod("show", "asAnscombe", 
    function(object){
        cat(paste("An object of class", dQuote(class(object)), "\n"))
        cat("risk type:\t", object@type, "\n")
        cat("ARE in the ideal model:\t", object@eff, "\n")
    })
