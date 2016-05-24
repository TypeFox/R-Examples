setMethod("show", "RandVariable", 
    function(object){
        cat(paste("An object of class", dQuote(class(object)), "\n"))
        cat("length of Map:\t", length(object), "\n")
        cat("Domain:\t")                
        if(is.null(object@Domain))
            print(NULL)
        else{
            cat(object@Domain@name, "with dimension ")
            cat(object@Domain@dimension, "\n")
        }
        cat("Range:\t")
        if(is.null(object@Range))
            print(NULL)
        else{
            cat(object@Range@name, "with dimension ")
            cat(object@Range@dimension, "\n")
        }
    })

setMethod("show", "EuclRandMatrix", 
    function(object){
        cat(paste("An object of class", dQuote(class(object)), "\n"))
        cat("Dim of Map:\t", Dim(object), "\n")
        cat("Domain:\t")                
        if(is.null(object@Domain))
            print(NULL)
        else{
            cat(object@Domain@name, "with dimension ")
            cat(object@Domain@dimension, "\n")
        }
        cat("Range:\t")
        cat(object@Range@name, "with dimension ")
        cat(object@Range@dimension, "\n")
    })
setMethod("show", "EuclRandVarList", 
    function(object){
        cat(paste("An object of class", dQuote(class(object)), "\n"))
        cat("Domain:\t")                
        if(is.null(object[[1]]@Domain))
            print(NULL)
        else{
            cat(object[[1]]@Domain@name, "with dimension ")
            cat(object[[1]]@Domain@dimension, "\n")
        }
        for(i in 1:length(object)){
            cat("[[", i, "]]\n", sep = "")
            if(is(object[[i]], "EuclRandMatrix"))
                cat("Dim of Map:\t", Dim(object[[i]]), "\n")
            else
                cat("length of Map:\t", length(object[[i]]), "\n")
            cat("Range:\t")
            cat(object[[i]]@Range@name, "with dimension ")
            cat(object[[i]]@Range@dimension, "\n")
        }
    })
