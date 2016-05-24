"bugs.inits" <- 
function (inits, n.chains, digits,
          inits.files = paste("inits", 1:n.chains, ".txt", sep = "")){
    if(!is.null(inits)) {
        for(i in 1:n.chains) {
            if(is.function(inits))
                    write.datafile(lapply(inits(), formatC, digits = digits, format = "E"),
                        inits.files[i])
            else 
                write.datafile(lapply(inits[[i]], formatC, digits = digits, format = "E"), 
                    inits.files[i])
        }
    }
    return(inits.files)
}
