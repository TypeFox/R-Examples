`createString` <-
function(mydata, collapse="*", uplow=FALSE, use.tilde=FALSE) {
    
    requireNamespace("QCA")
    
    mydata <- changemydata <- as.matrix(mydata)
    conditions <- colnames(mydata)
    if (uplow) {
        changemydata[mydata == 0] <- tolower(rep(conditions, each=nrow(mydata))[mydata == 0])
        changemydata[mydata == 1] <- toupper(rep(conditions, each=nrow(mydata))[mydata == 1])
    }
    else if (use.tilde) {
        changemydata[mydata == 0] <- paste("~", toupper(rep(conditions, each=nrow(mydata))[mydata == 0]), sep="")
        changemydata[mydata == 1] <- toupper(rep(conditions, each=nrow(mydata))[mydata == 1])
    }
    else {
        for (i in sort(unique(as.vector(mydata)))) {
            changemydata[mydata == i] <- paste(rep(conditions, each=nrow(mydata))[mydata == i], "{", i, "}", sep="")
        }
    }
    
    input <- rep(NA, nrow(mydata))
    
    for (i in 1:nrow(mydata)) {
        input[i] <- paste(changemydata[i, ], collapse = collapse)
    }
    return(input)
}

