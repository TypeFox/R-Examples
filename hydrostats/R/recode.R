recode <- function(x, oldvalue, newvalue) {
    
    # convert any factors to characters
    
    if (is.factor(x)) 
        x <- as.character(x)
    if (is.factor(oldvalue)) 
        oldvalue <- as.character(oldvalue)
    if (is.factor(newvalue)) 
        newvalue <- as.character(newvalue)
    
    # create the return vector
    
    newvec <- x
    
    # put recoded values into the correct position in the return vector
    
    for (i in unique(oldvalue)) newvec[x == i] <- newvalue[oldvalue == i]
    
    return(newvec)
    
} 
