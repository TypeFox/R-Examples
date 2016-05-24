define_data <- function(x){
    for(i in 1:ncol(x)) x[,i] <- as.character(x[,i])
    var_names <- as.character(unlist(x[1,]))
    nna <- sum(var_names=="" | is.na(var_names))
    new_data <- x[-1,]
    if(nna == 0){
        names(new_data) <- var_names    
    }
    if(nna != 0){
        var_names <- var_names[var_names!="" & !is.na(var_names)]
        new_data <- new_data[,1:length(var_names)]
        names(new_data) <- var_names    
    }  
    new_data
}

