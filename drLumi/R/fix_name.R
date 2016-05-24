fix_name <- function(x, values = c("&","'","@","*")){
    for (i in 1:length(values)){
        x <- gsub(values[i], "_", x,fixed=TRUE)
    }
    return(x)
}




