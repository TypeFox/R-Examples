showEnv <- function(x){
    if(!is.environment(x))
      x <- as.environment(x)
    if(is.null(attr(x, "name"))) print(x)
    else message("<environment: ", attr(x, "name"), ">")
    invisible(x)
} 
