## Function written to match MATLAB function
## Author: Andrew Hooker

fullfile <- function(pathname,filename){
    if(is.null(pathname)){
        return(filename)
    } else {
        return(file.path(pathname,filename))
    }
}
