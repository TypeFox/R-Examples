

#' @method str editmatrix
#' @rdname editmatrix
#' @export
str.editmatrix <- function(object,...){
    vars <- paste(getVars(object),collapse=", ")
    if (nchar(vars) > 20 ) vars <-  paste(strtrim(vars,16),"...") 
    cat(paste("editmatrix with", nrow(object), "edits containing variables",vars,"\n"))
}


#' @method str editarray
#' @export
str.editarray <- function(object,...){
    vars <- paste(getVars(object),collapse=", ")
    if (nchar(vars) > 20 ) vars <-  paste(strtrim(vars,16),"...") 
    cat(paste("editarray with", nrow(object), "edits containing variables",vars,"\n"))
}





