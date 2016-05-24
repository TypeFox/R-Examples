## Paste an Unspecified Number Of Text Columns
## 
## \code{paste2} - Paste unspecified columns or a list of vectors together.
paste2 <-
function(multi.columns, sep=".", handle.na=TRUE, trim=TRUE){
    if (is.matrix(multi.columns)) {
        multi.columns <- data.frame(multi.columns)
    }
    if (trim) multi.columns <- lapply(multi.columns, function(x) {
            gsub("^\\s+|\\s+$", "", x)
        }
    )
    if (!is.data.frame(multi.columns) & is.list(multi.columns)) {
        multi.columns <- do.call('cbind', multi.columns)
    } 
    m <- if (handle.na){
                 apply(multi.columns, 1, function(x){
                     if (any(is.na(x))){
                         NA
                     } else {
                         paste(x, collapse = sep)
                     }
                 }
             )   
         } else {
             apply(multi.columns, 1, paste, collapse = sep)
    }
    names(m) <- NULL
    return(m)
}




## Remove Leading/Trailing White Space
## 
## Remove leading/trailing white space.
Trim <-
function (x) gsub("^\\s+|\\s+$", "", x)

hijack <- function(FUN, ...){

    .FUN <- FUN

    args <- list(...)
    invisible(lapply(seq_along(args), function(i) {
        formals(.FUN)[[names(args)[i]]] <<- args[[i]]
    }))
    
    .FUN
}
