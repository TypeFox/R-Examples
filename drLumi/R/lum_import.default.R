#' @export
lum_import.default <- function(x, ...){
    if(!file.exists(x)) stop("Not found 'x' file")
    n <- list.files(x)
    nch <- nchar(basename(x))
    extension <- substring(basename(x), first = nch-3, last=nch)
    if(length(n)==0 & extension==".csv"){
        ans <- lum_import_fluor(x)
        ans$type_raw_data <- "Fluorescence"
    }    
    if(length(n)>0 | extension==".zip"){
        ans <- lum_import_bead(x)
        ans$type_raw_data <- "Bead"
    }
    if(length(n)==0 & extension%nin%c(".csv",".zip") ){
        stop("'x' must be a '.csv' path file")
    } 
    class(ans) <- "lum_import"
    ans
}

