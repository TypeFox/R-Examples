#' @export
print.lum_import <- function(x, ...){
    cat(paste("The file imported is:",x$type_raw_data,"type","\n",sep=" ") )
    if(x$type_raw_data=="Fluorescence"){
        cat(paste("Identified data type:", names(x$dtblock),"\n",sep=" "
                  , collapse="") )
    }
    if(x$type_raw_data=="Bead"){
        cat(paste("Number of unique wells files", 
            length(unique(x[[1]]$well)),"\n",sep=" ") )
    }  
    cat("\n")
}

