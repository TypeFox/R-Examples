localRFiles <- function(directory = getwd(), suffix = "[.][RSq]$", ask = FALSE) {
    files <- list.files(directory)
    which <- grep(suffix, files)
    if(ask) {
        files <- files[which]
        which <- menu(c(files, "Other?"), graphics = TRUE)
        if(which > length(files)) {
            message("Enter a file name:")
            return(scan(stdin(), "",n=1))
        }
    }
    files[which]
}

    
    
    
