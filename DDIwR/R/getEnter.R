getEnter <- function(OS) {

    detectedOS <- Sys.info()[['sysname']]
    
    if (OS == "Windows" | OS == "windows" | OS == "Win" | OS == "win") {
        enter <- ifelse(detectedOS == "Windows", "\n", "\r\n")
    }
    else if (OS == "Linux" | OS == "linux") {
        enter <- "\n"
    }
    else if (OS == "Darwin" | OS == "MacOS" | OS == "Apple" | OS == "Mac" | OS == "mac") {
        enter <- ifelse(detectedOS == "Darwin", "\n", "\r")
    }
    else {
        cat("\n")
        stop("The specified OS is not supported.\n\n", call. = FALSE)
    }
    
    return(enter)
}
