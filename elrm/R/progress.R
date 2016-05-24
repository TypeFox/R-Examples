`progress` <-
function (value, max.value = NULL) 
{
    if (!is.numeric(value)) 
        stop("`value' must be numeric!")
    if (is.null(max.value)) {
        max.value <- 100
        percent <- TRUE
    }
    else percent <- FALSE
    if (!is.numeric(max.value)) 
        stop("`max.value' must be numeric or NULL!")
    erase.only <- (value > max.value)
    max.value <- as.character(round(max.value))
    l <- nchar(max.value)
    value <- formatC(round(value), width = l)
    if (percent) {
        backspaces <- paste(rep("\b", l + 14), collapse = "")
        if (erase.only) 
            message <- ""
        else message <- paste("Progress: ", value, "%  ", sep = "")
        cat(backspaces, message, sep = "")
    }
    else {
        backspaces <- paste(rep("\b", 2 * l + 16), collapse = "")
        if (erase.only) 
            message <- ""
        else message <- paste("Progress: ", value, " on ", max.value, 
            "  ", sep = "")
        cat(backspaces, message, sep = "")
    }
    if (.Platform$OS.type == "windows") 
        flush.console()
    invisible(NULL)
}

