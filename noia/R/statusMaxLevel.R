statusMaxLevel <-
function (effect, max.level = NULL) 
{
    if (!is.null(max.level) && nchar(effect) > max.level) {
        if (sum(strsplit(effect, "")[[1]] != noia::effectsNames[1]) > 
            max.level) 
            return(FALSE)
    }
    return(TRUE)
}
