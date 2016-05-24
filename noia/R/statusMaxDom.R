statusMaxDom <-
function (effect, max.dom = NULL) 
{
    if (!is.null(max.dom) && nchar(effect) > max.dom) {
        if (sum(strsplit(effect, "")[[1]] == noia::effectsNames[3]) > 
            max.dom) 
            return(FALSE)
    }
    return(TRUE)
}
