.escapeRegex <- function(x) 
{
    return(gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", x))
}
