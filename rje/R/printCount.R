printCount <-
function (i, first = 1, prev = i - 1, last = NULL) 
{
    out = ""
    disp = round(i)
    if (prev >= first) 
        prev.disp = round(prev)
    else prev.disp = ""
    if (disp > prev.disp) {
        nc = nchar(prev.disp)
        if (i != first) {
            out = paste(out, paste(rep("\b", nc), collapse = ""), 
                sep = "")
        }
        out = paste(out, disp, sep = "")
    }
    if (!is.null(last) && i == last) 
        out = paste(out, "\n", sep = "")
    cat(out)
    return(NULL)
}
