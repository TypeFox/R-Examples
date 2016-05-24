printPercentage <-
function (i, n, dp = 0, first = 1, last = n, prev = i - 1) 
{
    out = ""
    disp = round(100 * i/n, dp)
    if (prev >= first) 
        prev.disp = round(100 * prev/n, dp)
    else prev.disp = ""
    if (disp > prev.disp) {
        nc = nchar(prev.disp)
        if (i != first) {
            out = paste(out, paste(rep("\b", nc + 1), collapse = ""), 
                sep = "")
        }
        out = paste(out, disp, "%", sep = "")
    }
    if (i == last) 
        out = paste(out, "\n", sep = "")
    cat(out)
    return(NULL)
}
