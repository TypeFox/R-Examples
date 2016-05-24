Split.CH <-
function (ch = ch) 
{
    chmat = t(sapply(strsplit(ch[, 1], ""), function(x) as.numeric(x)))
    return(chmat)
}
