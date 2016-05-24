Collapse.CH <-
function (ch = ch$ch_split) 
{
    ch_vec = apply(ch, 1, paste, collapse = "")
    ch = data.frame(ch = ch_vec)
    ch["ch"] = as.character(ch$ch)
    return(ch)
}
