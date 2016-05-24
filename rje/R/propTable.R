propTable <-
function (x, margin = NULL) 
{
    if (length(margin)) 
        .sweep0(x, margin, marginTable(x, margin), "/")
    else x/sum(x)
}
