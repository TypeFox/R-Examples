top <-
function (O, neighbors, n) 
{
    temp = as.matrix(apply(neighbors, 1, sum))
    out = rbind(O, temp)
    out = as.matrix(out)
    rowsn = as.numeric(rownames(out))
    rowsno = rowsn[order(out, decreasing = TRUE)]
    out.sort = as.matrix(out[order(-out)])
    rownames(out.sort) = rowsno
    outliers = as.matrix(out.sort[1:n, ])
    return(outliers)
}
