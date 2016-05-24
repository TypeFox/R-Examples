extract.x = function (data, xorder) 
{
    x = data$x
    y = data$y
    if (length(xorder) < 2) {
	    stop("xorder must be multiple x variable")
    }
    if (length(xorder) >= 2) {
	index = vector(, length(xorder))
        for (i in 1:length(xorder)) {
            index[i] = which(as.numeric(rownames(y)) == xorder[i])
        }
    }
    newdata = as.matrix(y[index, ])
    newx = x[index]
    ftsdata = fts(newx, newdata, start = as.numeric(rownames(y))[1], xname = data$xname, yname = data$yname)
    rownames(ftsdata$y) = as.character(xorder)
    class(ftsdata) = class(data)
    return(ftsdata)
}