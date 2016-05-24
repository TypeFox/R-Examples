extract.time = function (data, timeorder) 
{
    x = data$x
    y = data$y
    if (length(timeorder) == 1) {
        index = which(as.numeric(colnames(y)) == timeorder)
	if(length(index) ==0)
	{
	    stop("timeorder must have the same names as the colnames of y")
	}
    }
    if (length(timeorder) > 1) {
	index = vector(, length(timeorder))
        for (i in 1:length(timeorder)) {
            index[i] = which(as.numeric(colnames(y)) == timeorder[i])
        }
    }
    newdata = as.matrix(y[, index])
    ftsdata = fts(x, newdata, start=as.numeric(colnames(y))[1], xname = data$xname, yname = data$yname)
    colnames(ftsdata$y) = as.character(timeorder)
    class(ftsdata) = class(data)
    return(ftsdata)
}
