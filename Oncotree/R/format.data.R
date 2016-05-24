"formatData" <-
function(dataset){
    keep <- apply(dataset, 2, function(x) diff(range(x))>0)
    dropcol <- apply(dataset, 2, function(x) diff(range(x))==0)
    droplist <- colnames(dataset[,dropcol])
    if (length(droplist)>0) {
        cat("The following events had no observed occurances,")
        cat("so they will not be included in the construction of the tree:\n")
        cat(droplist, "\n")}
    datab <- dataset[, keep,drop=FALSE]  
    Root <- rep(1, nrow(datab))
    data <- cbind(Root, datab)
    if (!is.matrix(data)) data <- as.matrix(data)
  	if (!(all(data %in% c(0,1))))
	    stop("All the events should be coded 0 - 1.")
  	if (is.null(colnames(data)))
	    stop("The data set should have column names.")
    return(data)
}

