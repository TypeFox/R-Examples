## Function from the slouch package. Author : Jason Pienaar

'.ouch2slouch.mvsl'<-function (tree) 
{
    if (!inherits(tree, "ouchtree")) 
	stop(sQuote("tree"), " must be of class ", sQuote("ouchtree"))
    N <- length(tree@nodes)
    tmp <- as(tree, "data.frame")
    tmp$ancestors <- as.character(tmp$ancestors)
    tmp$ancestors <- as.numeric(tmp$ancestors)
    tmp$times <- as.character(tmp$times)
    tmp$times <- as.numeric(tmp$times)
    tmp$nodes <- as.character(tmp$nodes)
    tmp$nodes <- as.numeric(tmp$nodes)
    tmp$ancestors[1] <- 0
    slouch_node <- 1:N
    ancestor <- rep(NA, times = N)
    for (i in 1:N) {
        if (tmp$labels[i] == "" || is.na(tmp$labels[i])) 
		tmp$labels[i] = NA
    }
    rownames(tmp) <- 1:nrow(tmp)
    names(tmp) = c("nodes", "ancestor", "time", "species")
    return(tmp)
}
