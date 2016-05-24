`spearman` <-
function(x, y, importance, weights)
{   
    # Inputs:   x - lists to be combined
    #           y - candidate lists
    #           importance - the weight factors indicating the importance of ordered lists
    #           weights - weight matrix if weights to be used
    
	if(is.vector(x))
		x <- matrix(x,1)
	if(is.vector(y))
		y <- matrix(y,1)
	
    k <- ncol(x)
    N <- nrow(y)
   
    res <- rep(0, N)
	
    if(!is.null(weights))
        .C("spearmanCandsW", as.integer(y), as.integer(x), as.integer(N),
				as.integer(k), as.integer(nrow(x)), as.double(importance), 
				as.double(weights), as.double(res), PACKAGE="RankAggreg")[[8]]
    else
		.C("spearmanCands", as.integer(y), as.integer(x), as.integer(N),
				as.integer(k), as.integer(nrow(x)), as.double(importance), 
				as.double(res), PACKAGE="RankAggreg")[[7]]
}

