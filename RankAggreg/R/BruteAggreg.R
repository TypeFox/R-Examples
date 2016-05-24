`BruteAggreg` <-
function(x, k, weights=NULL, distance=c("Spearman", "Kendall"),
         importance=rep(1,nrow(x))){
    distance <- match.arg(distance, c("Spearman", "Kendall"))
    x <- x[,1:k]
    orig.x <- x

    orig.imp <- importance
    importance <- importance/sum(importance) #rescale importance weights

    distinct <- apply(x, 1, function(y) ifelse(length(unique(y)) < k, 1, 0))
    if (sum(distinct) >= 1)
        stop("Elements of Each Row Must Be Unique")
    if (nrow(x)<2)
        stop("X must have more than 1 row") 
    
    if(!is.null(weights)){
        weights <- weights[,1:k]
        #standardize weights:
        weights <- t(apply(weights,1,function(z){if(max(z)==min(z)) rep(0, length(z))
        	else (z-min(z))/(max(z)-min(z))}))
	  for(i in 1:nrow(weights))
        	if(weights[i,k]!=0)
            	weights[i,] <- 1-weights[i,]
        if(dim(x)[1] != dim(weights)[1] || dim(x)[2] != dim(weights)[2])
        	stop("Dimensions of x and weights matrices have to be the same")
    }

       
    comp.list <- unique(sort(as.vector(x)))
    n <- length(comp.list)
	
    x <- t(apply(x,1, function(xx) match(xx,comp.list)))

    if (k > n)
        stop("k must be smaller or equal to n") 

    perms <- permutations(n,k,1:n)
    if(distance=="Spearman")
        f.y <- spearman(x, perms, importance, weights)
    else
        f.y <- kendall(x, perms, importance, weights)
    
	tl <- comp.list[perms[which.min(f.y),]]
    res <- list(top.list=tl,optimal.value=min(f.y), distance=distance,
		method="BruteForce", importance=orig.imp, lists=orig.x, weights=weights, 
		sample=f.y, sample.size=length(f.y), summary=matrix(c(min(f.y),median(f.y)),1,2))
    class(res) <- "raggr"
    res
}

