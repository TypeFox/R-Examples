grtrnorm <- 
    function(n, 
        np = 2, 
        means = list(rep(0,np), rep(0,np)), 
        covs = diag(rep(1,np)), 
        clip.sd = Inf, 
        tol = 1e-6,
        empirical = TRUE, 
        seed = NULL,
        response.acc = NULL)
{
    if(length(n)!= np) n <- rep(n, np)[1:np]
    if(!is.list(means))	means <- list(means)
    if(length(means) != np) means <- rep(means, np)[1:np]
    dimen <- length(means[[1L]])
    if(!is.list(covs)) covs <- list(covs)
    if(length(covs) != np) covs <- rep(covs, np)[1:np]
    if(length(clip.sd) > 1) clip.sd <- abs(clip.sd[1L])
    
    #Set up random seed if not NULL
    if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)                     # initialize the RNG if necessary
    if(is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
	set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    X <- NULL
    for(i in 1:np)
    {
	X[[i]] <- MASS::mvrnorm(n=n[[i]], mu=means[[i]], Sigma=covs[[i]], tol = tol, empirical=empirical)
	if(is.finite(clip.sd))
	{
	    tmp <- scale(X[[i]])
	    tmp[tmp > clip.sd] <- clip.sd
	    tmp[tmp < -clip.sd] <- -clip.sd
	    X[[i]] <- unscale(tmp)
	}
	X[[i]] <- cbind(rep(i,n[[i]]),X[[i]])
    }
    X <- as.data.frame(do.call("rbind",X))
    label <- paste("x",1:dimen,sep="")
    colnames(X) <- c("category", label)
    if(!is.null(response.acc) & (np > 1)){
        X$response <- X$category
        acc <- sample(1:nrow(X), floor(nrow(X)*(1-response.acc)))
        X$response[acc] <- c(2:np,1)[X$category[acc]]
    }
    X
}

grtMeans <- 
    function(covs, 
         centroid = c(1,1), 
         optldb = c(1,-1,0), 
         p.correct = .85, 
         initd = 5, 
         stepsize = 1)
{
    #if(!is.list(covs)) covs <- replicate(g, list(covs))
    if(!is.vector(centroid)) stop("centroid is not a vector")
    dimen = length(centroid)
    if(!is.vector(optldb)) stop("optldb is not a vector")
    if(inherits(optldb, "glcStruct"))
       optldb <- unlist(optldb[c("coeffs","bias")])
    if(length(optldb)!= dimen+1) 
        stop("length(optldb) must be length(centroid) + 1")
    len <- sqrt(sum(optldb[1:(dimen)]^2))
    a <- optldb[1:(dimen)]/len
    b <- optldb[(dimen+1)]/len
    d <- initd
    d <- d - stepsize
    pc <- 0
    while(pc < p.correct)
    {
	d <- d + stepsize
	means <- list((centroid - d*a), (centroid + d*a))
	pc <- ldb.p.correct(means, covs, noise = 0)
    }
    res <- NULL
    res$means <- means
    res$covs <- covs
    res$p.correct <- pc
    res
}