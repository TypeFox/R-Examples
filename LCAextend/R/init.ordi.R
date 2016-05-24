init.ordi <-
function(y,K,x=NULL,var.list=NULL)
{
    if(nrow(y)<K) stop("there are a few individuals to perform a model with ",K," classes\n")
    classes <- rep(1:K,each=nrow(y)/K)
    if(length(classes)<nrow(y)) classes[(length(classes)+1):nrow(y)] <- K
    alpha <- list()
    for(j in 1:ncol(y))
    {
        level <- as.numeric(levels(factor(y[,j])))
		level <- 1:max(level)

		alpha[[j]] <- matrix(NA,nrow=K,ncol=length(level)-1)
        for(k in 1:K)
        {
            per.level <- as.vector(table(c(level,y[classes==k,j]))-1)
            per.level <- per.level/sum(per.level)
            alpha[[j]][k,] <- alpha.compute(per.level)
        }
        if(!is.null(var.list[[j]])) alpha[[j]] <- cbind(alpha[[j]],matrix(0,nrow=K,ncol=length(var.list[[j]])))
    }
    res <- list("alpha"=alpha)
    res
}

