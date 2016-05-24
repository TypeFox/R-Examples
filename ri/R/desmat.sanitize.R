desmat.sanitize <-
function(  Z,
                        blockvar = NULL,
                        clustvar = NULL)
                        {
        # To make blockvar sequence of integers from 1 to number of blocks
        if (!is.null(blockvar)){
            if(is.character(blockvar)){blockvar <- as.numeric(as.factor(blockvar))}
            if(is.factor(blockvar)){blockvar <- as.numeric(blockvar)}
            blockvar <- match(blockvar,unique(blockvar))
        }
        if (is.null(blockvar)) blockvar <- rep(1,length(Z))
        B <- max(blockvar)

        # To make clustvar sequence of integers from 1 to number of clusters
        if (!is.null(clustvar)){
            if(is.character(clustvar)){clustvar <- as.numeric(as.factor(clustvar))}
            if(is.factor(clustvar)){clustvar <- as.numeric(clustvar)}
            clustvar.sc <- clustvar
            clus.count <- 1
            for(i in 1:B){
                clus.b.i <- unique(clustvar[blockvar==i])
                for(j in clus.b.i){
                    clustvar.sc[blockvar==i&clustvar==j] <- clus.count
                    clus.count <- clus.count+1
                    }
                }
            clustvar <- clustvar.sc
        }
        if (is.null(clustvar)) clustvar <- c(1:length(Z))

        Z <- aggregate(Z,list(clustvar),mean)[,2]
        #ave(Z, clustvar)
        blockvar <- aggregate(blockvar,list(clustvar),mean)[,2]
        #ave(blockvar, clustvar)

        desmat <- data.frame(n=NA,m=NA, rand=NA)
        for(b in 1:B){
            desmat[b,"n"]    <- length(Z[blockvar==b])
            desmat[b,"m"]    <- sum(Z[blockvar==b])
            desmat[b,"rand"] <- choose(length(Z[blockvar==b]),sum(Z[blockvar==b]))
        }
    return(list(desmat=desmat,Z=Z,blockvar=blockvar,clustvar=clustvar))
}
