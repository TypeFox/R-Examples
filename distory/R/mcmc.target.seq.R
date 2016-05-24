mcmc.target.seq <- function(data, x, F, n) 
{
    N = ncol(data)
    D = rep(1, N)

    likvals = c()

    a = runif(n)

    old.tree = F(lookup.samples(data, list(convert.table.to.idx(D)))[[1]])
    old.dist = as.matrix(dist.multiPhylo(list(x, old.tree)))[1,2]

    best.tree = old.tree
    best.dist = old.dist

    for(i in 1:n)
    {
        Dp = D
        good = FALSE

        r = sample.int(N, 2)
        while(D[r[1]] == N || D[r[2]] == 0)
        {    
            r = sample.int(N, 2)
        }

        # propose the swap
        Dp[r[1]] = Dp[r[1]] + 1;
        Dp[r[2]] = Dp[r[2]] - 1;

        new.tree = F(lookup.samples(data, list(convert.table.to.idx(Dp)))[[1]])
        new.dist = as.matrix(dist.multiPhylo(list(x, new.tree)))[1,2]

        v = old.dist/new.dist

        if(v >= 1) 
        { 
            v = 1 
            best.tree = new.tree
            best.dist = new.dist
        } 

        if(v^(100-exp(1/i)) > a[i])
        {
            D = Dp
            old.tree = new.tree
            old.dist = new.dist
        }

        likvals = c(likvals, old.dist)
    }
    
    list(data = D, tree = best.tree, distance = best.dist, vals = likvals)
}

# modified from ape code for boot.phylo
boot.samples.idxs <- function(data, B = 100, block = 1)
{
    if (is.list(data) && !is.data.frame(data)) {
        if (inherits(data, "DNAbin")) 
            data <- as.matrix(data)
        else {
            nm <- names(data)
            n <- length(data)
            data <- unlist(data)
            nL <- length(data)
            data <- matrix(data, n, nL/n, byrow = TRUE)
            rownames(data) <- nm
        }
    }
    boot.smpls <- vector("list", B)
    for (i in 1:B) {
        if (block > 1) {
            y <- seq(block, ncol(data), block)
            boot.i <- sample(y, replace = TRUE)
            boot.samp <- numeric(ncol(data))
            boot.samp[y] <- boot.i
            for (j in 1:(block - 1)) boot.samp[y - j] <- boot.i - 
                j
        }
        else boot.samp <- sample(ncol(data), replace = TRUE)

        boot.smpls[[i]] <- boot.samp
    }

    boot.smpls
}

lookup.samples <- function(data, idxs)
{
    if (is.list(data) && !is.data.frame(data)) {
        if (inherits(data, "DNAbin")) 
            data <- as.matrix(data)
        else {
            nm <- names(data)
            n <- length(data)
            data <- unlist(data)
            nL <- length(data)
            data <- matrix(data, n, nL/n, byrow = TRUE)
            rownames(data) <- nm
        }
    }

    lapply(idxs, function(i) { data[,i] } )
}

phylo.diff <- function(x, y, ...)
{
    uniqT1 <- distinct.edges(x, y)
    uniqT2 <- distinct.edges(y, x)
    
    treeA.cs <- rep("black", dim(x$edge)[1]) 
    treeA.cs[uniqT1] = "red"

    treeB.cs <- rep("black", dim(y$edge)[1]) 
    treeB.cs[uniqT2] = "red"

    par(mfrow=c(1,2))
    plot(x, edge.color=treeA.cs, ...)
    plot(y, edge.color=treeB.cs, ...)

    invisible()
}

convert.table.to.idx <- function(T)
{
    dat = c()
    for(i in 1:length(T))
    {
        dat = c(dat, rep(i, T[i]))
    }

    dat
}

