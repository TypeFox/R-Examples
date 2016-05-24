position.leverage <- function(data, F, to = NULL, rep = 50, by = 1)
{
    N = ncol(data)
    D = rep(1, N)

    leverage.vals = c()

    if(is.null(to))
    {
        to = F(lookup.samples(data, list(convert.table.to.idx(D)))[[1]])
    }

    for(i in seq(1,N,by))
    {
        Dp = D;
        S = sample(N, N/rep);
        while(i %in% S)
        {
            S = sample(N,N/rep);
        }

        Dp[S] = 0;
        Dp[i] = N/rep;

        new.tree = F(lookup.samples(data, list(convert.table.to.idx(Dp)))[[1]])
        leverage.vals = c(leverage.vals, as.matrix(dist.multiPhylo(list(to,
                            new.tree)))[1,2])
    }

    leverage.vals
}


