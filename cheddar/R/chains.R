.TrophicChainsSize <- function(community)
{
    # Returns two integers: n.chains, longest

    if(!is.Community(community)) stop('Not a Community')

    max.queue <- .maxQueueOption()
    stopifnot(max.queue>=0)

    # Get the food web as an adjancency list
    alist <- .CAdjacencyList(community, ConsumersByNode(community))

    # 1 if basal, 0 otherwise
    is.basal <- as.integer(IsBasalNode(community))

    # Outputs to be filled by the C function
    n.chains <- as.integer(0)
    longest <- as.integer(0)
    status <- as.integer(-1)

    res <- .C('trophic_chains_size', 
              as.integer(alist), 
              as.integer(length(alist)), 
              as.integer(is.basal), 
              as.integer(nrow(alist)), 
              as.integer(0),  # test_overflow
              as.integer(max.queue),
              n.chains=n.chains,
              longest=longest,
              status=status, 
              PACKAGE='cheddar', NAOK=TRUE)

    if(-1==res$status)
    {
        stop('Unexpected error')
    }
    else if(0==res$status)
    {
        stopifnot((0==res$n.chains && 0==res$longest) || 
                  (0<res$n.chains && 0<res$longest))
        stopifnot(res$longest <= NumberOfNodes(community))
        return (c(n.chains=res$n.chains, longest=res$longest))
    }
    else if(1==res$status)
    {
        stop('Problem with an input parameter')
    }
    else
    {
        stop(paste('Unknown status [', res$status, ']', sep=''))
    }
}

.PrintChains <- function(community)
{
    if(!is.Community(community)) stop('Not a Community')

    max.queue <- .maxQueueOption()

    # Get the food web as an adjancency list
    alist <- .CAdjacencyList(community, ConsumersByNode(community))

    # 1 if basal, 0 otherwise
    is.basal <- as.integer(IsBasalNode(community))

    # Outputs to be filled by the C function
    status <- as.integer(-1)

    res <- .C('print_chains', 
              as.integer(alist), 
              as.integer(length(alist)), 
              as.integer(is.basal), 
              as.integer(nrow(alist)), 
              as.integer(max.queue),
              status=status, 
              PACKAGE='cheddar', NAOK=TRUE)

    if(-1==res$status)
    {
        stop('Unexpected error')
    }
    else if(0==res$status)
    {
        # Normal exit
    }
    else if(1==res$status)
    {
        stop('Problem with an input parameter')
    }
    else
    {
        stop(paste('Unknown status [', res$status, ']', sep=''))
    }
}

TrophicChainsStats <- function(community)
{
    if(!is.Community(community)) stop('Not a Community')

    # Small performance increase by removing cannibalistic links
    community <- RemoveCannibalisticLinks(community)
    if(0==NumberOfTrophicLinks(community))
    {
        return (NULL)
    }

    max.queue <- .maxQueueOption()
    stopifnot(max.queue>=0)

    # Get the food web as an adjancency list
    alist <- .CAdjacencyList(community, ConsumersByNode(community))

    # 1 if basal, 0 otherwise
    is.basal <- as.integer(IsBasalNode(community))

    # Outputs to be filled by the C function
    # Compute the number of chains and the longest chain
    chains_size <- .TrophicChainsSize(community)
    n.chains <- chains_size['n.chains']
    longest <- chains_size['longest']
    if(n.chains>0)
    {
        n <- NumberOfNodes(community)
        node.pos.counts <- as.integer(rep(NA, n*longest))
        chain.lengths <- as.integer(rep(NA, n.chains))
        status <- as.integer(-1)

        res <- .C('trophic_chains_stats', 
                  as.integer(alist), 
                  as.integer(length(alist)), 
                  as.integer(is.basal), 
                  as.integer(n),
                  as.integer(n.chains),
                  as.integer(longest), 
                  as.integer(max.queue),
                  node.pos.counts=node.pos.counts,
                  chain.lengths=chain.lengths,
                  status=status, 
                  PACKAGE='cheddar', NAOK=TRUE)

        if(-1==res$status)
        {
            stop('Unexpected error')
        }
        else if(0==res$status)
        {
            dim(res$node.pos.counts) <- c(longest, n)
            res$node.pos.counts <- t(res$node.pos.counts)
            rownames(res$node.pos.counts) <- unname(NP(community, 'node'))
            return (list(chain.lengths=res$chain.lengths, 
                         node.pos.counts=res$node.pos.counts))
                         
        }
        else if(1==res$status)
        {
            stop('Problem with an input parameter')
        }
        else
        {
            stop(paste('Unknown status [', res$status, ']', sep=''))
        }
    }
    else
    {
        return (NULL)
    }
}

TrophicChains <- function(community, node.properties=NULL, 
                          chain.properties=NULL)
{
    # Returns a data.frame
    # One row per unique trophic chain in the food web. 
    # Shorter chains are suffixed with ''. 
    # Cannibalism and loops are ignored, i.e. each node appears no more 
    # than once in a chain.
    # Isolated nodes are excluded so all chains are of length 2 or greater.

    max.queue <- .maxQueueOption()
    stopifnot(max.queue>=0)

    trace <- cat
    trace <- function(...) {}

    start <- Sys.time()

    if(!is.Community(community)) stop('Not a Community')

    # Small performance increase by removing cannibalistic links
    community <- RemoveCannibalisticLinks(community)
    if(0==NumberOfTrophicLinks(community))
    {
        return (NULL)
    }

    # Delegate to the C++ implementation.
    #    On iguana, ubuntu 10.04, R 2.14.1 (2011-12-22)
    #    library(cheddar)
    #    data(Benguela, BroadstoneStream, TL84, TL86, SkipwithPond,YthanEstuary)
    #    for(community in list(Benguela, BroadstoneStream, TL84, TL86, 
    #                          YthanEstuary))
    #    {
    #        print(community)
    #        print(system.time(chains <- TrophicChains(community)))
    #    }
    #   print(SkipwithPond)
    #   print(system.time(chains <- TrophicChains(SkipwithPond,max.chains=3e6)))

    # C++ implementation with dynamic queue
    #    Benguela containing 29 nodes.
    #       user  system elapsed 
    #      0.108   0.036   0.147 
    #    Broadstone Stream containing 37 nodes.
    #       user  system elapsed 
    #      0.088   0.048   0.136 
    #    Tuesday Lake sampled in 1984 containing 56 nodes.
    #       user  system elapsed 
    #      0.132   0.068   0.199 
    #    Tuesday Lake sampled in 1986 containing 57 nodes.
    #       user  system elapsed 
    #      0.124   0.076   0.197 
    #    Ythan Estuary containing 92 nodes.
    #       user  system elapsed 
    #      0.221   0.100   0.320 
    #    Skipwith Pond containing 37 nodes.
    #       user  system elapsed 
    #     12.608   3.436  16.345 

    # Old C implementation with fixed-size queue
    #    Tuesday Lake sampled in 1984 containing 56 nodes.
    #       user  system elapsed 
    #      0.092   0.024   0.115 
    #    Tuesday Lake sampled in 1986 containing 57 nodes.
    #       user  system elapsed 
    #      0.084   0.024   0.109 
    #    Ythan Estuary containing 92 nodes.
    #       user  system elapsed 
    #      0.112   0.052   0.168 
    #    Benguela containing 29 nodes.
    #       user  system elapsed 
    #      0.064   0.024   0.086 

    # Get the food web as an adjacency list
    alist <- .CAdjacencyList(community, ConsumersByNode(community))

    # 1 if basal, 0 otherwise
    is.basal <- as.integer(IsBasalNode(community))

    # Outputs
    # Compute the number of chains and the longest chain
    chains_size <- .TrophicChainsSize(community)

    trace('After .TrophicChainsSize', as.double(Sys.time()-start), '\n')

    ncols <- chains_size['n.chains']
    nrows <- chains_size['longest']

    # chains will be filled by the C function
    chains <- as.integer(rep(NA, ncols*nrows))
    status <- as.integer(-1)

    trace('Before trophic_chains', as.double(Sys.time()-start), '\n')

    res <- .C('trophic_chains', 
              as.integer(alist), 
              as.integer(length(alist)), 
              as.integer(is.basal), 
              as.integer(nrow(alist)), 
              as.integer(ncols), 
              as.integer(nrows), 
              as.integer(max.queue),
              chains=chains, 
              status=status, 
              PACKAGE='cheddar', NAOK=TRUE)

    trace('After trophic_chains', as.double(Sys.time()-start), '\n')

    if(-1==res$status)
    {
        stop('Unexpected error')
    }
    else if(0==res$status)
    {
        trace('Before dim', as.double(Sys.time()-start), '\n')
        dim(res$chains) <- c(nrows, ncols)

        # 1-indexed
        trace('Before 1+', as.double(Sys.time()-start), '\n')
        res$chains <- 1+res$chains

        # Prefer tall and narrow
        trace('Before t(chains)', as.double(Sys.time()-start), '\n')
        res$chains <- t(res$chains)

        trace('Before set names', as.double(Sys.time()-start), '\n')
        res$chains <- unname(NP(community, 'node'))[res$chains]
        trace('Before clear NA names', as.double(Sys.time()-start), '\n')
        res$chains[is.na(res$chains)] <- ''
        trace('Before dim', as.double(Sys.time()-start), '\n')
        dim(res$chains) <- c(ncols, nrows)

        trace('Before data.frame', as.double(Sys.time()-start), '\n')
        colnames(res$chains) <- paste('Node', 1:ncol(res$chains), sep='.')
        res$chains <- as.data.frame(res$chains, stringsAsFactors=FALSE)

        trace('Before .Chains', as.double(Sys.time()-start), '\n')
        res <- .ChainsWithProperties(res$chains, community, node.properties, 
                                     chain.properties)
        trace('After .Chains', as.double(Sys.time()-start), '\n')
        return (res)
    }
    else if(1==res$status)
    {
        stop('Problem with an input parameter')
    }
    else if(2==res$status)
    {
        stop(paste('Limit of', ncols, 'chains exceeded'))
    }
    else
    {
        stop(paste('Unknown status [', res$status, ']', sep=''))
    }
}

ThreeNodeChains <- function(community, exclude.loops=FALSE, 
                            node.properties=NULL, chain.properties=NULL)
{
    # Returns a matrix of 3 cols and a row per tri-trophic chain
    # Each node appears no more than once in a chain except for the 
    # case where a chain goes R > C > R, i.e. bottom==top. 
    # Setting exclude.loops=TRUE excludes these links
    # Cannibalism is ignored
    if(!is.Community(community)) stop('Not a Community')
    if(0==NumberOfTrophicLinks(community))
    {
        return (NULL)
    }

    # Get consumers of each node as indices
    consumers <- ConsumersByNode(community)
    chains <- NULL
    for(bottom in NP(community, 'node'))
    {
        for(intermediate in consumers[[bottom]])
        {
            top <- consumers[[intermediate]]
            if(length(top)>0)
            {
                new.chains <- cbind(bottom, intermediate, top)

                # Remove rows with duplicated node
                dup <- new.chains[,1] == new.chains[,2] | 
                       new.chains[,2] == new.chains[,3]
                if(any(dup))
                {
                    new.chains <- new.chains[!dup,,drop=FALSE]
                }

                if(nrow(new.chains)>0)
                {
                    chains <- rbind(chains, new.chains)
                }
            }
        }
    }

    if(is.null(chains))
    {
        return (NULL)
    }
    else
    {
        if(exclude.loops)
        {
            # Remove chains for which bottom==top
            # Chains with bottom==intermediate or intermediate==top were never 
            # added.
            chains <- chains[!apply(chains, 1, function(r) r[1]==r[3]),,
                             drop=FALSE]
        }

        chains <- as.data.frame(chains, stringsAsFactors=FALSE)
        return (.ChainsWithProperties(chains, community, node.properties, 
                                      chain.properties))
    }
}

.ChainsWithProperties <- function(chains, community, node.properties, 
                                  chain.properties)
{
    # Make a note of these columns before adding columns containing node 
    # properties
    chain.columns <- colnames(chains)

    if(!is.null(chain.properties))
    {
        cp <- lapply(chain.properties, function(p) 
        {
            do.call(p, args=list(community=community, chains=chains))
        })

        names(cp) <- chain.properties
        chains <- cbind(chains, do.call('cbind', cp))
    }

    chains <- .AddNodePropertiesToChains(community, chains=chains, 
                                         node.properties=node.properties, 
                                         chain.columns=chain.columns)

    return (chains)
}

.AddNodePropertiesToChains <- function(community, chains, node.properties, 
                                       chain.columns=colnames(chains))
{
    # chains - a data.frame
    # node.properties - names of functions that will be passed community and 
    # chain.columns - 
    if(!is.null(node.properties))
    {
        stopifnot(!is.null(colnames(chains)))
        stopifnot(all(chain.columns %in% colnames(chains)))

        np <- NPS(community, node.properties)

        for(column in chain.columns)
        {
            v <- lapply(colnames(np), 
                        function(p) np[chains[,column], p, drop=FALSE])
            v <- data.frame(v, stringsAsFactors=FALSE)
            colnames(v) <- paste(column, colnames(v), sep='.')
            chains <- cbind(chains, v)
        }
    }
    return (chains)
}
