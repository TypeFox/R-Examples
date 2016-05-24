# Community functions that operate on trophic links
.ResolveTrophicLinksToNodeNames <- function(community, links)
{
    # Returns a matrix with the columns 'resource' and 'consumer', both 
    # containing node names. 
    # links should be a matrix or data.frame with columns 'resource' and 
    # 'consumer'. If either column in links contains characters, they should be 
    # node names. If not characters, they should be integer node indices.
    # The function does not check to make sure that the links are present in 
    # TLPS(community)
    
    if(!is.Community(community)) stop('Not a Community')

    # Get appropriate columns from links matrix
    if(!is.null(colnames(links)))
    {
        stopifnot(all(is.element(c('resource','consumer'), colnames(links))))
        links <- links[,c('resource','consumer'),drop=FALSE]
    }
    else
    {
        # Assume first two columns
        stopifnot(ncol(links)>=2)
        links <- links[,c(1,2),drop=FALSE]
    }

    node <- unname(NP(community, 'node'))

    # Get resource and consumer as a node name
    f <- function(values)
    {
        if(is.factor(values))
        {
            return (as.character(values))
        }
        else if(!is.character(values))
        {
            stopifnot(all(values>0 & values<=length(node)))
            return (node[values])
        }
        else
        {
            return (values)
        }
    }

    links <- cbind(resource=f(links[,1]), consumer=f(links[,2]))
    rownames(links) <- NULL

    bad <- !links %in% node
    if(any(bad))
    {
        stop(paste('The names [', paste(unique(links[bad]), collapse=','), 
                   '] are not nodes', sep=''))
    }

    return (links)
}

.ResolveTrophicLinksToRowIndices <- function(community, links)
{
    # Returns a vector of integer indices in TLPS(community) for each row in 
    # links; NA for links that do not exist.

    # links should be a matrix or data.frame with columns 'resource' and 
    # 'consumer'. If either column in links contains characters, they should be 
    # node names. If not characters, they should be integer node indices.

    if(!is.Community(community)) stop('Not a Community')

    links <- .ResolveTrophicLinksToNodeNames(community, links)

    tlp <- TLPS(community)
    match <- apply(links, 1, function (row)
    {
        l <- which(row[1]==tlp[,1] & row[2]==tlp[,2])

        if(1==length(l))
        {
            return (l)
        }
        else if(0==length(l))
        {
            return (NA)
        }
        else
        {
            # If more than one row matches this link, something has gone really 
            # wrong.
            stop(paste('Link from [', row[1], '] to [', row[2], '] appears ', 
                       length(l), ' times in trophic links', sep=''))
        }
    })

    return (match)
}

PredationMatrix <- function(community, weight=NULL)
{
    # Returns a predation matrix. Columns and rows are named by node. 
    # Columns are consumers, rows are resource. If weight is NULL, elements 
    # will be either 0 (no trophic link) or 1 (trophic link).
    # If not NULL, weight should be the names of a link property; non-zero 
    # elements of the matrix will contain the value of the weight property.
    if(!is.Community(community)) stop('Not a Community')

    tlp <- TLPS(community, link.properties=weight)
    np <- NPS(community)
    n.nodes <- nrow(np)

    # Create predation matrix
    # Casts to as.integer keep pm as a matrix of integer rather than 
    # numeric, which I think what we want.
    pm <- matrix(as.integer(0), ncol=n.nodes, nrow=n.nodes)
    colnames(pm) <- rownames(pm) <- np$node

    if(!is.null(tlp))
    {
        stopifnot(all(tlp$resource %in% np$node))
        stopifnot(all(tlp$consumer %in% np$node))

        # Get numerical indices of resources and consumers
        res <- sapply(tlp$resource, function(n) return (which(n==np$node)))
        con <- sapply(tlp$consumer, function(n) return (which(n==np$node)))

        pm[res + n.nodes * (con-1)] <- as.integer(1)

        if(!is.null(weight))
        {
            pm[res + n.nodes * (con-1)] <- tlp[,weight]
        }
    }

    return (pm)
}

InDegree <- TrophicGenerality <- NumberOfResources <- function(community)
{
    # The number of resources of each node
    tlps <- TLPS(community)
    if(is.null(tlps))
    {
        res <- rep(0, NumberOfNodes(community))
        names(res) <- unname(NP(community, 'node'))
        return (res)
    }
    else
    {
        return (with(tlps, sapply(NP(community, 'node'), 
                                  function(n) sum(consumer==n))))
    }
}

NormalisedTrophicGenerality <- function(community)
{
    # The number of resources of each node, normalised with 
    # L/S. Williams and Martinez 2000 Nature.
    return ( (1/LinkageDensity(community)) * TrophicGenerality(community))
}

OutDegree <- TrophicVulnerability <- NumberOfConsumers <- function(community)
{
    # The number of consumers of each node
    tlps <- TLPS(community)
    if(is.null(tlps))
    {
        res <- rep(0, NumberOfNodes(community))
        names(res) <- unname(NP(community, 'node'))
        return (res)
    }
    else
    {
        return (with(tlps, sapply(NP(community, 'node'), 
                                  function(n) sum(resource==n))))
    }
}

NormalisedTrophicVulnerability <- function(community)
{
    # The number of consumers of each node, normalised with 
    # L/S. Williams and Martinez 2000 Nature.
    return ( (1/LinkageDensity(community)) * TrophicVulnerability(community))
}

IsBasalNode <- function(community)
{
    # TRUE for those nodes that consume no others, excluding 
    # themselves; will therefore return TRUE for nodes that are cannibalistic, 
    # have no resources and some consumers. Biologically, this should not 
    # happen, but mathematically, it should be accounted for.
    community <- RemoveCannibalisticLinks(community)
    return (0==NumberOfResources(community) & 0<NumberOfConsumers(community))
}

BasalNodes <- function(community)
{
    # Names of nodes that consume no others
    return (names(which(IsBasalNode(community))))
}

IsNonBasalNode <- function(community)
{
    # TRUE for those nodes that are not basal
    return (!IsBasalNode(community))
}

NonBasalNodes <- function(community)
{
    # Names of nodes that consume at least one other
    return (names(which(IsNonBasalNode(community))))
}

IsTopLevelNode <- function(community)
{
    # TRUE for those node that are consumed by no others, excluding 
    # themselves; will therefore return TRUE for cannibalistic top level 
    # predators.
    community <- RemoveCannibalisticLinks(community)
    return (0==NumberOfConsumers(community) & 0<NumberOfResources(community))
}

TopLevelNodes <- function(community)
{
    # Names of nodes that are consumed by no others, excluding themselves.
    return (names(which(IsTopLevelNode(community))))
}

IsNonTopLevelNode <- function(community)
{
    # Cannibalism ignored
    return (!IsTopLevelNode(community))
}

NonTopLevelNodes <- function(community)
{
    # Cannibalism ignored
    return (names(which(IsNonTopLevelNode(community))))
}

IsIntermediateNode <- function(community)
{
    # Neither basal nor top level. Cannibalism ignored.
    community <- RemoveCannibalisticLinks(community)
    return (0<NumberOfResources(community) & 0<NumberOfConsumers(community))
}

IntermediateNodes <- function(community)
{
    # Neither basal nor top level
    return (names(which(IsIntermediateNode(community))))
}

IsIsolatedNode <- function(community)
{
    # TRUE for those node with no resource and no consumers
    # Cannibalism ignored
    community <- RemoveCannibalisticLinks(community)
    return (0==NumberOfResources(community) & 0==NumberOfConsumers(community))
}

IsolatedNodes <- function(community)
{
    # Nodes with no resources and no consumers
    return (names(which(IsIsolatedNode(community))))
}

IsConnectedNode <- function(community)
{
    # Cannibalism ignored
    community <- RemoveCannibalisticLinks(community)
    return (!IsIsolatedNode(community))
}

ConnectedNodes <- function(community)
{
    return (names(which(IsConnectedNode(community))))
}

ResourcesByNode <- function(community)
{
    # A named list of length NumberOfNodes(). Elements are vectors containing 
    # those nodes that are resources.
    tlps <- TLPS(community)
    if(is.null(tlps))
    {
        return (lapply(NP(community, 'node'), 
                       function(n) vector(mode='character')))
    }
    else
    {
        return (with(tlps, lapply(NP(community, 'node'), 
                                  function(n) resource[consumer==n])))
    }
}

ResourcesOfNodes <- function(community, nodes)
{
    # Node should be either a vector of names or a vector of integer indices
    if(!is.Community(community)) stop('Not a Community')

    nodes <- .ResolveToNodeIndices(community, nodes)
    r <- ResourcesByNode(community)
    if(1==length(nodes))
    {
        # A vector of resources
        return (r[[nodes]])
    }
    else
    {
        # A list of resource by node
        return (r[nodes])
    }
}

ConsumersByNode <- function(community)
{
    # A named list of length NumberOfNodes(). Elements are vectors containing 
    # those nodes that are consumers.
    if(!is.Community(community)) stop('Not a Community')

    tlps <- TLPS(community)
    if(is.null(tlps))
    {
        return (lapply(NP(community, 'node'), 
                       function(n) vector(mode='character')))
    }
    else
    {
        return (with(tlps, lapply(NP(community, 'node'), 
                                  function(n) consumer[resource==n])))
    }
}

ConsumersOfNodes <- function(community, nodes)
{
    if(!is.Community(community)) stop('Not a Community')

    nodes <- .ResolveToNodeIndices(community, nodes)
    c <- ConsumersByNode(community)
    if(1==length(nodes))
    {
        # A vector of consumers
        return (c[[nodes]])
    }
    else
    {
        # A list of consumers by node
        return (c[nodes])
    }
}

ResourcesAndConsumersByNode <- function(community)
{
    # A named list of length NumberOfNodes(). Elements are vectors containing 
    # those nodes that are either resources or consumers.
    rc <- mapply(c, ResourcesByNode(community), ConsumersByNode(community), 
                 SIMPLIFY=FALSE)
    return (lapply(rc, unique))
}

NumberOfTrophicLinks <- function(community)
{
    # The total number of links in the web.
    tlps <- TLPS(community)
    if(is.null(tlps))   return (0)
    else                return (nrow(tlps))
}

LinkageDensity <- function(community)
{
    # The number of links per node
    return (NumberOfTrophicLinks(community) / NumberOfNodes(community))
}

Degree <- function(community)
{
    # A vector of length NumberOfNodes() containing the number of links
    # for each node. Food web links are directed so cannibalistic links 
    # count twice.
    return (InDegree(community) + OutDegree(community))
}

DegreeDistribution <- function(community, cumulative=FALSE)
{
    # A vector of proportions of nodes with 0:max(degree) trophic 
    # links. 
    degree <- Degree(community)
    dd <- tabulate(factor(degree, levels=0:max(degree)))
    names(dd) <- 0:max(degree)
    stopifnot(sum(dd)==NumberOfNodes(community))
    dd <- dd / sum(dd)

    if(cumulative)
    {
        dd <- rev(cumsum(rev(dd)))
    }

    return (dd)
}

DirectedConnectance <- function(community)
{
    # Martinez 1992 Ecological Monographs
    return (NumberOfTrophicLinks(community) / NumberOfNodes(community)^2)
}

IsCannibal <- function(community)
{
    # TRUE if the node is a consumer of itself
    tlps <- TLPS(community)
    if(is.null(tlps))
    {
        res <- rep(FALSE, NumberOfNodes(community))
        names(res) <- unname(NP(community, 'node'))
        return (res)
    }
    else
    {
        return (with(tlps, sapply(NP(community, 'node'), 
                                function(n) 1==sum(resource==n & consumer==n))))
    }
}

Cannibals <- function(community)
{
    # The names of those nodes which consume themselves
    return (names(which(IsCannibal(community))))
}

FractionCannibalistic <- function(community)
{
    n <- IsCannibal(community)
    return (sum(n)/length(n))
}

.MatrixInversionTL <- function(community, weight.by, include.isolated)
{
    # Returns a vector of trophic levels. If the trophic levels can not be 
    # computed then there is a problem with the community's network topology 
    # (see comment from Rich below), in which case all elements of the 
    # returned vector will be NA. 

    # Use the method described in Williams and Martinez (2004) Am Nat, Limits 
    # to Trophic Levels and Omnivory in Complex Food Webs: Theory and Data. 
    # Algorithm proposed by Levine 1980 J Theor Biol, Several measures of 
    # trophic structure applicable to complex food webs

    # weight.by should be either NULL or the name of a trophic link property 
    # or function that gives diet fraction or similar. 
    # If weight.by is NULL, the returned values are 'prey-averaged TL' 
    # (Williams and Martinez (2004) Am Nat. If weight.by is not NULL, 
    # returned values are 'flow-based TL'.

    # This method is quick and easy and accounts for flows through loops. 

    if(0==NumberOfTrophicLinks(community))
    {
        tl <- rep(NA, NumberOfNodes(community))
        names(tl) <- unname(NP(community, 'node'))
    }
    else
    {
        # If requested by the caller, weighting is added to the matrix by 
        # PredationMatrix().
        # In the absence of weighting / diet fraction information weight 
        # each link will be equally weighted.
        pm <- PredationMatrix(community, weight=weight.by)

        # The algorithm needs consumers on rows
        W <- t(pm)

        # Make each row sum to one.
        rs <- rowSums(W)
        W <- W / matrix(rs, ncol=ncol(W), nrow=nrow(W))
        W[0==rs,] <- 0      # Fix NA resulting from div by zero

        # Identity matrix
        I <- diag(ncol(W))

        # Invert the matrix. From Rich 2011-08-17: "If this matrix inversion 
        # fails, there is an important problem with the network topology. 
        # For a food web to be energetically feasible, every node must be 
        # connected to a basal node. When the inversion fails it is because 
        # there is at least one node that has no connection to a basal node."
        result <- tryCatch(solve(I-W), error = function(e) e)

        if('error' %in% class(result))
        {
            tl <- rep(NA, ncol(pm))
            names(tl) <- colnames(pm)
        }
        else
        {
            # Returned vector is named by node
            tl <- rowSums(result)   # Isolated nodes have tl of 1

            if(!include.isolated)
            {
                tl[IsolatedNodes(community)] <- NA
            }
        }
    }

    return (tl)
}

PreyAveragedTrophicLevel <- function(community, include.isolated=TRUE)
{

    # Williams and Martinez (2004) Am Nat, p 460
    # "Prey-averaged TL is equal to 1 + the mean TL of all the consumer's 
    #  trophic resources TLj =  1 + sum{i=1 to s} l_{ij} TL_i / n_j
    #  where nj is the number of prey species in the diet of species j. 
    #  This equation is equivalent to equation (1) [flow-based trophic level], 
    #  with each nonzero link strength pij p 1/nj, which assumes that a consumer 
    #  consumes all its prey species equally. "
    return (.MatrixInversionTL(community, weight.by=NULL, 
                               include.isolated=include.isolated))
}

FlowBasedTrophicLevel <- function(community, weight.by, include.isolated=TRUE)
{
    # Williams and Martinez (2004) Am Nat, p 460
    # Flow-based trophic level
    stopifnot(!missing(weight.by))
    return (.MatrixInversionTL(community, weight.by=weight.by, 
                               include.isolated=include.isolated))
}
TrophicLevels <- function(community, weight.by=NULL, include.isolated=TRUE)
{
    # Returns a matrix of either six or seven different measures of trophic 
    # level: the matrix-inversion method of PreyAveragedTrophicLevel, 
    # the matrix-inversion method of FlowBasedTrophicLevel (if weight is not 
    # NULL) and the five chain-position methods.
    patl <- PreyAveragedTrophicLevel(community, 
                                     include.isolated=include.isolated)
    fbtl <- NULL
    if(!is.null(weight.by))
    {
        fbtl <- FlowBasedTrophicLevel(community, weight.by=weight.by, 
                                      include.isolated=include.isolated)
    }

    # The following measures of trophic level are computed using the 
    # position of nodes in each unique food chain.
    chain.stats <- TrophicChainsStats(community)

    S <- NumberOfNodes(community)
    if(is.null(chain.stats))
    {
        # If no trophic links or if NumberOfTrophicLinks()>0 but chains is 
        # NULL (will happen if all links are cannibalistic), all trophic levels 
        # are NA
        stl <- swtl <- ltl <- lwtl <- catl <- rep(NA, S)
    }
    else
    {
        # Williams and Martinez (2004) Am Nat, p 460
        # "Shortest TL is equal to 1 + the shortest chain length from a 
        #  consumer to a basal species."
        stl <- apply(chain.stats$node.pos.counts, 1, function(row)
        {
            if(all(0==row)) NA
            else            min(which(0!=row))
        })

        # Williams and Martinez (2004) Am Nat, p 460
        # "Short-weighted TL is the average of shortest TL and prey-averaged 
        #  TL. This gives a measure biased toward shorter food chains."
        swtl <- sapply(1:S, function(n) mean(c(stl[n],patl[n])))

        # Williams and Martinez (2004) Am Nat, p 460
        # "Longest TL is equal to 1 + the longest chain length from a consumer 
        #  to a basal species."
        ltl <- apply(chain.stats$node.pos.counts, 1, function(row)
        {
            if(all(0==row)) NA
            else            max(which(0!=row))
        })

        # Williams and Martinez (2004) Am Nat, p 460
        # "Long-weighted TL is the average of longest TL and prey-averaged TL. 
        #  This gives a measure biased toward longer food chains. "
        lwtl <- sapply(1:S, function(n) mean(c(ltl[n], patl[n])))

        # I think that Williams and Martinez's Chain-averaged TL is the same 
        # as Jonsson et al's Trophic Height:

        # Williams and Martinez (2004) Am Nat, p 460
        # "Chain-averaged TL is equal to 1 + the average chain length of all 
        #  paths from a species to a basal species (Martinez 1991; Polis 1991; 
        #  Fussman and Heber 2002)."

        # Jonsson et al (2005) AER, p 8
        # "The trophic position of a species in a food chain is 1 + the 
        #  number of species preceding it in the ordered list of species in the 
        #  chain... Trophic height is the average trophic position of a species 
        #  in all food chains of which it is a part."

        # A additional comment on flow-based vs chain-averaged TL
        # Williams and Martinez (2004) Am Nat, p 460
        # "Flow-based TL and prey-averaged TL are both computed using the 
        #  matrix algebra method of Levine (1980) based on summing an infinite 
        #  geometric series that includes the contributions from all loops. 
        #  In contrast, the computation of chain-averaged TL maintains 
        #  tractability in complex food webs by only passing through a loop 
        #  once (Martinez 1991)."
        catl <- apply(chain.stats$node.pos.counts, 1, function(row)
        {
            if(all(0==row)) NA
            else            weighted.mean(1:length(row), row)
        })

        if(include.isolated)
        {
            # Isolated species will have been assigned a TL of NA in the 
            # following measures. Set these to 1.
            i <- IsIsolatedNode(community)
            stl[i] <- swtl[i] <- ltl[i] <- lwtl[i] <- catl[i] <- 1
        }
    }

    return (cbind(ShortestTL=stl, 
                  ShortWeightedTL=swtl, 
                  LongestTL=ltl, 
                  LongWeightedTL=lwtl, 
                  ChainAveragedTL=catl, 
                  PreyAveragedTL=patl, 
                  FlowBasedTL=fbtl))
}

# Convenience functions to access just one of the chain-averaged measures
ShortestTrophicLevel <- function(community, include.isolated=TRUE)
{
    tl <- TrophicLevels(community, weight.by=NULL, 
                        include.isolated=include.isolated)
    res <- tl[,'ShortestTL']
    names(res) <- rownames(tl)
    return (res)
}

ShortWeightedTrophicLevel <- function(community, include.isolated=TRUE)
{
    tl <- TrophicLevels(community, weight.by=NULL, 
                        include.isolated=include.isolated)
    res <- tl[,'ShortWeightedTL']
    names(res) <- rownames(tl)
    return (res)
}

LongestTrophicLevel <- function(community, include.isolated=TRUE)
{
    tl <- TrophicLevels(community, weight.by=NULL, 
                        include.isolated=include.isolated)
    res <- tl[,'LongestTL']
    names(res) <- rownames(tl)
    return (res)
}

LongWeightedTrophicLevel <- function(community, include.isolated=TRUE)
{
    tl <- TrophicLevels(community, weight.by=NULL, 
                        include.isolated=include.isolated)
    res <- tl[,'LongWeightedTL']
    names(res) <- rownames(tl)
    return (res)
}

ChainAveragedTrophicLevel <- 
TrophicHeight <- function(community, include.isolated=TRUE)
{
    tl <- TrophicLevels(community, weight.by=NULL, 
                        include.isolated=include.isolated)
    res <- tl[,'ChainAveragedTL']
    names(res) <- rownames(tl)
    return (res)
}

TrophicSpecies <- function(community, include.isolated=TRUE)
{
    # Returns a vector containing the trophic species number of each 
    # node in the community. Nodes with identical sets of resources and 
    # consumers are assigned to the same trophic species.

    # If include.isolated is TRUE, isolated node are assigned a trophic 
    # node. If include.isolated is FALSE isolated nodes have a trophic 
    # node of NA.

    #    library(cheddar)
    #    data(Benguela, BroadstoneStream, TL84, TL86, SkipwithPond,YthanEstuary)
    #    for(community in list(Benguela, BroadstoneStream, TL84, TL86, 
    #                          SkipwithPond, YthanEstuary))
    #    {
    #        print(community)
    #        print(system.time(ts <- TrophicSpecies(community)))
    #    }

    #    Benguela containing 29 nodes.
    #       user  system elapsed 
    #      0.008   0.000   0.009 
    #    Broadstone Stream containing 37 nodes.
    #       user  system elapsed 
    #      0.008   0.000   0.008 
    #    Tuesday Lake sampled in 1984 containing 56 nodes.
    #       user  system elapsed 
    #      0.012   0.000   0.012 
    #    Tuesday Lake sampled in 1986 containing 57 nodes.
    #       user  system elapsed 
    #      0.012   0.000   0.012 
    #    Skipwith Pond containing 37 nodes.
    #       user  system elapsed 
    #      0.012   0.000   0.012 
    #    Ythan Estuary containing 92 nodes.
    #       user  system elapsed 
    #      0.032   0.000   0.030 

    if(!is.Community(community)) stop('Not a Community')
    trophic.species = rep(NA, NumberOfNodes(community))
    names(trophic.species) <- unname(NP(community, 'node'))

    remaining <- 1:NumberOfNodes(community)

    if(!include.isolated)
    {
        remaining <- setdiff(remaining, which(IsIsolatedNode(community)))
    }

    resources <- ResourcesByNode(community)
    consumers <- ConsumersByNode(community)

    trophic.species.number <- 1
    while(length(remaining)>0)
    {
        # Consider the first node that has not yet been assigned a trophic 
        # species
        node <- remaining[1]

        # Remove from remaining and assign trophic species number
        remaining <- setdiff(remaining, node)

        trophic.species[node] <- trophic.species.number

        for(t in remaining)
        {
            if(length(resources[[node]])==length(resources[[t]]) && 
               length(consumers[[node]])==length(consumers[[t]]) && 
               all(resources[[node]]==resources[[t]]) && 
               all(consumers[[node]]==consumers[[t]]))
            {
                # Remove from remaining and assign trophic species number
                remaining <- setdiff(remaining, t)
                trophic.species[t] <- trophic.species.number
            }
        }

        trophic.species.number <- trophic.species.number + 1
    }

    return (trophic.species)
}

TrophicLinksForNodes <- function(community, nodes, node.properties=NULL, 
                                 link.properties=NULL)
{
    # Returns a matrix of links in-to and out-of the given nodes. 
    # nodes should be either node numbers or node names. 
    if(!is.Community(community)) stop('Not a Community')

    stopifnot((is.character(nodes) && all(nodes %in% NP(community, 'node'))) || 
              all(0<nodes & nodes<=NumberOfNodes(community)))
    if(!is.character(nodes))
    {
        nodes <- unname(NP(community, 'node')[nodes])
    }

    tlp <- TLPS(community, node.properties=node.properties, 
                link.properties=link.properties)
    links <- tlp[tlp$resource %in% nodes | tlp$consumer %in% nodes,,drop=FALSE]
    rownames(links) <- NULL
    return (links)
}

.DoLinksExist <- function(community, links)
{
    # Returns TRUE or FALSE for each row in links.
    # links should be either a matrix or a data.frame. 
    # links should contain column called resource and consumer.
    # Values should be either node numbers or node names.
    l <- .ResolveTrophicLinksToRowIndices(community, links)
    return (!is.na(l))
}

RemoveCannibalisticLinks <- function(community, 
                                     title=paste(CP(community, 'title'), 
                                               '(cannibalistic links removed)'))
{
    # Returns a copy of community with cannibalistic links removed
    if(!is.Community(community)) stop('Not a Community')

    tlp <- TLPS(community)
    if(is.null(tlp))
    {
        return (community)
    }
    else
    {
        remove <- tlp[,1] == tlp[,2]
        tlp <- tlp[!remove,,drop=FALSE]

        properties <- CPS(community)
        properties$title <- title

        return (Community(NPS(community), trophic.links=tlp, 
                          properties=properties))
    }
}

FractionBasalNodes <- function(community)
{
    n <- IsBasalNode(community)
    return (length(which(n)) / length(n))
}

FractionNonBasalNodes <- function(community)
{
    n <- IsNonBasalNode(community)
    return (length(which(n)) / length(n))
}

FractionTopLevelNodes <- function(community)
{
    n <- IsTopLevelNode(community)
    return (length(which(n)) / length(n))
}

FractionNonTopLevelNodes <- function(community)
{
    n <- IsNonTopLevelNode(community)
    return (length(which(n)) / length(n))
}

FractionIntermediateNodes <- function(community)
{
    n <- IsIntermediateNode(community)
    return (length(which(n)) / length(n))
}

FractionIsolatedNodes <- function(community)
{
    n <- IsIsolatedNode(community)
    return (length(which(n)) / length(n))
}

FractionConnectedNodes <- function(community)
{
    n <- IsConnectedNode(community)
    return (length(which(n)) / length(n))
}

ShortestPaths <- function(community, weight.by=NULL)
{
    # Returns a matrix of NumberOfNodes square. 
    # Dijkstra's algorithm is used to compute path lengths. 

    #    library(cheddar)
    #    data(Benguela, BroadstoneStream, TL84, TL86, SkipwithPond,YthanEstuary)
    #    for(community in list(Benguela, BroadstoneStream, TL84, TL86, 
    #                          SkipwithPond, YthanEstuary))
    #    {
    #        print(community)
    #        print(system.time(y <- ShortestPaths(community)))
    #    }
    #    Benguela containing 29 nodes.
    #       user  system elapsed 
    #      0.024   0.000   0.024 
    #    Broadstone Stream containing 37 nodes.
    #       user  system elapsed 
    #      0.020   0.000   0.018 
    #    Tuesday Lake sampled in 1984 containing 56 nodes.
    #       user  system elapsed 
    #      0.032   0.000   0.034 
    #    Tuesday Lake sampled in 1986 containing 57 nodes.
    #       user  system elapsed 
    #      0.032   0.000   0.032 
    #    Skipwith Pond containing 37 nodes.
    #       user  system elapsed 
    #      0.032   0.000   0.034 
    #    Ythan Estuary containing 92 nodes.
    #       user  system elapsed 
    #      0.060   0.000   0.067 

    if(!is.Community(community)) stop('Not a Community')

    # Get the food web as two adjancency lists
    consumers <- .CAdjacencyList(community, ConsumersByNode(community))
    resources <- .CAdjacencyList(community, ResourcesByNode(community))
    pm <- PredationMatrix(community, weight=weight.by)

    n <- NumberOfNodes(community)
    lengths <- as.double(rep(NA, n*n))
    status <- as.integer(-1)
    res <- .C('shortest_paths', 
              as.integer(consumers), 
              as.integer(length(consumers)), 
              as.integer(resources), 
              as.integer(length(resources)), 
              as.double(pm), 
              as.integer(n), 
              lengths=lengths, 
              status=status, 
              PACKAGE='cheddar', 
              NAOK=TRUE)

    if(0==res$status)
    {
        # Take the subset of the matrix that we are interested in
        sp <- res$lengths
        dim(sp) <- c(n,n)
        colnames(sp) <- rownames(sp) <- unname(NP(community, 'node'))
    }
    else if(1==res$status)
    {
        stop('Problem with an input parameter')
    }
    else
    {
        stop(paste('Unknown status [', res$status, ']', sep=''))
    }

    return (sp)
}

CharacteristicPathLength <- function(community)
{
    # From Williams and Martinez 2002 PNAS:
    # "...characteristic path length, D quantifies the average number of
    # links necessary for information or an effect to propagate along
    # the shortest paths between nodes in a network"
    spl <- ShortestPaths(community)
    spl[is.infinite(spl)] <- NA
    return (mean(spl, na.rm=TRUE))
}

.SumPredationMatrixGap <- function(community, diet.gap)
{
    # Stouffer et al (2006) PNAS. Delegate to the C implementation not for 
    # speed but for consistency. 
    if(!is.Community(community)) stop('Not a Community')
    pm <- PredationMatrix(community)
    if(!diet.gap)
    {
        pm <- t(pm)
    }

    n <- as.integer(ncol(pm))
    sum_diet_gaps <- as.integer(0)
    status <- as.integer(-1)

    res <- .C('sum_diet_gaps', 
              as.integer(pm), 
              as.integer(n), 
              as.integer((1:n) - 1),   # Row order, 0-indexed
              sum_diet_gaps=sum_diet_gaps, 
              status=status)
    if(-1==res$status)
    {
        stop('Unexpected error')
    }
    else if(0==res$status)
    {
        # C 0-index back to R 1-index
        return (res$sum_diet_gaps)
    }
    else if(1==res$status)
    {
        stop('Problem with an input parameter')
    }
    else(1==res$status)
    {
        stop(paste('Unknown status [', res$status, ']', sep=''))
    }
}

SumDietGaps <- function(community)
{
    # Stouffer et al (2006) PNAS. Delegate to the C implementation not for 
    # speed but for consistency. 
    return (.SumPredationMatrixGap(community, diet.gap=TRUE))
}

SumConsumerGaps <- function(community)
{
    # Zook et al (2011) J Theor Biol. Delegate to the C implementation not for 
    # speed but for consistency. 
    return (.SumPredationMatrixGap(community, diet.gap=FALSE))
}

.MinimisePredationMatrixGaps <- function(community, T.start, T.stop, c, 
                                         swaps.per.T, trace.anneal, 
                                         diet.gap, n, best)
{
    # Using Stouffer et al's (2006) PNAS simulated annealing learning method.
    # If diet.gap is TRUE then the sum diet gap is minimized otherwise the 
    # sum in each species' consumers is minimized (ZookEtAl2011JTheorEcol).
 
    #    library(cheddar)
    #    data(Benguela, BroadstoneStream, TL84, TL86, SkipwithPond,YthanEstuary)
    #    for(community in list(Benguela, BroadstoneStream, TL84, TL86, 
    #                          SkipwithPond,YthanEstuary))
    #    {
    #        print(community)
    #        print(system.time(i <- MinimiseSumDietGaps(community)))
    #        print(paste('Diet gap:',i$sum.gaps))
    #    }

    #    Benguela containing 29 nodes.
    #       user  system elapsed 
    #      0.144   0.000   0.144 
    #    [1] "Diet gap: 27"
    #    Broadstone Stream containing 37 nodes.
    #       user  system elapsed 
    #      0.080   0.004   0.083 
    #    [1] "Diet gap: 7"
    #    Tuesday Lake sampled in 1984 containing 56 nodes.
    #       user  system elapsed 
    #      0.204   0.000   0.204 
    #    [1] "Diet gap: 8"
    #    Tuesday Lake sampled in 1986 containing 57 nodes.
    #       user  system elapsed 
    #      0.184   0.004   0.185 
    #    [1] "Diet gap: 9"
    #    Skipwith Pond containing 37 nodes.
    #       user  system elapsed 
    #      0.156   0.000   0.159 
    #    [1] "Diet gap: 31"
    #    Ythan Estuary containing 92 nodes.
    #       user  system elapsed 
    #      0.512   0.000   0.513 
    #    [1] "Diet gap: 262"

    if(!is.Community(community)) stop('Not a Community')

    stopifnot(T.start>T.stop)
    stopifnot(T.stop>0)
    stopifnot(c>0 && c<1)
    stopifnot(swaps.per.T>0)
    if(is.null(TLPS(community)))
    {
        stop('The community has no trophic links')
    }
    stopifnot(n>0)

    pm <- PredationMatrix(community)
    if(!diet.gap)
    {
        pm <- t(pm)
    }

    F <- function()
    {
        # Out-params
        sum.gaps <- as.integer(-1)
        best <- as.integer(rep(NA, NumberOfNodes(community)))
        status <- as.integer(-1)

        res <- .C('minimise_sum_diet_gaps', 
                  as.integer(pm), 
                  as.integer(NumberOfNodes(community)), 
                  as.double(T.start), 
                  as.double(T.stop), 
                  as.double(c), 
                  as.integer(swaps.per.T), 
                  as.integer(trace.anneal), 
                  sum.gaps=sum.gaps, 
                  best=best, 
                  status=status, 
                  PACKAGE='cheddar', 
                  NAOK=TRUE)

        if(-1==res$status)
        {
            stop('Unexpected error')
        }
        else if(0==res$status)
        {
            # C 0-index back to R 1-index
            best <- res$best + 1
            stopifnot(all(best>0 & best<=NumberOfNodes(community)))
            
            reordered <- OrderCommunity(community, new.order=best, 
                                title=paste(CP(community, 'title'), 
                                            'sum', 
                                             ifelse(diet.gap, 'diet', 'consumer'), 
                                            'gap minimised to',
                                            res$sum.gaps))
            return (list(sum.gaps=res$sum.gaps, 
                         order=unname(NP(community, 'node'))[best], 
                         reordered=reordered))
        }
        else if(1==res$status)
        {
            stop('Problem with an input parameter')
        }
        else(1==res$status)
        {
            stop(paste('Unknown status [', res$status, ']', sep=''))
        }
    }

    res <- replicate(n, F(), simplify=FALSE)
    res <- res[order(sapply(res, '[[', 'sum.gaps'))]
    if(best || 1==n)
    {
        return (res[[1]])
    }
    else
    {
        return (res)
    }
}

MinimiseSumDietGaps <- function(community, T.start=10, T.stop=0.1, c=0.9, 
                                swaps.per.T=1000, trace.anneal=FALSE, n=1, 
                                best=TRUE)
{
    return (.MinimisePredationMatrixGaps(community, T.start, T.stop, c, 
                                         swaps.per.T, trace.anneal, 
                                         diet.gap=TRUE, n, best))

}

MinimiseSumConsumerGaps <- function(community, T.start=10, T.stop=0.1, c=0.9, 
                                    swaps.per.T=1000, trace.anneal=FALSE, 
                                    n=1, best=TRUE)
{
    return (.MinimisePredationMatrixGaps(community, T.start, T.stop, c, 
                                         swaps.per.T, trace.anneal, 
                                         diet.gap=FALSE, n, best))
}

TrophicSimilarity <- function(community)
{
    # Martinez 1991 Ecol. Monog. 61, 367-392 p 370.
    # I = c/(a+b+c)
    # where c = number of predators and prey common to the two taxa, 
    # a = number of predators and prey unique to one taxon, and b = number 
    # of predators and prey unique to one taxon. When the two taxa have the 
    # same set of predators and prey, I = 1, When the two taxa have no 
    # common predators or common prey, I = 0.

    rc <- ResourcesAndConsumersByNode(community)
    nodes <- names(rc)
    I <- matrix(NA, ncol=length(nodes), nrow=length(nodes), 
                  dimnames=list(nodes,nodes))
    diag(I) <- 1

    if(1<NumberOfNodes(community))
    {
        # Inner loop depends upon the community having more than one node
        for(x in 1:(length(nodes)-1))
        {
            for(y in (x+1):length(nodes))
            {
                a <- length(setdiff(rc[[x]], rc[[y]]))
                b <- length(setdiff(rc[[y]], rc[[x]])) 
                c <- length(intersect(rc[[x]], rc[[y]]))
                I[x,y] <- I[y,x] <- c/(a+b+c)
            }
        }
    }

    I[IsolatedNodes(community),IsolatedNodes(community)] <- NA
    stopifnot(all.equal(I, t(I)))
    stopifnot(all(is.na(I)) || 1>=max(I, na.rm=TRUE))
    stopifnot(1==diag(I)[-which(IsIsolatedNode(community))])
    return (I)
}

MeanMaximumTrophicSimilarity <- function(community)
{
    # Williams and Martinez 2000 Nature 404 180-182, p 181.
    # MxSim = 1/S \sum_1^S max(s_{ij}) (i!=j)
    I <- TrophicSimilarity(community)
    diag(I) <- NA    # Exclude diagonal
    return (mean(apply(I, 2, max, na.rm=TRUE)))
}

IsOmnivore <- function(community, level=PreyAveragedTrophicLevel)
{
    # Species that consume two or more species and have a non-integer 
    # trophic level
    # Polis, G. A. Complex desert food webs: an empirical critique of food web 
    # theory. Am. Nat. 138, 123-155 (1991).
    n.resources <- sapply(ResourcesByNode(community), length)
    tl <- level(community)
    return (n.resources>=2 & 0!=(tl %% 1))
}

Omnivores <- function(community, ...)
{
    return (names(which(IsOmnivore(community, ...))))
}

FractionOmnivorous <- Omnivory <- function(community, ...)
{
    # The fraction of species that are omnivores
    return (length(Omnivores(community, ...))/NumberOfNodes(community))
}

NodeQuantitativeDescriptors <- function(community, weight)
{
    # Bersier et al (2002) Ecology

    if(!is.Community(community)) stop('Not a Community')
    .RequireTrophicLinks(community)

    # A common operation
    vlog2v <- function(v)   v*log2(v)

    b <- PredationMatrix(community, weight)
    bIn <- colSums(b)
    bOut <- rowSums(b)

    # Diversity of inflows and outflows
    HN <- -rowSums(vlog2v(t(b)/bIn), na.rm=TRUE)    # p 2397, eq 5
    HP <- -rowSums(vlog2v(b/bOut), na.rm=TRUE)      # p 2397, eq 6

    # Equivalent numbers of prey and predators
    nN <-2^HN          # p 2397, eq 7
    nN[0==bIn] <- 0
    nP <-2^HP          # p 2397, eq 8
    nP[0==bOut] <- 0

    d.prime <- nN / (nN + nP)            # p 2397, eq 9
    d <- bIn*nN / (bIn*nN + bOut*nP)     # p 2397, eq 10

    # p 2396:
    # Trophic level of a taxon (Yodzis 1989); the characterization we adopt 
    # here is 'one plus the length of the longest chain from the focal taxon
    # to a basal taxon.'
    t <- LongestTrophicLevel(community)
    res <- ResourcesByNode(community)
    o.prime <- -1 + 2^sapply(1:NumberOfNodes(community), function(k)
    {
        # Compute the number of prey taxa located at each trophic level t
        r <- res[[k]]     # Resources of k
        n <- table(t[r])  # The number of resources at each trophic level
        return (-sum(vlog2v(n/sum(n))))
    })

    o <- -1 + 2^sapply(1:NumberOfNodes(community), function(k)
    {
        r <- res[[k]]     # Resources of k
        # Biomass going into k summed by trophic level of resource
        bt <- tapply(b[r,k], t[r], sum)
        return (-sum(vlog2v(bt/sum(bt))))
    })

    # Standardised quantitative G and V, p 2400, eq 28-29
    g.prime <- nN*NumberOfNodes(community) / sum(nN)
    g <- bIn*nN*NumberOfNodes(community) / sum(bIn*nN)
    v.prime <- nP*NumberOfNodes(community) / sum(nP)
    v <- bOut*nP*NumberOfNodes(community) / sum(bOut*nP)

    # Bersier et al table 1
    res <- cbind(NResources=NumberOfResources(community), 
                 NConsumers=NumberOfConsumers(community), 
                 bIn,
                 bOut,
                 nN, 
                 nP, 
                 d.prime, 
                 d,
                 o.prime, 
                 o,
                 g.prime, 
                 g, 
                 v.prime, 
                 v)
    rownames(res) <- unname(NP(community, 'node'))
    return (res)
}

QuantitativeDescriptors <- function(community, weight, top.level.threshold=0.99)
{
    # A common operation
    vlog2v <- function(v)   v*log2(v)

    # Community-level descriptors are derived from the node-level decsriptors

    # TODO prevent double-computation of chains

    np <- NodeQuantitativeDescriptors(community, weight)

    b <- PredationMatrix(community, weight)
    sumb <- sum(b)

    # Diversity of inflows and outflows
    HN <- -rowSums(vlog2v(t(b)/np[,'bIn']), na.rm=TRUE)    # p 2397, eq 5
    HP <- -rowSums(vlog2v(b/np[,'bOut']), na.rm=TRUE)      # p 2397, eq 6

    # See comment at bottom of 2397 about about top-level threshold
    fracT.q.prime <- mean(np[,'d.prime']>=top.level.threshold)
    fracI.q.prime <- mean(0<np[,'d.prime'] & np[,'d.prime']<top.level.threshold)
    fracB.q.prime <- mean(0==np[,'d.prime'])

    fracT.q <- mean(np[,'d']>=top.level.threshold)
    fracI.q <- mean(0<np[,'d'] & np[,'d']<top.level.threshold)
    fracB.q <- mean(0==np[,'d'])

    # Ratios of resources to consumers - p 2398, eq 11 and 12
    NP.q.prime <- 2^(-sum(vlog2v(np[,'nP']/sum(np[,'nP'])), na.rm=TRUE)) / 
                  2^(-sum(vlog2v(np[,'nN']/sum(np[,'nN'])), na.rm=TRUE))

    NP.q <- 2^(-sum(vlog2v((np[,'bOut']*np[,'nP'])/sum(np[,'bOut']*np[,'nP'])), na.rm=TRUE)) / 
            2^(-sum(vlog2v((np[,'bIn']*np[,'nN'])/sum(np[,'bIn']*np[,'nN'])), na.rm=TRUE))

    # Link properties
    # p 2398, eq 13
    LD.q.prime <- (sum(np[,'nP']) + sum(np[,'nN'])) / (2*NumberOfNodes(community))
    # p 2398, eq 14
    LD.q <- (sum(np[,'bOut']*np[,'nP']/sumb, na.rm=TRUE) + 
             sum(np[,'bIn']*np[,'nN']/sumb, na.rm=TRUE))/2

    # p 2398, col 2
    C.q.prime <- LD.q.prime/NumberOfNodes(community)
    C.q <- LD.q/NumberOfNodes(community)

    # ...the sum of the diversity of outflows weighted by the total outflows, 
    # of the diversity of inflows weighted by the total inflows. Phi can be 
    # thought of as the average amount of choice in trophic pathways (Ulanowicz 
    # and Wolff 1991)
    Phi <- sum(HP*np[,'bOut']/sumb) + sum(HN*np[,'bIn']/sumb) # p 2398, eq 16
    m <- 2^(Phi/2) # Effective connectance per node p 2398, eq 15

    # Page 2398, eq 18
    PhiAB <- function(A, B)
    {
        # A and B should be functions that take a community and return node 
        # names or indices.
        A <- A(community)
        B <- B(community)
        bOut <- rowSums(b)[A]
        bIn <- colSums(b)[B]
        return (sum((bOut/sumb) * -vlog2v(b[A,B]/bOut),  na.rm=TRUE) + 
                sum((bIn/sumb)  * -vlog2v(t(b[A,B])/bIn), na.rm=TRUE))
    }

    fracTI.q <- PhiAB(IntermediateNodes, TopLevelNodes) / Phi
    fracTB.q <- PhiAB(BasalNodes, TopLevelNodes) / Phi
    fracII.q <- PhiAB(IntermediateNodes, IntermediateNodes) / Phi
    fracIB.q <- PhiAB(BasalNodes, IntermediateNodes) / Phi

    # Page 2399, top left
    PhiPrime <- sum(HP/NumberOfNodes(community)) + sum(HN/NumberOfNodes(community))
    PhiABPrime <- function(A, B)
    {
        # A and B should be functions that take a community and return node 
        # names or indices.
        A <- A(community)
        B <- B(community)
        bOut <- rowSums(b)[A]
        bIn <- colSums(b)[B]
        s <- NumberOfNodes(community)
        return (sum((1/s) * -vlog2v(b[A,B]/bOut),  na.rm=TRUE) + 
                sum((1/s) * -vlog2v(t(b[A,B])/bIn), na.rm=TRUE))
    }

    fracTI.q.prime <- PhiABPrime(IntermediateNodes, TopLevelNodes) / PhiPrime
    fracTB.q.prime <- PhiABPrime(BasalNodes, TopLevelNodes) / PhiPrime
    fracII.q.prime <- PhiABPrime(IntermediateNodes, IntermediateNodes) / PhiPrime
    fracIB.q.prime <- PhiABPrime(BasalNodes, IntermediateNodes) / PhiPrime

    # Weighted chain lengths, p 2399, eq 19 and 20
    tlps <- TLPS(community, link.properties=weight)
    chains <- TrophicChains(community)

    # The biomass flowing through every link in every chain
    bc <- apply(as.matrix(chains), 1, function(r)
    {
        links <- r[which(""!= r)]
        sapply(2:length(links), function(index)
        {
            tlps[tlps$resource==links[index-1] & tlps$consumer==links[index],3]
        })
    })
    cl.q.prime <- 2^(-sapply(bc, function(b.chain) sum(vlog2v(b.chain/sum(b.chain)))))
    bc.mean <- sapply(bc, sum)/cl.q.prime
    cl.q <- cl.q.prime*bc.mean*nrow(chains) / sum(bc.mean)

    # Quantitative generality and vulnerability, p 2400, eq 24-27
    nT <- length(TopLevelNodes(community))
    nI <- length(IntermediateNodes(community))
    nB <- length(BasalNodes(community))
    G.q.prime <- sum(np[,'nN'])/(nT+nI)
    G.q <- sum(np[,'nN']*np[,'bIn']/sumb, na.rm=TRUE)
    V.q.prime <- sum(np[,'nP'])/(nI+nB)
    V.q <- sum(np[,'nP']*np[,'bOut']/sumb, na.rm=TRUE)

    # Bersier et al table 2
    tlps <- TLPS(community, node.properties=c('IsTopLevelNode', 'IsIntermediateNode', 'IsBasalNode'))
    fracTI <- with(tlps, sum(resource.IsIntermediateNode & consumer.IsTopLevelNode))/nrow(tlps)
    fracTB <- with(tlps, sum(resource.IsBasalNode & consumer.IsTopLevelNode))/nrow(tlps)
    fracII <- with(tlps, sum(resource.IsIntermediateNode & consumer.IsIntermediateNode))/nrow(tlps)
    fracIB <- with(tlps, sum(resource.IsBasalNode & consumer.IsIntermediateNode))/nrow(tlps)
    # This is far faster than ChainLengths
    chains.length <- TrophicChainsStats(community)$chain.lengths
    Qualitative <- c(FractionTopLevelNodes(community), 
                     FractionIntermediateNodes(community), 
                     FractionBasalNodes(community), 
                     sum(NumberOfConsumers(community)>0) / sum(NumberOfResources(community)>0), 
                     LinkageDensity(community), 
                     DirectedConnectance(community), 
                     fracTI, 
                     fracTB, 
                     fracII, 
                     fracIB, 
                     mean(chains.length), 
                     median(chains.length), 
                     sd(chains.length), 
                     max(chains.length), 
                     Omnivory(community, level=ChainAveragedTrophicLevel),
                     mean(TrophicGenerality(community)[NumberOfResources(community)>0]), 
                     mean(TrophicVulnerability(community)[NumberOfConsumers(community)>0]), 
                     sd(NormalisedTrophicGenerality(community)), 
                     sd(NormalisedTrophicVulnerability(community)))

    Unweighted <- c(fracT.q.prime, 
                    fracI.q.prime, 
                    fracB.q.prime, 
                    NP.q.prime, 
                    LD.q.prime,
                    C.q.prime, 
                    fracTI.q.prime, 
                    fracTB.q.prime, 
                    fracII.q.prime, 
                    fracIB.q.prime, 
                    mean(cl.q.prime), 
                    median(cl.q.prime), 
                    sd(cl.q.prime), 
                    max(cl.q.prime), 
                    mean(np[,'o.prime']), 
                    G.q.prime, 
                    V.q.prime, 
                    sd(np[,'g.prime']), 
                    sd(np[,'v.prime']))

    Weighted <- c(fracT.q, 
                  fracI.q, 
                  fracB.q, 
                  NP.q, 
                  LD.q,
                  C.q, 
                  fracTI.q, 
                  fracTB.q, 
                  fracII.q, 
                  fracIB.q, 
                  mean(cl.q), 
                  median(cl.q), 
                  sd(cl.q), 
                  max(cl.q), 
                  mean(np[,'o']), 
                  G.q, 
                  V.q, 
                  sd(np[,'g']), 
                  sd(np[,'v']))
    res <- cbind(Qualitative, Unweighted, Weighted)

    rownames(res) <- c('Fraction top level', 
        'Fraction intermediate', 'Fraction basal', 'Ratio resources:consumers', 
        'Link density', 'Connectance', 'Fraction links top:intermediate', 
        'Fraction links top:basal', 'Fraction links intermediate:intermediate', 
        'Fraction links intermediate:basal', 'Mean chain length', 
        'Median chain length', 'SD chain length', 'Max chain length', 
        'Degree of omnivory', 'Generality', 'Vulnerability', 
        'SD standardised generality', 'SD standardised vulnerability')
    return (res)
}
