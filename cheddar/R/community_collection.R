# Collections of Community objects
is.CommunityCollection <- function(x)
{
    return (inherits(x, 'CommunityCollection'))
}

"[<-.CommunityCollection" <- function(x, i, value)
{
    if(!is.CommunityCollection(x)) stop('Not a CommunityCollection')
    stop("Can't assign to a CommunityCollection")
}

"[[<-.CommunityCollection" <- function(x, i, value)
{
    if(!is.CommunityCollection(x)) stop('Not a CommunityCollection')
    stop("Can't assign to a CommunityCollection")
}

'$<-.CommunityCollection' <- function(object, x, value)
{
    if(!is.CommunityCollection(object)) stop('Not a CommunityCollection')
    stop("Can't assign to a CommunityCollection")
}

"names<-.CommunityCollection" <- function(x, value)
{
    if(!is.CommunityCollection(x)) stop('Not a CommunityCollection')
    stop("Can't assign to a CommunityCollection")
}

"length<-.CommunityCollection" <- function(x, value)
{
    if(!is.CommunityCollection(x)) stop('Not a CommunityCollection')
    stop("Can't change length of a CommunityCollection")
}

"levels<-.CommunityCollection" <- function(x, value)
{
    if(!is.CommunityCollection(x)) stop('Not a CommunityCollection')
    stop("Can't change levels of a CommunityCollection")
}

"dim<-.CommunityCollection" <- function(x, value)
{
    if(!is.CommunityCollection(x)) stop('Not a CommunityCollection')
    stop("Can't change dim of a CommunityCollection")
}

CommunityCollection <- function(communities)
{
    if(is.CommunityCollection(communities))
    {
        return (communities)
    }
    else if('list'==class(communities) && length(communities)>0)
    {
        # Are all the elements of communities Community objects?
        ok <- sapply(communities, is.Community)
        if(any(!ok))
        {
            stop(paste('The list elements [', paste(which(!ok), collapse=','), 
                       '] are not community objects', sep=''))
        }

        # Are any communities duplicated?
        names(communities) <- sapply(communities, function(c) CP(c, 'title'))

        bad <- duplicated(names(communities))
        if(any(bad))
        {
            stop(paste('More than one community has the same name [', 
                       paste(names(communities)[bad], collapse=','), ']', 
                       sep=''))
        }

        # Units should be consistent
        M.units <- sapply(communities, CP, property='M.units')
        M.units <- unique(M.units)
        if(length(M.units)>1)
        {
            stop(paste('The communities have inconsistent M.units [', 
                       paste(M.units, collapse=','), ']', 
                       sep=''))
        }

        N.units <- unique(unlist(sapply(communities, CP, property='N.units')))
        if(length(N.units)>1)
        {
            stop(paste('The communities have inconsistent N.units [', 
                       paste(N.units, collapse=','), ']', 
                       sep=''))
        }

        self <- communities
        class(self) <- c('CommunityCollection', 'list')
        return (self)
    }
    else if('list'==class(communities) && 0==length(communities))
    {
        stop('The list of communities is empty')
    }
    else
    {
        stop(paste('The communities must be in a list. Unable to create ', 
                   'community collection from [', class(communities), ']', 
                   sep=''))
    }
}

print.CommunityCollection <- function(x, ...)
{
    if(!is.CommunityCollection(x)) stop('Not a CommunityCollection')
    cat(paste('A collection of', length(x), 'communities\n'))
    invisible(x)
}

summary.CommunityCollection <- function(object, ...)
{
    if(!is.CommunityCollection(object)) stop('Not a CommunityCollection')
    return (do.call('.SimpleRBindFill', lapply(object, function(o) as.data.frame(summary(o)))))
}

"[.CommunityCollection" <- function(x, i)
{
    if(!is.CommunityCollection(x)) stop('Not a CommunityCollection')
    res <- NextMethod(x,i)
    if(0==length(res))  return (NULL)
    else                return (CommunityCollection(res))
}

CollectionNPS <- function(collection, properties=NULL)
{
    # A data.frame containing a row for every node in every community
    if(!is.CommunityCollection(collection)) stop('Not a CommunityCollection')

    # 'node' must be in the output. Can't use union() - it nukes any names 
    # that properties has.
    if(!is.null(properties) && !'node' %in% properties)
    {
        properties <- c('node', properties)
    }

    # Should be OK to call rbind here as I can't think of a case where the 
    # number of cols returned by NPS() would be different for different 
    # communities in the collection. Call .SimpleRBindFill just in case.
    p <- do.call('.SimpleRBindFill', lapply(collection, function(community)
    {
        np <- NPS(community, properties=properties)
        stopifnot(!'community' %in% colnames(np))
        return (cbind(community=CP(community, 'title'), 
                      np, 
                      stringsAsFactors=FALSE))
    }))
    rownames(p) <- NULL
    return (p)
}

CollectionTLPS <- function(collection, node.properties=NULL, 
                                            link.properties=NULL)
{
    # A data.frame containing a row for every trophic link in every community
    if(!is.CommunityCollection(collection)) stop('Not a CommunityCollection')
    p <- lapply(collection, function(community)
    {
        tlp <- TLPS(community, node.properties=node.properties, 
                    link.properties=link.properties)
        if(is.null(tlp))
        {
            return (NULL)
        }
        else
        {
            stopifnot(!'community' %in% colnames(tlp))
            return (cbind(community=CP(community, 'title'), tlp, 
                          stringsAsFactors=FALSE))
        }
    })

    # Elements will be NULL for communities that lack trophic links
    p <- p[!sapply(p, is.null)]

    # Should be OK to call rbind here as I can't think of a case where the 
    # number of cols returned by TLPS() would be different for different 
    # communities in the collection. Call .SimpleRBindFill just in case.
    
    if(length(p)>0)
    {
        p <- do.call('.SimpleRBindFill', p)
        rownames(p) <- NULL
        return (p)
    }
    else
    {
        return (NULL)
    }
}

CollectionCPS <- function(collection, properties=NULL)
{
    # Returns a data.frame containing a row for each community in collection
    # properties - a vector of community properties
    if(!is.CommunityCollection(collection)) stop('Not a CommunityCollection')
    if(is.null(properties))
    {
        # This solution also works for the non-NULL case but gets column orders 
        # wrong if some properties are functions that return > 1 item 
        # and the number of items is different for different communities in 
        # the collection.
        res <- lapply(collection, function(community)
        {
            return (data.frame(CPS(community), stringsAsFactors=FALSE, 
                               check.names=FALSE))
        })

        res <- do.call('.SimpleRBindFill', res)
    }
    else
    {
        # A more complicated solution than the one above but one that gets 
        # column orders (mostly) correct.
        res <- NULL

        for(index in 1:length(properties))
        {
            if(is.list(properties)) p <- properties[index]
            else                    p <- list(properties[index])

            if(!is.null(names(properties)))
            {
                names(p) <- names(properties)[index]
            }
            cols <- lapply(collection, CPS, properties=p)
            cols <- lapply(cols, data.frame, stringsAsFactors=FALSE, 
                           check.names=FALSE)

            # Order by number of cols - not guaranteed to get the 'correct' 
            # order in all cases but works for some.
            o <- order(sapply(cols, ncol), decreasing=TRUE)
            cols <- cols[o]
            cols <- do.call('.SimpleRBindFill', cols)
            cols <- cols[order(o),,drop=FALSE]
            if(is.null(res)) res <- cols
            else             res <- cbind(res, cols)
        }
    }
    return (res)
}

AggregateCommunitiesBy <- function(collection, aggregate.by, ...)
{
    # Returns a new CommunityCollection object
    if(!is.CommunityCollection(collection)) stop('Not a CommunityCollection')

    aggregate <- CollectionCPS(collection, aggregate.by)[,aggregate.by]
    stopifnot(length(aggregate)==length(collection))

    aggregated <- NULL
    for(ab in unique(aggregate))
    {
        a <- AggregateCommunities(collection[aggregate==ab])
        aggregated[[1+length(aggregated)]] <- a
    }

    aggregated <- do.call('list', aggregated)
    if(1==length(aggregated))   return (aggregated[[1]])
    else                        return (CommunityCollection(aggregated))
}

AggregateCommunities <- function(collection, 
                                 aggregate=names(collection), 
                                 weight.by='N', 
                                 title=NULL)
{
    # Returns a new Community object

    if(!is.CommunityCollection(collection)) stop('Not a CommunityCollection')

    if(is.null(title))
    {
        title <- paste('Aggregation of', paste(aggregate, collapse=','))
    }

    if(is.character(aggregate))
    {
        bad <- !aggregate %in% names(collection)
        if(any(bad))
        {
            stop(paste('The names [', paste(aggregate[bad], collapse=','), 
                       '] are not names of communities in the collection', 
                       sep=''))
        }
    }
    else
    {
        bad <- aggregate<1 | aggregate>length(collection)
        if(any(bad))
        {
            stop(paste('The indices [', paste(aggregate[bad], collapse=','), 
                       '] are not in the collection', 
                       sep=''))
        }
        aggregate <- names(collection)[aggregate]
    }

    # Aggregate nodes
    # A data.frame containing all node properties in across all contained webs
    node.properties <- CollectionNPS(collection)

    # Take the subset that we are interested in
    rows.to.keep <- node.properties$community %in% aggregate
    node.properties <- node.properties[rows.to.keep,]

    # A union of all nodes in aggregate
    new.nodes <- data.frame(node=unique(node.properties$node), 
                            stringsAsFactors=FALSE)

    # Aggregate node properties
    if(ncol(node.properties)>2)
    {
        # Add a row for species X where X is not in community A
        setdiff.data.frame <- function(A, B)
        {
            # setdiff for data.frames
            # https://stat.ethz.ch/pipermail/r-devel/2007-December/047708.html
            a <- do.call("paste", c(A, sep = "\r"))
            b <- do.call("paste", c(B, sep = "\r"))
            return (A[match(setdiff(a, b),a), ])
        }

        possible <- expand.grid(community=unique(node.properties$community), 
                                node=unique(node.properties$node), 
                                stringsAsFactors=FALSE)
        actual <- node.properties[,c('community', 'node')]
        missing <- setdiff.data.frame(possible, actual)

        if(nrow(missing)>0)
        {
            # Add rows with missing numeric and integer values to 0
            before <- nrow(node.properties)
            node.properties <- .SimpleRBindFill(node.properties, missing)
            after <- nrow(node.properties)
            classes <- sapply(node.properties, class)
            node.properties[(1+before):after,classes %in% c('numeric', 'integer')] <- 0
#            node.properties[(1+before):after,classes %in% c('character')] <- ''
        }

        select.cols <- !colnames(node.properties) %in% c('community', 'node')
        new.nodes <- cbind(new.nodes, 
                       .AggregateDataFrame(data=node.properties[,select.cols,drop=FALSE],
                                           aggregate.by=node.properties$node, 
                                           weight.by=weight.by))
    }

    # Aggregate trophic links
    new.tl <- CollectionTLPS(collection)
    if(!is.null(new.tl))
    {
        new.tl <- new.tl[new.tl$community %in% aggregate,]
        select.cols <- !colnames(new.tl) %in% 'community'
        new.tl <- .AggregateDataFrame(data=new.tl[,select.cols,drop=FALSE], 
                  aggregate.by=paste(new.tl[,'resource'],new.tl[,'consumer']),
                  weight.by=NULL)
        stopifnot(!any(duplicated(new.tl)))
    }

    # Aggregate community properties
    new.properties <- CollectionCPS(collection)
    new.properties <- new.properties[new.properties$title %in% aggregate,,
                                     drop=FALSE]
    if(ncol(new.properties)>1)
    {
        select.cols <- !colnames(new.properties) %in% 'title'
        new.properties <- .AggregateDataFrame(data=new.properties[,select.cols,drop=FALSE],
                              aggregate.by=rep(1, nrow(new.properties)), 
                              weight.by=NULL)
        stopifnot(1==nrow(new.properties))
        new.properties <- do.call(list, new.properties)
        new.properties$title <- title
    }
    else
    {
        new.properties <- list(title=title)
    }

    return (Community(nodes=new.nodes, 
                      trophic.links=new.tl, 
                      properties=new.properties))
}

plot.CommunityCollection <- function(x, 
                                     ncol=min(length(x),5), 
                                     by.col=TRUE, 
                                     plot.fn=plot, 
                                     ...)
{
    nrow <- ceiling(length(x)/ncol)
    if(by.col)
    {
        mfrow <- par(mfrow=c(nrow, ncol))
        on.exit(par(mfrow))
    }
    else
    {
        mfcol <- par(mfcol=c(nrow, ncol))
        on.exit(par(mfcol))
    }

    junk <- sapply(x, function(community)
    {
        tryCatch(do.call(plot.fn, c(list(community), list(...))))
    })
}

LoadCollection <- function(dir, ...)
{
    # Load each community directory
    path <- file.path(dir, 'communities')
    if(!file.exists(path) || !file.info(path)$isdir)
    {
        stop(paste('The community collection directory [', path, 
                   '] does not exist', sep=''))
    }
    else
    {
        files <- list.files(path)
        collection <- lapply(file.path(path, files), LoadCommunity, ...)

        return (CommunityCollection(collection))
    }
}

SaveCollection <- function(collection, dir, ...)
{
    if(!is.CommunityCollection(collection)) stop('Not a CommunityCollection')
    if(file.exists(dir))
    {
        stop(paste('The directory [', dir, '] already exists', sep=''))
    }
    else
    {
        dir.create(file.path(dir, 'communities'), recursive=TRUE)

        # Assignment to junk prevents result of mapply() being returned
        junk <- mapply(SaveCommunity, community=collection, 
                       dir=file.path(dir, 'communities', names(collection)), 
                       MoreArgs=list(...))
    }
}

OrderCollection <- function(collection, ..., decreasing=FALSE)
{
    # Returns a CommunityCollection containing communities in collection 
    # ordered by the community properties given in ...
    if(!is.CommunityCollection(collection)) stop('Not a CommunityCollection')

    # Get the properties that we will be sorting on
    properties <- unlist(list(...))
    ccp <- CollectionCPS(collection, properties=properties)
    sort.by <- as.list(ccp)[properties]

    # Delegate to order()
    params <- c(sort.by, list(decreasing=decreasing))
    o <- do.call('order', params)
    return (CommunityCollection(collection[o]))
}

subset.CommunityCollection <- function(x, subset, properties=NULL, ...)
{
    # Returns a subset of the collection
    # select and drop are ignored
    if(!is.CommunityCollection(x)) stop('Not a CommunityCollection')

    # dots can be used to fetch properties used by subset
    ccp <- CollectionCPS(x, properties=properties)

    subset <- substitute(subset)
    s <- with(ccp, eval(subset))
    if(!is.logical(s))
    {
        stop('subset should evaluate to a logical')
    }

    return (x[s])
}

CollectionApply <- function(collection, f, ...)
{
    if(!is.CommunityCollection(collection)) stop('Not a CommunityCollection')
    return(CommunityCollection(lapply(collection, f, ...)))
}

SiteBySpeciesMatrix <- function(collection, abundance=NULL, na.missing=FALSE)
{
    # Returns a matrix with a column per community and a row per unique node 
    # within all communities in the collection.
    if (!is.CommunityCollection(collection)) 
        stop("Not a CommunityCollection")
    # The names of all nodes in the collection
    cps <- CollectionNPS(collection, c('node', abundance))
    nodes <- sort(unique(unname(cps[,'node'])))

    res <- matrix(ifelse(na.missing, NA, 0), 
                  ncol=length(collection), nrow=length(nodes), 
                  dimnames=list(nodes, names(collection)))

    decode.node <- 1:length(nodes)
    names(decode.node) <- nodes

    decode.site <- 1:length(collection)
    names(decode.site) <- names(collection)

    indices <- decode.node[cps$node] + 
               length(nodes)*(decode.site[cps$community]-1)
    if(is.null(abundance)) res[indices] <- 1
    else                   res[indices] <- cps[,abundance]
    return (res)
}

NvMTriTrophicTable <- function(collection)
{
    # Cohen et al 2009 PNAS Table 1
    if(!is.CommunityCollection(collection)) stop('Not a CommunityCollection')

    UnsafeMean <- function(x)
    {
        ifelse(is.null(x), NA, mean(x))
    }

    res <- lapply(collection, function(community)
    {
        # NvMTriTrophicStatistics() removes both nodes lacking M and/or N and 
        # cannibalistic links. Remove these here too so that the network statistics 
        # (L, S^2, L/S^2 and L/S) are consistent with the rest of the values in the 
        # table.
        unadjusted <- community
        community <- RemoveNodes(community, remove=with(NPS(community), node[is.na(M) | is.na(N)]))
        community <- RemoveCannibalisticLinks(community)
        community <- RemoveIsolatedNodes(community)

        tts <- NvMTriTrophicStatistics(community)
        lp <- tts[['links']]
        tncp <- tts[['three.node.chains']]
        tcp <- tts[['trophic.chains']]
        
        community.span <- diff(range(Log10M(community))) +
                          diff(range(Log10N(community)))
                          
        wiggling <- mean(tcp$sum.chain.length) / mean(tcp$chain.span)

        return (c(
         'Mean link length'=mean(lp$length),
         'Mean L upper'=UnsafeMean(tncp$Lupper),
         'Mean L lower'=UnsafeMean(tncp$Llower),
         '2 x mean link length'=2*mean(lp$length),
         'Mean 2-span'=UnsafeMean(tncp$two.span),
         'Mean L upper + L lower'=UnsafeMean(tncp$Lupper + tncp$Llower),
         '2 x mean link length / mean 2-span'=2 * mean(lp$length) / UnsafeMean(tncp$two.span),
         'Mean L upper + L lower/ mean 2-span'=UnsafeMean(tncp$Lupper + tncp$Llower) / UnsafeMean(tncp$two.span),
         'Mean count chain length'=mean(tcp$count.chain.length),
         'Mean count chain length x mean link length'=mean(tcp$count.chain.length)*mean(lp$length),
         'Community span'=community.span,
         'Mean count chain length x mean link length / community span'=mean(tcp$count.chain.length)*mean(lp$length)/community.span,
         'Mean sum chain lengths'=mean(tcp$sum.chain.length),
         'Mean chain span'=mean(tcp$chain.span),
         'Mean chain span / community span'=mean(tcp$chain.span) / community.span,
         'Mean sum chain lengths / mean chain span'=mean(tcp$sum.chain.length) / mean(tcp$chain.span),
         'Mean sum chain lengths / community span'=mean(tcp$sum.chain.length) / community.span,
         'L'=NumberOfTrophicLinks(community),
         'S^2'=NumberOfNodes(community)^2,
         'L/S^2'=DirectedConnectance(community),
         'L/S'=LinkageDensity(community), 
         'L'=NumberOfTrophicLinks(unadjusted),
         'S^2'=NumberOfNodes(unadjusted)^2,
         'L/S^2'=DirectedConnectance(unadjusted),
         'L/S'=LinkageDensity(unadjusted)))
    })
    res <- do.call('cbind', res)
    colnames(res) <- names(collection)
    return (res)
}
