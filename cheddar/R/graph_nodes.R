# High-level and helper functions for plotting community nodes
.PlotTrophicLinksAsNetwork <- function(community, x, y, 
                                       link.colour.by,link.colour.spec,link.col,
                                       link.line.type.by, link.line.type.spec, 
                                       link.lty, link.lwd, highlight.links, ...)
{
    # A private helper function that assembles parameters and calls .PlotWeb().
    # The food web is plotted using link.col, link.lwd and link.lty, which are 
    # resolved according to specs.
    if(is.function(highlight.links))
    {
        highlight.links <- highlight.links(community)
    }

    if(length(highlight.links)>0)
    {
        highlight.links <- .ResolveTrophicLinksToRowIndices(community, 
                                                            highlight.links)
        if(any(is.na(highlight.links)))
        {
            stop(paste("One or more of the links given in", 
                       "'highlight.links' is not present in the community"))
        }
    }
    else
    {
        highlight.links <- NULL
    }

    link.col <- .TrophicLinkGraphParamFromSpec(community, 'link.col', 
                                               link.col, 
                                               link.colour.by, 
                                               link.colour.spec, 
                                               default=DefaultLinkColour())
    link.col <- .HighlightColours(NumberOfTrophicLinks(community), 
                                  link.col, highlight.links, 
                                  NULL)

    link.lty <- .TrophicLinkGraphParamFromSpec(community, 'link.lty', 
                                               link.lty, 
                                               link.line.type.by, 
                                               link.line.type.spec)

    # Get resource and consumer as node indices
    network <- TLPS(community, link.properties=c('resource', 'consumer'))
    network <- as.matrix(network)
    dim(network) <- NULL
    network <- .ResolveToNodeIndices(community, network)
    dim(network) <- c(length(network)/2, 2)

    # Strip out the names that we are using
    args <- list(x=x, y=y, network=network, col=link.col, lty=link.lty, 
                 lwd=link.lwd, highlight=highlight.links)
    dots <- list(...)
    args <- c(args, dots[setdiff(names(dots), names(args))])
    do.call(.PlotNetwork, args)
}

.PlotNodes <- function(community, x, y, 
                       colour.by=NULL, colour.spec=NULL, col=NULL, 
                       symbol.by=NULL, symbol.spec=NULL, pch=NULL, 
                       bg.by=NULL, bg.spec=NULL, bg=NULL, 
                       cex.by=NULL, cex.spec=NULL, cex=NULL, 
                       label.colour.by=NULL, label.colour.spec=NULL, 
                       label.colour=NULL, 
                       highlight.nodes=NULL, lowlight.nodes=NULL, 
                       show.nodes.as='points', 
                       node.labels=NULL, label.cex=0.6, 
                       ...)
{
    # A private helper that assembles graphical params and calls 
    # .PlotHighlightedPoints()
    # Resolves node col, pch and bg according to specs.
    # show.nodes.as should be either 'symbols', 'labels' or 'both'
    # If show.nodes.as is=='labels' and node.labels is NULL, labels are 
    # 1:NumberOfNodes().
    # If show.nodes.as=='labels' labels are plotted using label.cex and 
    # label.colour. 
    # Nodes to be highlighted or lowlighted have col, pch, bg and cex 
    # adjusted accordingly.
    stopifnot(show.nodes.as %in% c('points', 'labels', 'both'))

    # col
    if(missing(colour.by))
    {
        colour.by <- NULL
        if('category' %in% NodePropertyNames(community))
        {
            colour.by <- 'category'
        }
    }

    if(missing(colour.spec))
    {
        colour.spec <- NULL
        if(!is.null(colour.by) && 'category'==colour.by)
        {
            colour.spec <- DefaultCategoryColours()
        }
    }

    # bg
    if(missing(bg.by))
    {
        bg.by <- NULL
        if('category' %in% NodePropertyNames(community))
        {
            bg.by <- 'category'
        }
    }

    if(missing(bg.spec))
    {
        bg.spec <- NULL
        if(!is.null(bg.by) && 'category'==bg.by)
        {
            bg.spec <- DefaultCategoryColours()
        }
    }

    # pch
    if(missing(symbol.by))
    {
        symbol.by <- NULL
        if('category' %in% NodePropertyNames(community))
        {
            symbol.by <- 'category'
        }
    }

    if(missing(symbol.spec))
    {
        symbol.spec <- NULL
        if(!is.null(symbol.by) && 'category'==symbol.by)
        {
            symbol.spec <- DefaultCategorySymbols()
        }
    }

    # Nodes
    if(is.logical(highlight.nodes) && 1==length(highlight.nodes) && 
       !highlight.nodes)
    {
        # The use probably meant to set highlight.nodes=NULL
        highlight.nodes <- NULL
    }

    if(is.function(highlight.nodes))
    {
        highlight.nodes <- highlight.nodes(community)
    }
    highlight.nodes <- .ResolveToNodeIndices(community, highlight.nodes)

    if(is.logical(lowlight.nodes) && 1==length(lowlight.nodes) && 
       !lowlight.nodes)
    {
        # The use probably meant to set lowlight.nodes=NULL
        lowlight.nodes <- NULL
    }
    if(is.function(lowlight.nodes))
    {
        lowlight.nodes <- lowlight.nodes(community)
    }
    lowlight.nodes <- .ResolveToNodeIndices(community, lowlight.nodes)
    lowlight.nodes <- setdiff(lowlight.nodes, highlight.nodes)

    # col can not be NULL
    col <- .NodeGraphParamFromSpec(community, 'col', col, colour.by, colour.spec,
                                   default=par('col'))
    col <- .HighlightColours(NumberOfNodes(community), col, highlight.nodes, 
                             lowlight.nodes)

    cex <- .NodeGraphParamFromSpec(community, 'cex', cex, cex.by, cex.spec)
    cex <- .HighlightCex(NumberOfNodes(community), cex, highlight.nodes, 
                         lowlight.nodes)

    # Labels
    if(show.nodes.as %in% c('labels', 'both'))
    {
        node.labels <- .GraphParNodeLabels(community, node.labels)

        if('both'==show.nodes.as)
        {
            # Only if plotting points and labels together do we use 
            # label.colour.

            # col can not be NULL
            label.colour <- .NodeGraphParamFromSpec(community, 'label.colour', 
                                                    label.colour,label.colour.by,
                                                    label.colour.spec, 
                                                    default=par('col'))
            label.colour <- .HighlightLabelColours(NumberOfNodes(community), 
                                                   label.colour, highlight.nodes,
                                                   lowlight.nodes)
        }
        else if('labels'==show.nodes.as)
        {
            # If we are showing only labels, use the normal cex and colours
            label.colour <- col
            label.cex <- cex
        }
    }

    # Points
    if(show.nodes.as %in% c('points', 'both'))
    {
        pch <- .NodeGraphParamFromSpec(community, 'pch', pch, symbol.by, 
                                      symbol.spec)
        pch <- .HighlightSymbols(NumberOfNodes(community), pch, highlight.nodes, 
                                 lowlight.nodes)

        bg <- .NodeGraphParamFromSpec(community, 'bg', bg, bg.by, bg.spec)
        bg <- .HighlightBG(NumberOfNodes(community), bg, highlight.nodes, 
                           lowlight.nodes)
    }

    .PlotHighlightedPoints(x, y, col=col, pch=pch, bg=bg, cex=cex, 
                           labels=node.labels, label.col=label.colour, 
                           label.cex=label.cex, 
                           show.points=show.nodes.as %in% c('points','both'), 
                           show.labels=show.nodes.as %in% c('labels','both'), 
                           highlight=highlight.nodes, ...)
}

PlotNPS <- function(community, 
                    X,
                    Y,
                    main=CPS(community)$title, 
                    xlab, 
                    ylab, 
                    xlim=NULL, 
                    ylim=NULL, 
                    colour.by,
                    colour.spec,
                    col=NULL, 
                    symbol.by,
                    symbol.spec,
                    pch=NULL, 
                    bg.by,
                    bg.spec,
                    bg=NULL, 
                    cex.by=NULL,
                    cex.spec=NULL,
                    cex=NULL, 
                    label.colour.by=NULL,
                    label.colour.spec=NULL,
                    label.colour=NULL, 
                    link.colour.by=NULL, 
                    link.colour.spec=NULL, 
                    link.col=NULL, 
                    link.line.type.by=NULL, 
                    link.line.type.spec=NULL, 
                    link.lty=NULL, 
                    link.lwd=NULL, 
                    highlight.links=NULL,
                    highlight.nodes=Cannibals,
                    lowlight.nodes, 
                    show.na=FALSE,
                    show.web=TRUE,
                    show.nodes.as='points', 
                    node.labels=NULL, 
                    label.cex=0.6, 
                    are.values=FALSE, 
                    frame.plot=TRUE, 
                    ...)
{
    if(!is.Community(community)) stop('Not a Community')
    if(are.values)
    {
        if(missing(xlab))   xlab <- ''
        if(missing(ylab))   ylab <- ''

        stopifnot(length(X)==NumberOfNodes(community))
        stopifnot(length(Y)==NumberOfNodes(community))
        names(X) <- names(Y) <- NP(community, 'node')
    }
    else
    {
        if(missing(xlab))   xlab <- X
        if(missing(ylab))   ylab <- Y

        X <- NPS(community, X)[,X]
        Y <- NPS(community, Y)[,Y]
    }

    if(missing(lowlight.nodes))
    {
        lowlight.nodes <- which(is.na(X) | is.na(Y))
    }

    if(show.na)
    {
        points <- PlaceMissingPoints(X, xlim, Y, ylim)
        X <- points[,1]
        Y <- points[,2]
    }

    plot(X, Y, type='n', main=main, xlab=xlab, ylab=ylab, 
         xlim=xlim, ylim=ylim, frame.plot=frame.plot, ...)
    .AddAxisTicks(...)
    if(show.web && !is.null(TLPS(community)))
    {
        .PlotTrophicLinksAsNetwork(community=community, x=X, y=Y, 
                                   link.colour.by=link.colour.by, 
                                   link.colour.spec=link.colour.spec, 
                                   link.col=link.col, 
                                   link.line.type.by=link.line.type.by, 
                                   link.line.type.spec=link.line.type.spec, 
                                   link.lty=link.lty, link.lwd=link.lwd, 
                                   highlight.links=highlight.links, ...)
    }
    
    .PlotNodes(community=community, x=X, y=Y, 
               colour.by=colour.by, colour.spec=colour.spec, col=col, 
               symbol.by=symbol.by, symbol.spec=symbol.spec, pch=pch, 
               bg.by=bg.by, bg.spec=bg.spec, bg=bg, 
               cex.by=cex.by, cex.spec=cex.spec, cex=cex, 
               label.colour.by=label.colour.by, 
               label.colour.spec=label.colour.spec, 
               label.colour=label.colour, 
               highlight.nodes=highlight.nodes, 
               lowlight.nodes=lowlight.nodes, 
               show.nodes.as=show.nodes.as, 
               node.labels=node.labels, label.cex=label.cex, 
               ...)
}

PlotNvM <- function(community, 
                    xlab=Log10MLabel(community), 
                    ylab=Log10NLabel(community), 
                    ...)
{
    if(!is.Community(community)) stop('Not a Community')

    .RequireM(community)
    .RequireN(community)
    PlotNPS(community, 'Log10M', 'Log10N', xlab=xlab, ylab=ylab, 
            are.values=FALSE, ...)
}

PlotMvN <- function(community, 
                    xlab=Log10NLabel(community), 
                    ylab=Log10MLabel(community), 
                    ...)
{
    if(!is.Community(community)) stop('Not a Community')

    .RequireM(community)
    .RequireN(community)
    PlotNPS(community, 'Log10N', 'Log10M', xlab=xlab, ylab=ylab, 
            are.values=FALSE, ...)
}

PlotBvM <- function(community, 
                    xlab=Log10MLabel(community), 
                    ylab=Log10BLabel(community), 
                    ...)
{
    if(!is.Community(community)) stop('Not a Community')

    .RequireM(community)
    .RequireN(community)
    PlotNPS(community, 'Log10M', 'Log10Biomass', xlab=xlab, ylab=ylab, 
            are.values=FALSE, ...)
}

PlotMvB <- function(community, 
                    xlab=Log10BLabel(community), 
                    ylab=Log10MLabel(community), 
                    ...)
{
    if(!is.Community(community)) stop('Not a Community')

    .RequireM(community)
    .RequireN(community)
    PlotNPS(community, 'Log10Biomass', 'Log10M', xlab=xlab, ylab=ylab, 
            are.values=FALSE, ...)
}

PlotRankNPS <- function(community, property, rank.by=property, 
                        log10.rank=FALSE, xlab, ylab, show.web=FALSE, ...)
{
    # Plots values vs rank values.
    # col, pch, bg and label.colour are reordered for rank(values)
    values <- NPS(community, property)[,property]

    # The ordering of values
    rank <- order(NPS(community, rank.by)[,rank.by], decreasing=TRUE)
    rank <- order(rank)                      # Node order
    if(log10.rank)  rank <- log10(rank)

    if(missing(xlab))
    {
        if(log10.rank)  xlab <- as.expression(substitute(log[10](Rank~r), 
                                                         list(r=rank.by)))
        else            xlab <- paste('Rank', rank.by)
    }

    if(missing(ylab))
    {
        ylab <- property
    }

    PlotNPS(community, rank, values, xlab=xlab, ylab=ylab, show.web=show.web, 
            are.values=TRUE, ...)
}

PlotMvRankM <- function(community, 
                        log10.rank=FALSE, 
                        xlab, 
                        ylab, 
                        ...)
{
    if(!is.Community(community)) stop('Not a Community')

    .RequireM(community)

    if(missing(xlab))
    {
        if(log10.rank)  xlab <- ~log[10](Rank~italic(M))
        else            xlab <- ~Rank~italic(M)
    }

    if(missing(ylab))
    {
        ylab <- Log10MLabel(community)
    }

    PlotRankNPS(community, property='Log10M', rank.by='M', 
                log10.rank=log10.rank, xlab=xlab, ylab=ylab, ...)
}

PlotNvRankN <- function(community, 
                        log10.rank=FALSE, 
                        xlab, 
                        ylab, 
                        ...)
{
    if(!is.Community(community)) stop('Not a Community')

    .RequireN(community)

    if(missing(xlab))
    {
        if(log10.rank)  xlab <- ~log[10](Rank~italic(N))
        else            xlab <- ~Rank~italic(N)
    }

    if(missing(ylab))
    {
        ylab <- Log10NLabel(community)
    }

    PlotRankNPS(community, property='Log10N', rank.by='N', 
                log10.rank=log10.rank, xlab=xlab, ylab=ylab, ...)
}

PlotBvRankB <- function(community, 
                        log10.rank=FALSE, 
                        xlab, 
                        ylab, 
                        ...)
{
    if(!is.Community(community)) stop('Not a Community')

    .RequireM(community)
    .RequireN(community)

    if(missing(xlab))
    {
        if(log10.rank)  xlab <- ~log[10](Rank~italic(B))
        else            xlab <- ~Rank~italic(B)
    }

    if(missing(ylab))
    {
        ylab <- Log10BLabel(community)
    }

    PlotRankNPS(community, property='Log10Biomass', rank.by='Biomass', 
                log10.rank=log10.rank, xlab=xlab, ylab=ylab, ...)
}

PlotNPSDistribution <- function(community, property, main=CPS(community)$title, 
                                density.args=list(), ...)
{
    if(!is.Community(community)) stop('Not a Community')

    values <- NPS(community, property)[,property]
    values <- values[!is.na(values)]
    stopifnot(!'x' %in% names(density.args))
    density.args$x <- values
    plot(do.call('density', args=density.args), main=main, ...)
    rug(values)
    .AddAxisTicks(...)
}

PlotMDistribution <- function(community, xlab=Log10MLabel(community), ...)
{
    if(!is.Community(community)) stop('Not a Community')
    .RequireM(community)
    PlotNPSDistribution(community, 'Log10M', xlab=xlab, ...)
}

PlotNDistribution <- function(community, xlab=Log10NLabel(community), ...)
{
    if(!is.Community(community)) stop('Not a Community')
    .RequireN(community)
    PlotNPSDistribution(community, 'Log10N', xlab=xlab, ...)
}

PlotBDistribution <- function(community, xlab=Log10BLabel(community), ...)
{
    if(!is.Community(community)) stop('Not a Community')
    .RequireM(community)
    .RequireN(community)
    PlotNPSDistribution(community, 'Log10Biomass', xlab=xlab, ...)
}

PlotWebByLevel <- function(community, 
                           level='PreyAveragedTrophicLevel', 
                           max.nodes.per.row=20,
                           round.levels.to.nearest=0.2,
                           stagger=0.1, 
                           x.layout='wide', # 'skinny', 'narrow' or 'wide'
                           y.layout='compress',     # 'stagger' or 'compress'
                           show.level.labels=TRUE,
                           show.level.lines=FALSE, 
                           xaxt='n', yaxt='n', xlab='', ylab='', 
                           frame.plot=FALSE, 
                           ylim=NULL, 
                           ...)
{
    # Plots the food web with higher trophic level species towards the top, 
    # similar to Jonsson et al 2005 AER Fig. 1 but without the abundance and 
    # body mass bars.
    if(!is.Community(community)) stop('Not a Community')

    # .ResolveLevel assembles a vector of length NumberOfNodse and checks that 
    # values are valid
    level <- .ResolveLevel(community, level)

    stopifnot(max.nodes.per.row>2)
    stopifnot(round.levels.to.nearest>=0 && round.levels.to.nearest<1)
    stopifnot(stagger>=0 && stagger<1)

    stopifnot(x.layout %in% c('skinny', 'narrow', 'wide'))
    if('skinny'==x.layout)
    {
        x.layout <- function(nodes.at.level, max.nodes.per.row)
        {
            # Spaced one x unit apart
            return (1:length(nodes.at.level) - (1+length(nodes.at.level))/2)
        }
    }
    else if('narrow'==x.layout)
    {
        x.layout <- function(nodes.at.level, max.nodes.per.row)
        {
            if(length(nodes.at.level)<max.nodes.per.row)
            {
                # Spaced one x unit apart
                return (1:length(nodes.at.level) - (1+length(nodes.at.level))/2)
            }
            else
            {
                # Squashed into available x space
                return (seq(-max.nodes.per.row/2, max.nodes.per.row/2, 
                            length.out=length(nodes.at.level)))
            }
        }
    }
    else if('wide'==x.layout)
    {
        x.layout <- function(nodes.at.level, max.nodes.per.row)
        {
            if(length(nodes.at.level)<max.nodes.per.row/2)
            {
                # Space widely
                extent <- max.nodes.per.row/2 - 
                          max.nodes.per.row / length(nodes.at.level) / 2
                return (seq(-extent, extent, length.out=length(nodes.at.level)))
            }
            else
            {
                # Take up all available space
                return (seq(-max.nodes.per.row/2, max.nodes.per.row/2, 
                            length.out=length(nodes.at.level)))
            }
        }
    }
    else
    {
        stop(paste('Unknown x.layout [', x.layout, ']'))
    }

    # y.layout functions return values >= 0
    stopifnot(y.layout %in% c('stagger', 'compress'))
    if('stagger'==y.layout)
    {
        y.layout <- function(nodes.at.level, max.nodes.per.row)
        {
            if(length(nodes.at.level)>=max.nodes.per.row)
            {
                # Stagger across several rows
                return (0:(length(nodes.at.level) %/% max.nodes.per.row))
            }
            else
            {
                # Everything on the same row
                return (0)
            }
        }
    }
    else if('compress'==y.layout)
    {
        y.layout <- function(nodes.at.level, max.nodes.per.row)
        {
            # Everything always on the same row
            return (0)
        }
    }
    else
    {
        stop(paste('Unknown y.layout [', y.layout, ']'))
    }

    # Assign each species x values
    x <- rep(NA, NumberOfNodes(community))
    y <- level

    # Adjust x and y to avoid overprinting

    Round <- function(v, n)
    {
        # A helper that rounds to nearest fractional part
        if(0==round.levels.to.nearest)
        {
            return (v)
        }
        else
        {
            i <- floor(v)    # integer part
            f <- v - i       # fractional part
            f <- f / n
            f <- round(f, 0)
            f <- f * n
            return (i + f)
        }
    }

    # Assign each node an x value
    while(any(is.na(x)))
    {
        l <- level[min(which(is.na(x)))]
        nodes.at.level <- which(Round(level, round.levels.to.nearest) == 
                                Round(l, round.levels.to.nearest))

        x[nodes.at.level] <- x.layout(nodes.at.level, max.nodes.per.row)
        level.offset <- y.layout(nodes.at.level, max.nodes.per.row)
        y.offset <- rep(level.offset, length.out=length(nodes.at.level))
        y.offset <- y.offset - max(y.offset)/2
        y.offset <- y.offset * stagger
        stopifnot(length(y.offset)==length(nodes.at.level))
        y[nodes.at.level] <- y[nodes.at.level] + y.offset
    }

    stopifnot(all(!is.na(x)))

    PlotNPS(community, x, y, are.values=TRUE, xlab=xlab, ylab=ylab, 
            frame.plot=frame.plot, xaxt=xaxt, yaxt=yaxt, ylim=ylim, ...)
    if(show.level.labels)
    {
        if(!is.null(ylim))
        {
           labels <- floor(min(ylim)):floor(max(ylim))
        }
        else
        {
           labels <- floor(min(level)):floor(max(level))
        }

        # Label each whole number trophic level.
        axis(2, at=unique(labels), labels=unique(labels), las=1, tick=FALSE)
    }

    if(show.level.lines)
    {
        abline(h=unique(level))
    }
}

PlotCircularWeb <- function(community, 
                            clockwise=TRUE, 
                            origin.degrees=0, 
                            proportional.radius=1, 
                            frame.plot=FALSE, 
                            xlim=c(-1,1), 
                            ylim=c(-1,1), 
                            ...)
{
    # A plot of the food web in a circular arrangement, although will be an 
    # ellipse if the graph window is not square. 
    # origin.degrees=0 will plot the first species at the 12 o'clock position. 
    if(!is.Community(community)) stop('Not a Community')

    stopifnot(proportional.radius>0 && proportional.radius<=1)

    # Create a sequence of unit vectors
    z <- complex(modulus=proportional.radius, 
                 argument=seq(0, 2*pi,length.out=1+NumberOfNodes(community)))
    z <- head(z, -1)

    # Rotate origin
    # The two calls to ifelse(clockwise,...) ensure that origin.degrees 
    # works consistently regardless of the value of clockwise. 
    # -90 is so that origin.degrees of 0 puts the origin at 12 o'clock.
    z <- z * complex(argument=(ifelse(clockwise,1,-1)*origin.degrees-90 + 
                               ifelse(clockwise,0,1)*180)*2*pi/360)

    # Adjust direction
    if(clockwise)
    {
        z <- Conj(z)
    }

    # x and y of angles
    x <- Re(z)
    y <- Im(z)

    PlotNPS(community, x, y, xlab='', ylab='', are.values=TRUE, 
            frame.plot=frame.plot, xaxt='n', yaxt='n', xlim=xlim, ylim=ylim,...)
}

PlotWagonWheel <- function(community, focus,
                           clockwise=TRUE,
                           origin.degrees=0,
                           frame.plot=FALSE,
                           main=NULL,
                           ...)
{
    # A plot of the food web in a circular arrangement, with the focal species
    # at the centre with concentric circles radiating outwards to species that
    # are 1, 2, 3 etc links away.
    # origin.degrees=0 will plot the first species at each level at the 12
    # o'clock position. 

    if(!is.Community(community)) stop('Not a Community')
    if(!is.character(focus))
    {
        stopifnot(focus>0 && focus<=NumberOfNodes(community))
        focus <- unname(NP(community, 'node'))[focus]
    }
    stopifnot(1==length(focus) && focus %in% NP(community, 'node'))

    stopifnot(!focus %in% IsolatedNodes(community))

    x <- y <- rep(NA, NumberOfNodes(community))
    names(x) <- names(y) <- unname(NP(community, 'node'))

    # Focal species at the origin
    x[focus] <- y[focus] <- 0

    x[IsolatedNodes(community)] <- y[IsolatedNodes(community)] <- NA

    # Shortest distances between nodes
    d <- ShortestPaths(community)[,focus]
    d <- d[setdiff(names(d), focus)]

    # Distance to isolated species will be Inf
    d <- d[!is.infinite(d)]

    # Should never happen
    stopifnot(!any(is.na(d)))

    # Axes limits
    max.distance <- max(d)
    lim <- c(-max.distance,max.distance)

    # For each circle, create a sequence of vectors
    for(level in sort(unique(d)))
    {
        at.this.level <- d[d==level]
        z <- complex(modulus=level,
                     argument=seq(0, 2*pi,length.out=1+length(at.this.level)))
        names(z) <- names(at.this.level)
        z <- head(z, -1)

        # Rotate origin
        # The two calls to ifelse(clockwise,...) ensure that origin.degrees 
        # works consistently regardless of the value of clockwise. 
        # -90 is so that origin.degrees of 0 puts the origin at 12 o'clock.
        z <- z * complex(argument=(ifelse(clockwise,1,-1)*origin.degrees-90 + 
                                   ifelse(clockwise,0,1)*180)*2*pi/360)

        # Adjust direction
        if(clockwise)
        {
            z <- Conj(z)
        }

        x[names(at.this.level)] <- Re(z)
        y[names(at.this.level)] <- Im(z)
    }

    if(is.null(main))
    {
        main <- paste0(CPS(community)$title, ' (', focus, ')')
    }
    PlotNPS(community, x, y, xlab='', ylab='', are.values=TRUE,
            frame.plot=frame.plot, xaxt='n', yaxt='n', xlim=lim, ylim=lim,
            main=main, ...)

}
