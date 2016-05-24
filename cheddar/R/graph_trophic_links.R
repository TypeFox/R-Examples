# High-level and helper functions for plotting trophic links
PlotTLPS <- function(community, 
                     X, 
                     Y, 
                     xlab, 
                     ylab, 
                     axes.limits.equal=FALSE,
                     xlim=NULL, 
                     ylim=NULL, 
                     main=CPS(community)$title, 
                     highlight.links=NULL,
                     lowlight.links=NULL,
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
                     are.values=FALSE, 
                     ...)
{
    if(!are.values)
    {
        if(missing(xlab))   xlab <- .ResolveTrophicLinkAxisLabel(X)
        if(missing(ylab))   ylab <- .ResolveTrophicLinkAxisLabel(Y)

        X <- .ResolveTrophicLinkProperty(community, X)
        Y <- .ResolveTrophicLinkProperty(community, Y)
    }
    else
    {
        if(missing(xlab))   xlab <- ''
        if(missing(ylab))   ylab <- ''

        stopifnot(length(X)==NumberOfTrophicLinks(community))
        stopifnot(length(Y)==NumberOfTrophicLinks(community))
    }

    if(missing(colour.by))
    {
        colour.by <- NULL
        if('category' %in% NodePropertyNames(community))
        {
            colour.by <- 'resource.category'
        }
    }

    if(missing(colour.spec))
    {
        colour.spec <- NULL
        if(!is.null(colour.by) && 
           colour.by %in% c('resource.category','consumer.category'))
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
            bg.by <- 'resource.category'
        }
    }

    if(missing(bg.spec))
    {
        bg.spec <- NULL
        if(!is.null(bg.by) && 
           bg.by %in% c('resource.category','consumer.category'))
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
            symbol.by <- 'resource.category'
        }
    }

    if(missing(symbol.spec))
    {
        symbol.spec <- NULL
        if(!is.null(symbol.by) && 
           symbol.by %in% c('resource.category','consumer.category'))
        {
            symbol.spec <- DefaultCategorySymbols()
        }
    }

    if(is.function(highlight.links))
    {
        highlight.links <- highlight.links(community)
    }
    highlight.links <- .ResolveTrophicLinksToRowIndices(community, 
                                                        highlight.links)

    if(is.function(lowlight.links))
    {
        lowlight.links <- lowlight.links(community)
    }
    lowlight.links <- .ResolveTrophicLinksToRowIndices(community,lowlight.links)


    col <- .TrophicLinkGraphParamFromSpec(community, 'col', col, 
                                          colour.by, colour.spec, 
                                          default=par('col'))
    col <- .HighlightColours(NumberOfTrophicLinks(community), col, 
                             highlight.links, lowlight.links)

    pch <- .TrophicLinkGraphParamFromSpec(community, 'pch', pch, 
                                          symbol.by, symbol.spec, default=19)
    pch <- .HighlightSymbols(NumberOfTrophicLinks(community), pch, 
                             highlight.links, lowlight.links)

    bg <- .TrophicLinkGraphParamFromSpec(community, 'bg', bg, 
                                         bg.by, bg.spec)
    bg <- .HighlightBG(NumberOfTrophicLinks(community), bg, highlight.links, 
                       lowlight.links)

    cex <- .TrophicLinkGraphParamFromSpec(community, 'cex', cex, 
                                          cex.by, cex.spec)
    cex <- .HighlightCex(NumberOfTrophicLinks(community), cex, highlight.links, 
                         lowlight.links)

    if(axes.limits.equal && 
       (is.null(xlim) || missing(xlim)) && (is.null(ylim) || missing(ylim)))
    {
        xlim <- ylim <- range(c(X, Y))
    }

    plot(X, Y, xlim=xlim, ylim=ylim, main=main, xlab=xlab, ylab=ylab, type='n', 
         ...)
    .AddAxisTicks(...)
    .PlotHighlightedPoints(X, Y, col=col, pch=pch, bg=bg, cex=cex, 
                           labels=NULL, label.col=NULL, label.cex=NULL, 
                           show.points=TRUE, show.labels=FALSE, 
                           highlight=highlight.links, ...)
}

PlotPredationMatrix <- function(community, 
                                xlab='Consumer', 
                                ylab='Resource', 
                                resource.order, 
                                consumer.order,
                                ...)
{
    # Plots the predation matrix
    if(!is.Community(community)) stop('Not a Community')

    .RequireTrophicLinks(community)

    tlp <- TLPS(community)
    nodes <- unname(NP(community, 'node'))

    if(missing(resource.order))
    {
        resource.order <- nodes
    }
    else if(1==length(resource.order))
    {
        resource.order <- order(NPS(community, resource.order)[,resource.order])
        resource.order <- nodes[resource.order]
    }
    else if(is.numeric(resource.order))
    {
        stopifnot(all(0 < resource.order & 
                      resource.order <= NumberOfNodes(community)))
        stopifnot(length(unique(resource.order)) == NumberOfNodes(community))
        stopifnot(range(resource.order) == c(1, NumberOfNodes(community)))
        resource.order <- nodes[resource.order]
    }

    if(missing(consumer.order))
    {
        consumer.order <- nodes
    }
    else if(1==length(consumer.order))
    {
        consumer.order <- order(NPS(community, consumer.order)[,consumer.order])
        consumer.order <- nodes[consumer.order]
    }
    else if(is.numeric(consumer.order))
    {
        stopifnot(all(0 < consumer.order & 
                      consumer.order <= NumberOfNodes(community)))
        stopifnot(length(unique(consumer.order)) == NumberOfNodes(community))
        stopifnot(range(consumer.order) == c(1, NumberOfNodes(community)))
        consumer.order <- nodes[consumer.order]
    }

    stopifnot(sort(resource.order)==sort(nodes))
    stopifnot(sort(consumer.order)==sort(nodes))

    decode <- 1:length(nodes)
    names(decode) <- consumer.order
    X <- decode[tlp[,'consumer']]

    names(decode) <- resource.order
    Y <- decode[tlp[,'resource']]

    # Plot y values inverted
    n <- length(nodes)
    PlotTLPS(community, 
             X=X, 
             Y=1+n-Y,
             xlab=xlab, ylab=ylab, 
             xlim=c(1, n), ylim=c(1, n), 
             are.values=TRUE, 
             xaxt='n', yaxt='n', ...)

    if(all(resource.order==consumer.order))
    {
        # Points on this line are cannibals
        abline(a=n+1, b=-1, lty=2)
    }
}

PlotMRvMC <- function(community, 
                      xlab=Log10MLabel(community, name='italic(M)[consumer]'),
                      ylab=Log10MLabel(community, name='italic(M)[resource]'),
                      axes.limits.equal=TRUE,
                      ...)
{
    # Plots log10M of resources against log10M of consumers
    if(!is.Community(community)) stop('Not a Community')
    .RequireM(community)
    .RequireTrophicLinks(community)
    PlotTLPS(community, X='consumer.Log10M', 
             Y='resource.Log10M', xlab=xlab, ylab=ylab, 
             are.values=FALSE, axes.limits.equal=axes.limits.equal, ...)
    abline(a=0, b=1, lty=2)
}

PlotMCvMR <- function(community, 
                      xlab=Log10MLabel(community, name='italic(M)[resource]'),
                      ylab=Log10MLabel(community, name='italic(M)[consumer]'),
                      axes.limits.equal=TRUE,
                      ...)
{
    # Plots log10M of resources against log10M of consumers
    if(!is.Community(community)) stop('Not a Community')
    .RequireM(community)
    .RequireTrophicLinks(community)
    PlotTLPS(community, X='resource.Log10M', 
             Y='consumer.Log10M', xlab=xlab, ylab=ylab, 
             are.values=FALSE, axes.limits.equal=axes.limits.equal, ...)
    abline(a=0, b=1, lty=2)
}

PlotNRvNC <- function(community, 
                      xlab=Log10NLabel(community, name='italic(N)[consumer]'), 
                      ylab=Log10NLabel(community, name='italic(N)[resource]'), 
                      axes.limits.equal=TRUE,
                      ...)
{
    # Plots log10M of resources against log10M of consumers
    if(!is.Community(community)) stop('Not a Community')
    .RequireN(community)
    .RequireTrophicLinks(community)
    PlotTLPS(community, X='consumer.Log10N', 
             Y='resource.Log10N', xlab=xlab, ylab=ylab, 
             are.values=FALSE, axes.limits.equal=axes.limits.equal, ...)
    abline(a=0, b=1, lty=2)
}

PlotNCvNR <- function(community, 
                      xlab=Log10NLabel(community, name='italic(N)[resource]'), 
                      ylab=Log10NLabel(community, name='italic(N)[consumer]'), 
                      axes.limits.equal=TRUE,
                      ...)
{
    # Plots log10M of resources against log10M of consumers
    if(!is.Community(community)) stop('Not a Community')
    .RequireN(community)
    .RequireTrophicLinks(community)
    PlotTLPS(community, X='resource.Log10N', 
             Y='consumer.Log10N', xlab=xlab, ylab=ylab, 
             are.values=FALSE, axes.limits.equal=axes.limits.equal, ...)
    abline(a=0, b=1, lty=2)
}

PlotBRvBC <- function(community, 
                      xlab=Log10BLabel(community, name='italic(B)[consumer]'),
                      ylab=Log10BLabel(community, name='italic(B)[resource]'),
                      axes.limits.equal=TRUE,
                      ...)
{
    # Plots log10M of resources against log10M of consumers
    if(!is.Community(community)) stop('Not a Community')
    .RequireM(community)
    .RequireTrophicLinks(community)
    PlotTLPS(community, X='consumer.Log10Biomass', 
             Y='resource.Log10Biomass', xlab=xlab, ylab=ylab, 
             are.values=FALSE, axes.limits.equal=axes.limits.equal, ...)
    abline(a=0, b=1, lty=2)
}

PlotBCvBR <- function(community, 
                      xlab=Log10BLabel(community, name='italic(B)[resource]'),
                      ylab=Log10BLabel(community, name='italic(B)[consumer]'),
                      axes.limits.equal=TRUE,
                      ...)
{
    # Plots log10M of resources against log10M of consumers
    if(!is.Community(community)) stop('Not a Community')
    .RequireM(community)
    .RequireTrophicLinks(community)
    PlotTLPS(community, X='resource.Log10Biomass', 
             Y='consumer.Log10Biomass', xlab=xlab, ylab=ylab, 
             are.values=FALSE, axes.limits.equal=axes.limits.equal, ...)
}

PlotDegreeDistribution <- function(community, 
                                   xlab='Number of links', 
                                   ...)
{
    if(!is.Community(community)) stop('Not a Community')
    .RequireTrophicLinks(community)
    d <- DegreeDistribution(community)
    barplot(height=d, xlab=xlab, ...)
}

