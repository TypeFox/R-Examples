# Miscelaneous plot functions
.PlotPyramid <- function(values, xlim, xlab, ylab, col, main, 
                         show.level.labels, text.col, ...)
{
    # Plots a pyramid of values
    if(is.null(xlim))
    {
        xlim <- range(values, na.rm=TRUE)
    }

    # Values can be -ve so offset them
    plot.values <- values - xlim[1] + 1
    stopifnot(min(plot.values, na.rm=TRUE)>=0)

    # The largest value that we can show
    max.lxl <- xlim[2] - xlim[1] + 1

    if(!is.null(names(col)))
    {
        col <- col[names(values)]
    }

    plot(NA, NA, type='n', xlim=c(-max.lxl/2, max.lxl/2), xaxt='n', xlab=xlab, 
         ylim=c(1,1+length(plot.values)), yaxt='n', ylab=ylab, main=main, 
         frame.plot=FALSE, ...)

    rect(-abs(plot.values)/2,    1:length(plot.values), 
          abs(plot.values)/2, 1+(1:length(plot.values)), 
         col=col)

    to.print <- sprintf('%.2f', values)
    to.print[is.na(values)] <- ''
    text(0, y=1:length(values)+0.5, to.print, col=text.col, ...)
    if(show.level.labels)
    {
        axis(2, at=1:length(plot.values)+0.5, labels=names(plot.values), las=1, 
             tick=FALSE, ...)
    }
}

.PyramidLevels <- function(values, level, expected, fill.missing, 
                           order.by.expected)
{
    # Some checks and processing on values used for levels in pyramid plots
    if(missing(expected))
    {
        if(is.numeric(level))
        {
            expected <- as.character(floor(min(level)):ceiling(max(level)))
        }
        else if('category'==level)
        {
            expected <- c(.UnnamedString(),'producer','invertebrate',
                          'vert.ecto','vert.endo')
            expected <- intersect(expected, names(values))
        }
        else
        {
            expected <- names(values)
        }
    }

    duplicated.levels <- duplicated(expected)
    extra.levels <- setdiff(names(values), expected)
    if(any(duplicated.levels))
    {
        stop(paste('The levels [', 
                   paste(expected[duplicated.levels], collapse=','), 
                   '] appear in expected.levels more than once', sep=''))
    }
    else if(length(extra.levels)>0)
    {
        stop(paste('The levels [', paste(extra.levels, collapse=','), 
                   '] are not in expected.levels', sep=''))
    }
    else
    {
        if(fill.missing)
        {
            missing.levels <- setdiff(expected, names(values))
            values[missing.levels] <- NA
        }
        if(order.by.expected)
        {
            values <- values[intersect(expected, names(values))]
        }
    }
    return (values)
}

PlotNPyramid <- function(community, 
                         level=floor(PreyAveragedTrophicLevel(community)),
                         expected.levels,
                         fill.missing.levels=TRUE,
                         order.by.expected=TRUE,
                         show.level.labels=TRUE,
                         xlab=Log10NLabel(community, 
                                          name=expression(~sum(italic(N)))), 
                         ylab='', 
                         xlim=NULL, 
                         col=NULL, 
                         text.col='black',
                         main=CPS(community)$title, 
                         ...)
{
    if(!is.Community(community)) stop('Not a Community')

    .RequireN(community)
    values <- log10(SumNByClass(community, class=level, na.rm=TRUE))
    values[is.infinite(values)] <- NA
    values <- .PyramidLevels(values, level, expected.levels, 
                             fill.missing.levels, order.by.expected)
    .PlotPyramid(values=values, xlab=xlab, ylab=ylab, xlim=xlim, 
                 col=col, main=main, show.level.labels=show.level.labels, 
                 text.col=text.col, ...)
}

PlotBPyramid <- function(community,
                         level=floor(PreyAveragedTrophicLevel(community)),
                         expected.levels,
                         fill.missing.levels=TRUE,
                         order.by.expected=TRUE,
                         show.level.labels=TRUE,
                         xlab=Log10BLabel(community, 
                                          name=expression(~sum(italic(B)))), 
                         ylab='', 
                         xlim=NULL, 
                         col=NULL, 
                         text.col='black',
                         main=CPS(community)$title, 
                         ...)
{
    if(!is.Community(community)) stop('Not a Community')

    .RequireM(community)
    .RequireN(community)
    values <- log10(SumBiomassByClass(community, class=level, na.rm=TRUE))
    values[is.infinite(values)] <- NA
    values <- .PyramidLevels(values, level, expected.levels, 
                             fill.missing.levels, order.by.expected)
    .PlotPyramid(values=values, xlab=xlab, ylab=ylab, xlim=xlim, 
                 col=col, main=main, show.level.labels=show.level.labels, 
                 text.col=text.col, ...)
}

.PlotAbundanceSpectrum <- function(bins, binned.data, main, 
                                   xlab, ylab, xlim, ylim, pch, 
                                   show.bin.limits, show.bin.centres, ...)
{
    if(is.null(ylim))
    {
        ylim <- range(binned.data)
    }

    bin.centres <- attr(bins, 'bin.centres')
    breaks <- attr(bins, 'breaks')

    if(is.null(xlim))
    {
        xlim <- c(breaks[1], tail(breaks, 1))
    }

    plot(bin.centres[as.integer(names(binned.data))], binned.data, 
         xlab=xlab, xlim=xlim, ylab=ylab, ylim=ylim, 
         pch=pch, main=main, ..., xaxt='n')

    # Tick marks at top and right of plot
    axis(1, at=bin.centres, labels=round(bin.centres, 2), ...)

    dots <- list(...)
    if('n'!=par('xaxt') && !'n' %in% dots[['xaxt']])
    {
        axis(3, at=bin.centres, labels=FALSE, ...)
    }

    if('n'!=par('yaxt') && !'n' %in% dots[['yaxt']])
    {
        axis(4, labels=FALSE, ...)
    }

    if(show.bin.limits)
    {
        # Vertical lines delimiting bins
        abline(v=breaks, lty=3)
    }

    if(show.bin.centres)
    {
        # Vertical lines showing bin centres
        abline(v=bin.centres, lty=2)
    }

    # Plot lm
    x <- bin.centres[as.integer(names(binned.data))]
    y <- binned.data
    m <- lm(y~x)
    abline(m)

    to.return <- list(bins=bins, lm=m)

    invisible(to.return)
}

PlotNSpectrum <- function(community, 
                          lower=min(NP(community, 'M'), na.rm=TRUE), 
                          upper=max(NP(community, 'M'), na.rm=TRUE), 
                          n.bins=10, 
                          main=CPS(community)$title,  
                          xlab=Log10MLabel(community), 
                          ylab=Log10NLabel(community), 
                          xlim=NULL,
                          ylim=NULL, 
                          pch=19, 
                          show.bin.limits=TRUE, 
                          show.bin.centres=FALSE, 
                          ...)
{
    # The log10(sum(numerical abundance)) in equally-spaced log10(M) bins, 
    # together with a linear regression. If not NULL N.final is also binned 
    # and plotted together with it's own linear regression. The function 
    # returns a list(bin.centres, lm).

    if(!is.Community(community)) stop('Not a Community')

    .RequireM(community)
    .RequireN(community)

    # Bin totals
    bins <- BodyMassBins(community, lower=lower, upper=upper, n.bins=n.bins)
    binned.data <- log10(SumNByClass(community, class=bins))

    return (.PlotAbundanceSpectrum(bins, binned.data, main, xlab, ylab, 
                                   xlim, ylim, pch, show.bin.limits, 
                                   show.bin.centres, ...))
}

PlotBSpectrum <- function(community, 
                          lower=min(NP(community, 'M'), na.rm=TRUE), 
                          upper=max(NP(community, 'M'), na.rm=TRUE), 
                          n.bins=10, 
                          main=CPS(community)$title,
                          xlab=Log10MLabel(community), 
                          ylab=Log10BLabel(community), 
                          xlim=NULL,
                          ylim=NULL, 
                          pch=19, 
                          show.bin.limits=TRUE, 
                          show.bin.centres=FALSE, 
                          ...)
{
    # The log10(sum(biomass abundance)) in equally-spaced log10(M) bins, 
    # together with a linear regression. Returns a list(bin.centres, lm).

    if(!is.Community(community)) stop('Not a Community')

    .RequireM(community)
    .RequireN(community)

    # Bin totals
    bins <- BodyMassBins(community, lower=lower, upper=upper, n.bins=n.bins)
    binned.data <- log10(SumBiomassByClass(community, class=bins))

    return (.PlotAbundanceSpectrum(bins, binned.data, main, xlab, ylab, 
                                   xlim, ylim, pch, show.bin.limits, 
                                   show.bin.centres, ...))
}

PlotAuppervAlower <- function(community, 
                              main=CPS(community)$title,
                              xlab=~A[lower], 
                              ylab=~A[upper], 
                              xlim=c(-180, 180),
                              ylim=c(-180, 180), 
                              pch=19, 
                              ...)
{
    # Upper vs lower link angles. Cohen et al 2009 PNAS.

    # Ensure that we pick up the title of the unmanipulated community, when 
    # using the default value of main.
    main <- main

    community <- RemoveNodes(community, remove=with(NPS(community), 
                                                    node[is.na(M) | is.na(N)]))
    community <- RemoveCannibalisticLinks(community)
    tncp <- ThreeNodeChains(community, 
                            chain.properties='.NvMThreeNodeChainProperties')

    plot(tncp[,'Alower'], tncp[,'Aupper'], pch=pch, xlim=xlim, 
         ylim=ylim, main=main, xlab=xlab, ylab=ylab, ...)

    axis(3, labels=FALSE)
    axis(4, labels=FALSE)

    abline(v=median(tncp[,'Alower']))
    abline(h=median(tncp[,'Aupper']))
}

