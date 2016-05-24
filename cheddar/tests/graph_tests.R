.TestPlotFunction <- function(plot.fn, ...)
{
    # Calls fn with args, redirecting graphics output to a temporary file that
    # is unlinked after the test
    temp.path <- tempfile()
    png(temp.path)
    print(temp.path)
    tryCatch(do.call(plot.fn, args=list(...)),
        finally={dev.off(); unlink(temp.path)})
}

TestSimplePlotFunctions <- function()
{
    # Test those plot functions that take a community object as the first, and 
    # only mandatory, argument
    exclude <- c('PlotNPSDistribution', 'PlotRankNPS', 'PlotTLPS', 'PlotNPS',
                 'PlotWagonWheel', 'PlotLinearModels')
    for(fn in setdiff(ls('package:cheddar', pattern='Plot*'), exclude))
    {
        cat(fn, '\n')
        .TestPlotFunction(fn, TL84)
    }
}

TestPlotWagonWheel <- function()
{
    .TestPlotFunction(PlotWagonWheel, TL84, focus=30)
}
