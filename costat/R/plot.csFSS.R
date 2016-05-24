plot.csFSS <-
function (x, solno = NULL, ALLplotclust = TRUE, ALLplotscale = TRUE,
    sollabels=TRUE, 
    SNinc = 0, ...) 
{
    if (is.null(solno)) 
        plot(COEFbothscale(x), plotclust = ALLplotclust, plotscale = ALLplotscale, sollabels=sollabels, 
            ...)
    else {
        LCTSres(x, tsx = x$tsx, tsy = x$tsy, inc = SNinc, solno = solno, 
            ...)
    }
}
