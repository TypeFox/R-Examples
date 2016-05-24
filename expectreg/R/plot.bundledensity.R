plot.bundledensity <-
function (x, ...) 
{
    hist(x$random, breaks = seq(min(x$random), max(x$random), 
        length = 30), freq = F, xlim = range(x$x), main = "density", 
        ...)
    lines(x$x, x$density, col = "red")
}
