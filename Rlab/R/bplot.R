"bplot" <-
function (x, by,style = "tukey", outlier = TRUE, plot = TRUE, ...) 
{
    obj <- stats.bplot(x, style = style, outlier = outlier, by=by)
    if (plot) {
        bplot.obj(obj, ...)
    }
    else {
        return(obj)
    }
    invisible()
}
