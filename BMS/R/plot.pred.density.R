plot.pred.density <-
function (x, predict_index = NULL, addons = "eslz", realized.y = NULL, 
    addons.lwd = 1.5, ...) 
{
    if (!is(x, "pred.density")) 
        stop("x must be of class 'pred.density'!")
    x$plot(predict_index, realized.y = realized.y, addons = addons, 
        addons.lwd = addons.lwd, ...)
}
