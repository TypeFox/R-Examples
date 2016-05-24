#' @export
panel.lm <-
function (x, y, model, fits, ...) 
{
    if (missing(fits)) {
        if (missing(model)) {
            model <- lm(y ~ x, ...)
        }
        fits <- fitted(model)
    }
    panel.xyplot(x[order(x)], fits[order(x)], type = "l", lwd = 2, 
        ...)
    panel.xyplot(x, y, ...)
}
