#' @export
panel.xyplotsmooth <-
function (x, y, type = c("p", "smooth"), ...) 
{
    panel.xyplot(x, y, type = c("smooth"), ...)
    panel.xyplot(x, y, type = c("p"), ...)
}
