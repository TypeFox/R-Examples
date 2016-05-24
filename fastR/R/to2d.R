#' @export
to2d <-
function (x, y, z, type = NULL, xas = c(0.4, -0.3), yas = c(1, 
    0), zas = c(0, 1)) 
{
    if (type == "yz") {
        xas = c(0.4, -0.3)
        yas = c(1, 0)
        zas = c(0, 1)
    }
    if (type == "xy") {
        zas = c(0.4, 0.3)
        xas = c(1, 0)
        yas = c(0, 1)
    }
    if (type == "xz") {
        yas = c(0.4, 0.3)
        xas = c(1, 0)
        zas = c(0, 1)
    }
    if (length(x) != 3) {
        x <- c(x[1], y[1], z[1])
    }
    return(as.vector(x %*% rbind(xas, yas, zas)))
}
