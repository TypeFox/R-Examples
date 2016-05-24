tryCatch(utils::globalVariables(c('gp')),
         error=function(e) message('Looks like you should update R.'))

#' Geometric representation of linear model
#' 
#' \code{geolm} create a graphical representation of the fit of a linear model.
#' 
#' 
#' @aliases geolm to2d
#' @param formula a formula as used in \code{\link{lm}}.
#' @param data a data frame as in \code{\link{lm}}.
#' @param type character: indicating the type of projection to use to collapse
#' multi-dimensional data space into two dimensions of the display.
#' @param version an integer (currently \code{1} or \code{2}).  Which version
#' of the plot to display.
#' @param plot a logical: should the plot be displayed?
#' @param x,y,z numeric.
#' @param xas,yas,zas numeric vector of length 2 indicating the projection of
#' \code{c(1,0,0)}, \code{c(0,1,0)}, and \code{c(0,0,1)}.
#' @param \dots other arguments passed to \code{\link{lm}}
#' @author Randall Pruim
#' @seealso \code{\link{lm}}.
#' @keywords stats
#' @export
#' @examples
#' 
#' geolm(pollution ~ location, data=airpollution)
#' geolm(distance ~ projectileWt, data=trebuchet2)
#' 
geolm <-
function (formula, data = parent.env(), type = "xz", version = 1, 
    plot = TRUE, ...) 
{
    model <- lm(formula, data = data, ...)
    ef <- -1 * effects(model)
    w.int <- which(names(ef) == "(Intercept)")
    w.resid <- which(names(ef) == "")
    w.effect <- (1:length(ef))[-1 * c(w.int, w.resid)]
    s.int <- sum(ef[w.int]^2)
    s.effect <- sum(ef[w.effect]^2)
    s.resid <- sum(ef[w.resid]^2)
    l.int <- sign(ef[w.int]) * sqrt(s.int)
    l.effect <- sqrt(s.effect)
    l.resid <- sqrt(s.resid)
    if (version == 1) {
        transform <- function(x) {
            to2d(x, type = type)
        }
    }
    else {
        transform <- function(x) {
            x
        }
    }
    origin <- transform(c(0, 0, 0))
    y_bar_pt <- transform(c(l.int, 0, 0))
    effect_vec <- transform(c(0, l.effect, 0))
    y_hat_pt <- y_bar_pt + effect_vec
    resid_vec <- transform(c(0, 0, l.resid))
    y_pt <- y_hat_pt + resid_vec
    pts <- rbind(origin, y_bar_pt, origin, y_hat_pt, origin, 
        y_pt, y_bar_pt, y_hat_pt, y_bar_pt, y_pt, y_hat_pt, y_pt)
    if (version == 1) {
        pts <- cbind(pts, rep(0, nrow(pts)))
    }
    d <- data.frame(x = pts[, 1], y = pts[, 2], z = pts[, 3], gp = rep(letters[1:6], each = 2))
    aspect <- c(diff(range(d$y))/diff(range(d$x)), diff(range(d$z))/diff(range(d$x)))
    if (plot) {
        if (version == 1) {
          return(xyplot(y ~ x, data = d, groups = gp, type = "l", 
                        main=deparse(formula),
                        lwd = 4, lty = 1, aspect = "iso", 
                        col = c("orange", "blue", "black", "forestgreen", "purple", "red"), 
                        scales = list(draw = F), xlab = "", ylab = "")
                 )
        }
        else {
          return(cloud(z ~ x * y, data = d, groups = gp, aspect = aspect, 
                       main=deparse(formula),
                       col = c("orange", "blue", "black", "forestgreen",  "purple", "red"), 
                       xlab = "mean", ylab = "effect", 
                       zlab = "resid", scales = list(draw = F), type = "l")
          )
        }
    }
    return(model)
}


#' @rdname geolm
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
