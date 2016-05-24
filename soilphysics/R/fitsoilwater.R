fitsoilwater <-
function (theta, x, xlab = NULL, ylab = NULL, ...) 
{
    if (!inherits(c(theta, x), "numeric")) 
        stop("non-numeric arguments!")
    if (length(theta) != length(x)) 
        stop("incompatible dimensions!")
    dat <- data.frame(theta, x)
    if (is.null(xlab)) 
        xlab = "Matric potential"
    if (is.null(ylab)) 
        ylab = "Soil water content"
    f.graph <- function() {
        plot(theta ~ x, data = dat, las = 1, xlab = xlab, ylab = ylab, 
            main = "Soil Water Retention Curve", ...)
    }
    f.graph()
    theta_R <- theta_S <- alpha <- n <- NULL
    f.panel <- function(pan) {
        f.graph()
        with(pan, curve(soilwater(x, theta_R, theta_S, alpha, 
            n), add = TRUE, col = "red"))
        return(pan)
    }
    f.fit <- function(pan) {
        start <- with(pan, pan[c("theta_R", "theta_S", "alpha", 
            "n")])
        fit <- try(with(pan, nls(theta ~ soilwater(x, theta_R, 
            theta_S, alpha, n), data = dat, start = start)))
        if (inherits(fit, "try-error")) {
            rp.messagebox("No convergence... try other initial values.", 
                title = "Warning!")
        }
        else {
            f.graph()
            est <- coef(fit)
            curve(soilwater(x, est[1], est[2], est[3], est[4]), 
                add = TRUE, col = "blue")
            print(summary(fit))
            print(Rsq(fit))
        }
        return(pan)
    }
    panel <- rp.control("Interactive fit")
    rp.slider(panel, variable = theta_R, from = 0, to = max(theta), 
        resolution = 0.01, initval = 0.8 * min(theta), 
        title = "theta_R", action = f.panel)
    rp.doublebutton(panel, variable = theta_R, step = 0.01, title = "", 
        action = f.panel, showvalue = TRUE, foreground = "blue")
    rp.slider(panel, variable = theta_S, from = 0, to = max(theta), 
        resolution = 0.01, initval = 0.8 * max(theta), 
        title = "theta_S", action = f.panel)
    rp.doublebutton(panel, variable = theta_S, step = 0.01, title = "", 
        action = f.panel, showvalue = TRUE, foreground = "blue")
    rp.slider(panel, variable = alpha, from = 0, to = 2, resolution = 0.01, 
        initval = 0.01, title = "alpha", action = f.panel)
    rp.doublebutton(panel, variable = alpha, step = 0.01, title = "", 
        action = f.panel, showvalue = TRUE, foreground = "blue")
    rp.slider(panel, variable = n, from = 0, to = 15, resolution = 0.01, 
        initval = 2, title = "n", action = f.panel)
    rp.doublebutton(panel, variable = n, step = 0.01, title = "", 
        action = f.panel, showvalue = TRUE, foreground = "blue")
    rp.button(panel, title = "NLS estimates", action = f.fit,
        foreground = "white", background = "navy")
    rp.button(panel, title = "__________________ Quit __________________", 
        action = function(pan) return(pan), quitbutton = TRUE, 
        foreground = "red")
}
