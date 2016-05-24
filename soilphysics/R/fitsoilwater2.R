fitsoilwater2 <-
function (theta, x, x0 = 6.653, xlab = NULL, ylab = NULL, ...) 
{
    if (!inherits(c(theta, x), "numeric")) 
        stop("non-numeric arguments!")
    if (length(theta) != length(x)) 
        stop("incompatible dimensions!")
    dat <- data.frame(theta, x)
    if (is.null(xlab)) 
        xlab = "pF"
    if (is.null(ylab)) 
        ylab = "Soil water content"
    f.graph <- function() {
        plot(theta ~ x, data = dat, las = 1, xlab = xlab, ylab = ylab, 
            main = "Soil Water Retention Curve", ...)
    }
    f.graph()
    k0 <- k1 <- n <- NULL
    f.panel <- function(pan) {
        f.graph()
        with(pan, curve(soilwater2(x, x0, k0, k1, n), add = TRUE, 
            col = "red"))
        return(pan)
    }
    f.fit <- function(pan) {
        start <- with(pan, pan[c("k0", "k1", "n")])
        fit <- try(with(pan, nls(theta ~ soilwater2(x, x0, k0, 
            k1, n), data = dat, start = start)))
        if (inherits(fit, "try-error")) {
            rp.messagebox("No convergence... try other initial values.", 
                title = "Warning!")
        }
        else {
            f.graph()
            est <- coef(fit)
            curve(soilwater2(x, x0, est[1], est[2], est[3]), 
                add = TRUE, col = "blue")
            print(summary(fit))
            print(Rsq(fit))
        }
        return(pan)
    }
    panel <- rp.control("Interactive fit")
    ran.t <- 2 * range(theta)
    rp.slider(panel, variable = k0, from = 0, to = 10, resolution = 0.01, 
        initval = 1.9, title = "k0", action = f.panel)
    rp.doublebutton(panel, variable = k0, step = 0.01, title = "", 
        action = f.panel, showvalue = TRUE, foreground = "blue")
    rp.slider(panel, variable = k1, from = 0, to = 10, resolution = 0.01, 
        initval = 0.4, title = "k1", action = f.panel)
    rp.doublebutton(panel, variable = k1, step = 0.01, title = "", 
        action = f.panel, showvalue = TRUE, foreground = "blue")
    rp.slider(panel, variable = n, from = 0, to = 10, resolution = 0.01, 
        initval = 2.4, title = "n", action = f.panel)
    rp.doublebutton(panel, variable = n, step = 0.01, title = "", 
        action = f.panel, showvalue = TRUE, foreground = "blue")
    rp.button(panel, title = "NLS estimates", action = f.fit,
        foreground = "white", background = "navy")
    rp.button(panel, title = "__________________ Quit __________________", 
        action = function(pan) return(pan), quitbutton = TRUE, 
        foreground = "red")
}
