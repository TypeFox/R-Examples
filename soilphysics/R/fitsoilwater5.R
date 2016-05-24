fitsoilwater5 <-
function (theta, x, theta_S, xlab = NULL, ylab = NULL, ...) 
{
    if (!inherits(c(theta, x), c("numeric", "integer"))) 
        stop("non-numeric arguments!")
    if (length(theta) != length(x)) 
        stop("incompatible dimensions!")
    stopifnot(theta_S >= 0)
    dat <- data.frame(theta, x)
    if (is.null(ylab)) 
        ylab = "Soil water content"
    if (is.null(xlab)) 
        xlab = "Matric potential"
    f.graph <- function() {
        plot(theta ~ x, data = dat, las = 1, xlab = xlab, ylab = ylab, 
            main = "Soil Water Retention Curve", ...)
    }
    f.graph()
    theta_R <- alpha <- n <- b0 <- b1 <- b2 <- NULL
    f.panel <- function(pan) {
        f.graph()
        with(pan, curve(soilwater5(x, theta_R, theta_S = theta_S, alpha, 
            n, m = 1 - 1/n, b0, b1, b2), add = TRUE, col = "red"))
        return(pan)
    }
    f.fit <- function(pan) {
        start <- with(pan, pan[c("theta_R", "alpha", 
            "n", "b0", "b1", "b2")])
        fit <- try(with(pan, nls(theta ~ soilwater5(x, theta_R, 
            theta_S = theta_S, alpha, n, m = 1 - 1/n, b0, b1, b2), 
            data = dat, start = start)))
        if (inherits(fit, "try-error")) {
            rp.messagebox("No convergence... try other initial values.", 
                title = "Warning!")
        }
        else {
            f.graph()
            est <- coef(fit)
            with(dat, lines(x, soilwater5(x, theta_R = est[1], 
               theta_S = theta_S, alpha = est[2], n = est[3], b0 = est[4],
               b1 = est[5], b2 = est[6]), col = "blue"))
            print(summary(fit))
            print(Rsq(fit))
        }
        return(pan)
    }
    panel <- rp.control("Interactive fit")
    rp.slider(panel, variable = theta_R, from = 0, to = max(theta)*1.5, 
        resolution = 0.01, initval = 0.2, title = "theta_R", 
        action = f.panel)
    rp.doublebutton(panel, variable = theta_R, step = 0.01, title = "", 
        action = f.panel, showvalue = TRUE, foreground = "blue")
    rp.slider(panel, variable = alpha, from = 0, to = 2, resolution = 0.01, 
        initval = 0.05, title = "alpha", action = f.panel)
    rp.doublebutton(panel, variable = alpha, step = 0.01, title = "", 
        action = f.panel, showvalue = TRUE, foreground = "blue")
    rp.slider(panel, variable = n, from = 0, to = 30, resolution = 0.01, 
        initval = 10, title = "n", action = f.panel)
    rp.doublebutton(panel, variable = n, step = 0.01, title = "", 
        action = f.panel, showvalue = TRUE, foreground = "blue")
    rp.slider(panel, variable = b0, from = -2, to = 2, resolution = 0.01, 
        initval = 0.1, title = "b0", action = f.panel)
    rp.doublebutton(panel, variable = b0, step = 0.01, title = "", 
        action = f.panel, showvalue = TRUE, foreground = "blue")
    rp.slider(panel, variable = b1, from = -0.5, to = 0.5, resolution = 1e-04, 
        initval = -0.017, title = "b1", action = f.panel)
    rp.doublebutton(panel, variable = b1, step = 1e-04, title = "", 
        action = f.panel, showvalue = TRUE, foreground = "blue")
    rp.slider(panel, variable = b2, from = -1, to = 1, resolution = 1e-05, 
        initval = 1e-04, title = "b2", action = f.panel)
    rp.doublebutton(panel, variable = b2, step = 1e-05, title = "", 
        action = f.panel, showvalue = TRUE, foreground = "blue")
    rp.button(panel, title = "NLS estimates", action = f.fit, 
        foreground = "white", background = "navy")
    rp.button(panel, title = "__________________ Quit __________________", 
        action = function(pan) return(pan), quitbutton = TRUE, 
        foreground = "red")
}
