fitsoilwater3 <-
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
    theta_R <- a1 <- p1 <- a2 <- p2 <- NULL
    f.panel <- function(pan) {
        f.graph()
        with(pan, curve(soilwater3(x, theta_R, a1, p1, a2, p2), 
            add = TRUE, col = "red"))
        return(pan)
    }
    f.fit <- function(pan) {
        start <- with(pan, pan[c("theta_R", "a1", "p1", "a2", 
            "p2")])
        fit <- try(with(pan, nls(theta ~ soilwater3(x, theta_R, 
            a1, p1, a2, p2), data = dat, start = start)))
        if (inherits(fit, "try-error")) {
            rp.messagebox("No convergence... try other initial values.", 
                title = "Warning!")
        }
        else {
            f.graph()
            est <- coef(fit)
            curve(soilwater3(x, est[1], est[2], est[3], est[4], 
                est[5]), add = TRUE, col = "blue")
            print(summary(fit))
            print(Rsq(fit))
        }
        return(pan)
    }
    panel <- rp.control("Interactive fit")
    ran.t <- 2 * range(theta)
    rp.slider(panel, variable = theta_R, from = 0, to = max(theta), 
        resolution = 0.01, initval = 0.8 * min(theta), title = "theta_R", 
        action = f.panel)
    rp.doublebutton(panel, variable = theta_R, step = 0.01, title = "", 
        action = f.panel, showvalue = TRUE, foreground = "blue")
    rp.slider(panel, variable = a1, from = -0.5, to = 10, resolution = 0.01, 
        initval = 0.07, title = "a1", action = f.panel)
    rp.doublebutton(panel, variable = a1, step = 0.01, title = "", 
        action = f.panel, showvalue = TRUE, foreground = "blue")
    rp.slider(panel, variable = p1, from = 0, to = 15000, resolution = 5, 
        initval = 3670, title = "p1", action = f.panel)
    rp.doublebutton(panel, variable = p1, step = 1, title = "", 
        action = f.panel, showvalue = TRUE, foreground = "blue")
    rp.slider(panel, variable = a2, from = 0, to = 10, resolution = 0.01, 
        initval = 0.32, title = "a2", action = f.panel)
    rp.doublebutton(panel, variable = a2, step = 0.01, title = "", 
        action = f.panel, showvalue = TRUE, foreground = "blue")
    rp.slider(panel, variable = p2, from = 0, to = 1500, resolution = 5, 
        initval = 70, title = "p2", action = f.panel)
    rp.doublebutton(panel, variable = p2, step = 1, title = "", 
        action = f.panel, showvalue = TRUE, foreground = "blue")
    rp.button(panel, title = "NLS estimates", action = f.fit, 
        foreground = "white", background = "navy")
    rp.button(panel, title = "__________________ Quit __________________", 
        action = function(pan) return(pan), quitbutton = TRUE, 
        foreground = "red")
}
