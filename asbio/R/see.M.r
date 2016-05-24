see.M <- function () 
{

    if (!exists("slider.env")) 
        slider.env <- NULL
    suppressWarnings(rm(slider.env))
    slider.env <<- new.env()
    cval <- 0
    assign("cval", tclVar(cval), envir = slider.env)
    xmin <- -5
    assign("xmin", tclVar(xmin), envir = slider.env)
    xmax <- 5
    assign("xmax", tclVar(xmax), envir = slider.env)
    norm.refresh <- function(...) {
        cval <- as.numeric(evalq(tclvalue(cval), envir = slider.env))
        xmin <- as.numeric(evalq(tclvalue(xmin), envir = slider.env))
        xmax <- as.numeric(evalq(tclvalue(xmax), envir = slider.env))
        xx <- seq(xmin, xmax, length = 500)
        yy <- sapply(xx, function(x) {
            min(1, cval/abs(x))
        })
        dev.hold()
        p <- plot(xx, yy, type = "l", xlim = c(xmin, xmax), ylab = "Weight", 
            xlab = expression(paste(italic(x), " - ", hat(mu))), 
            lwd = 2, bty = "l")
        grid(p)
        legend("topright", legend = bquote(paste(italic(c)," = ", .(cval), sep = "")), box.col = "white", 
            bg = "white")
        dev.flush()
    }
    tclServiceMode(TRUE)
    m <- tktoplevel()
    tkwm.title(m, "Visualizing M-estimation of location")
    tkwm.geometry(m, "+0+0")
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "c", width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 0.01, 
        to = 5, orient = "horiz", resolution = 0.01, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = cval), envir = slider.env)
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "Xmin:", width = 6), side = "left")
    tkpack(e <- tkentry(fr, width = 8), side = "left")
    assign("e", e, envir = slider.env)
    evalq(tkconfigure(e, textvariable = xmin), envir = slider.env)
    tkpack(tklabel(fr, text = "Xmax:", width = 6), side = "left")
    tkpack(e <- tkentry(fr, width = 8), side = "left")
    assign("e", e, envir = slider.env)
    evalq(tkconfigure(e, textvariable = xmax), envir = slider.env)
    tkpack(tkbutton(m, text = "Refresh", command = norm.refresh), 
        side = "left")
    tkpack(tkbutton(m, text = "Exit", command = function() tkdestroy(m)), 
        side = "right")
}
