see.nbin.tck <-function () 
{

    if (!exists("slider.env")) 
        slider.env <- NULL
    suppressWarnings(rm(slider.env))
    slider.env <<- new.env()
    r <- 1
    p <- 0.5
    assign("r", tclVar(r), envir = slider.env)
    assign("p", tclVar(p), envir = slider.env)
    xmin <- 0
    assign("xmin", tclVar(xmin), envir = slider.env)
    xmax <- 20
    assign("xmax", tclVar(xmax), envir = slider.env)
    norm.refresh <- function(...) {
        r <- as.numeric(evalq(tclvalue(r), envir = slider.env))
        p <- as.numeric(evalq(tclvalue(p), envir = slider.env))
        xmin <- as.numeric(evalq(tclvalue(xmin), envir = slider.env))
        xmax <- as.numeric(evalq(tclvalue(xmax), envir = slider.env))
        xx <- seq(xmin, xmax)
        yy <- dnbinom(xx, r, p)
        plot(xx, yy, type = "h", xlab = expression(italic(x)), 
            ylab = expression(paste(italic(f), "(", italic(x), 
                ")", sep = "")))
        points(xx, yy, pch = 19)
        mtext(bquote(paste(italic(X), " ~ ", italic(NB), "(", 
            .(r), ", ", .(p), ")", sep = "")), line = 1, side = 3)
        dev.flush()
    }
    tclServiceMode(TRUE)
    m <- tktoplevel()
    tkwm.title(m, "NB(r, \u03C0)")
    tkpack(tklabel(m, text = "      Visualizing the Negative Binomial Distribution      "))
    tkwm.geometry(m, "+0+0")
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "r", font = c("Helvetica", "9", 
        "italic"), width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 1, 
        to = 20, orient = "horiz", resolution = 1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = r), envir = slider.env)
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "\u03C0", font = c("Helvetica", "9", 
        "italic"), width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 0.01, 
        to = 1, orient = "horiz", resolution = 0.01, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = p), envir = slider.env)
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








see.nbincdf.tck <- function () 
{

    if (!exists("slider.env")) 
        slider.env <- NULL
    suppressWarnings(rm(slider.env))
    slider.env <<- new.env()
    r <- 1
    p <- 0.5
    assign("r", tclVar(r), envir = slider.env)
    assign("p", tclVar(p), envir = slider.env)
    xmin <- 0
    assign("xmin", tclVar(xmin), envir = slider.env)
    xmax <- 20
    assign("xmax", tclVar(xmax), envir = slider.env)
    dev.new(height = 4, width = 8)
    par(mar = c(4.4, 4.5, 1, 0.5), cex = 0.85, oma = c(0, 0, 
        1.5, 0))
    layout(matrix(c(1, 2), 1, 2, byrow = TRUE))
    norm.refresh <- function(...) {
        r <- as.numeric(evalq(tclvalue(r), envir = slider.env))
        p <- as.numeric(evalq(tclvalue(p), envir = slider.env))
        xmin <- as.numeric(evalq(tclvalue(xmin), envir = slider.env))
        xmax <- as.numeric(evalq(tclvalue(xmax), envir = slider.env))
        xx <- seq(xmin, xmax)
        yy <- dnbinom(xx, r, p)
        y1 <- pnbinom(xx, r, p)
        plot(xx, yy, type = "h", xlab = expression(italic(x)), 
            ylab = expression(paste(italic(f), "(", italic(x), 
                ")", sep = "")))
        points(xx, yy, pch = 19)
        plot(xx, y1, type = "n", xlab = expression(italic(x)), 
            ylab = expression(paste(italic(F), "(", italic(x), 
                ")", sep = "")))
        points(xx, y1, pch = 19)
        segments(xx, y1, xx + 1, y1)
        points(xx + 1, y1, pch = 1)
        mtext(bquote(paste(italic(X), " ~ ", italic(NB), "(", 
            .(r), ", ", .(p), ")", sep = "")), outer = TRUE, 
            side = 3, cex = 0.9)
        dev.flush()
    }
    tclServiceMode(TRUE)
    m <- tktoplevel()
    tkwm.title(m, "NB(r, \u03C0)")
    tkpack(tklabel(m, text = "      Visualizing the Negative Binomial Distribution      "))
    tkwm.geometry(m, "+0+0")
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "r", font = c("Helvetica", "9", 
        "italic"), width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 1, 
        to = 20, orient = "horiz", resolution = 1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = r), envir = slider.env)
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "\u03C0", font = c("Helvetica", "9", 
        "italic"), width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 0.01, 
        to = 1, orient = "horiz", resolution = 0.01, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = p), envir = slider.env)
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
