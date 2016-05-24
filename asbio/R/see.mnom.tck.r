see.mnom.tck <- function () 
{


    if (!exists("slider.env")) 
        slider.env <- NULL
    suppressWarnings(rm(slider.env))
    slider.env <<- new.env()
    n <- 10
    p1 <- 0.5
    assign("n", tclVar(n), envir = slider.env)
    assign("p1", tclVar(p1), envir = slider.env)

    norm.refresh <- function(...) {
        n <- as.numeric(evalq(tclvalue(n), envir = slider.env))
        p1 <- as.numeric(evalq(tclvalue(p1), envir = slider.env))
        p2 <- 1 - p1
        x1 <- seq(0, n); x2 <- n - x1 
        z <- dbinom(x1, n, p1)
        dev.hold()
        par(mar = c(5,4,0,2))
        scatterplot3d(x = x1, y = x2, z = z, type = "h", zlab = expression(paste(italic(f),"(",italic(y)[1],",",italic(y)[2],")", sep = "")), xlab = expression(italic(y)[1]), ylab = expression(italic(y)[2]), xlim = c(0,20), ylim= c(0,20), pch=19, cex.symbols=.8, box = FALSE)
        mtext(bquote(paste(italic(Y)[1],", ",italic(Y)[2], " ~ ", italic(MNOM), "(", 
            .(n), ", ", .(p1), ", ", .(1-p1),")", sep = "")), line = -3, side = 3)
        dev.flush()
    }
    tclServiceMode(TRUE)
    m <- tktoplevel()
    tkwm.title(m, "MNOM(n, \u03C0\u2081, \u03C0\u2082)")
    tkpack(tklabel(m, text = "      Visualizing the Multinomial Distribution      "))
    tkwm.geometry(m, "+0+0")
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "n", font = c("Helvetica", "9", 
        "italic"), width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 1, 
        to = 20, orient = "horiz", resolution = 1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = n), envir = slider.env)
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "\u03C0\u2081", width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 0, 
        to = 1, orient = "horiz", resolution = 0.01, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = p1), envir = slider.env)
    
    tkpack(tkbutton(m, text = "Refresh", command = norm.refresh), 
        side = "left")
    tkpack(tkbutton(m, text = "Exit", command = function() tkdestroy(m)), 
        side = "right")
}
