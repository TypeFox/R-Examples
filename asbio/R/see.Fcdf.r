see.Fcdf.tck<-function () 
{

    if (!exists("slider.env")) 
    slider.env <- NULL; suppressWarnings(rm(slider.env)); slider.env <<- new.env()# Dummy to trick R CMD check 
    nu1 <- 5
    assign("nu1", tclVar(nu1), envir = slider.env)
    nu2 <- 5
    assign("nu2", tclVar(nu2), envir = slider.env)
    xmin <- 0
    assign("xmin", tclVar(xmin), envir = slider.env)
    xmax <- 10
    assign("xmax", tclVar(xmax), envir = slider.env)
        
    dev.new(height=4,width=8);par(mar=c(4.4,4.5,1,0.5),cex=.85, oma = c(0,0,1.5,0)); layout(matrix(c(1,2), 1, 2, byrow = TRUE))
    norm.refresh <- function(...) {
        nu1 <- as.numeric(evalq(tclvalue(nu1), envir = slider.env))
        nu2 <- as.numeric(evalq(tclvalue(nu2), envir = slider.env))
        xmin <- as.numeric(evalq(tclvalue(xmin), envir = slider.env))
        xmax <- as.numeric(evalq(tclvalue(xmax), envir = slider.env))
        xx <- seq(0, 10, length = 200)
        yy <- df(xx, nu1, nu2)
        y1 <- pf(xx, nu1, nu2)
        dev.hold()
        plot(xx, yy, type = "l", xlim = c(xmin, xmax),xlab=expression(italic(x)),ylab=expression(paste(italic(f),"(",italic(x),")", sep = "")))
        plot(xx, y1, type = "l", xlim = c(xmin, xmax),xlab=expression(italic(x)),ylab=expression(paste(italic(F),"(",italic(x),")", sep = "")))
        mtext(bquote(paste(italic(X), " ~ ", italic(F), "(", .(nu1), ", ", .(nu2),")", sep = "")), outer = TRUE, side = 3, cex = .9)
        dev.flush()
            }
    tclServiceMode(TRUE)
    m <- tktoplevel()
    tkwm.title(m, "F(\u03bd\u2081, \u03bd\u2082)")
    tkpack(tklabel(m,text="      Visualizing the F Distribution      "))
    tkwm.geometry(m, "+0+0")
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = '\u03bd\u2081',width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 1, 
        to = 30, orient = "horiz", resolution = 1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = nu1), envir = slider.env)
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = '\u03bd\u2082',width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 1, 
        to = 30, orient = "horiz", resolution = 1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    tkpack(fr <- tkframe(m), side = "top")
    evalq(tkconfigure(sc, variable = nu2), envir = slider.env)
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


see.F.tck<-function () 
{

    if (!exists("slider.env"))
    slider.env <- NULL; suppressWarnings(rm(slider.env)); slider.env <<- new.env()# Dummy to trick R CMD check 
    nu1 <- 5
    assign("nu1", tclVar(nu1), envir = slider.env)
    nu2 <- 5
    assign("nu2", tclVar(nu2), envir = slider.env)
        xmin <- 0
    assign("xmin", tclVar(xmin), envir = slider.env)
    xmax <- 10
    assign("xmax", tclVar(xmax), envir = slider.env)
        
    norm.refresh <- function(...) {
        nu1 <- as.numeric(evalq(tclvalue(nu1), envir = slider.env))
        nu2 <- as.numeric(evalq(tclvalue(nu2), envir = slider.env))
        xmin <- as.numeric(evalq(tclvalue(xmin), envir = slider.env))
        xmax <- as.numeric(evalq(tclvalue(xmax), envir = slider.env))
        xx <- seq(0, 10, length = 200)
        yy <- df(xx, nu1, nu2)
        dev.hold()
        plot(xx, yy, type = "l", xlim = c(xmin, xmax),xlab=expression(italic(x)),ylab=expression(paste(italic(f),"(",italic(x),")", sep = "")))
        mtext(bquote(paste(italic(X), " ~ ", italic(F), "(", .(nu1), ", ", .(nu2),")", sep = "")), line = 1, side = 3)
        dev.flush()
                    }
    tclServiceMode(TRUE)
    m <- tktoplevel()
    tkwm.title(m, "F(\u03bd\u2081, \u03bd\u2082)")
    tkpack(tklabel(m,text="      Visualizing the F Distribution      "))
    tkwm.geometry(m, "+0+0")
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = '\u03bd\u2081', width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 1, 
        to = 30, orient = "horiz", resolution = 1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = nu1), envir = slider.env)
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = '\u03bd\u2082', width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 1, 
        to = 30, orient = "horiz", resolution = 1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = nu2), envir = slider.env)
        tkpack(fr <- tkframe(m), side = "top")
    evalq(tkconfigure(sc, variable = nu2), envir = slider.env)
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