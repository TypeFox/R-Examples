see.exppower.tck<-function () 
{

    if (!exists("slider.env")) 
        slider.env <- NULL; suppressWarnings(rm(slider.env)); slider.env <<- new.env()# Dummy to trick R CMD check 
    m <- 1
    assign("m", tclVar(m), envir = slider.env)
    xmin <- -3
    assign("xmin", tclVar(xmin), envir = slider.env)
    xmax <- 3
    assign("xmax", tclVar(xmax), envir = slider.env)
        
    exppower<-function(x,m){exp(-abs(x)^m)}
    
    norm.refresh <- function(...) {
        m <- as.numeric(evalq(tclvalue(m), envir = slider.env))
        xmin <- as.numeric(evalq(tclvalue(xmin), envir = slider.env))
        xmax <- as.numeric(evalq(tclvalue(xmax), envir = slider.env))
        xx <- seq(xmin, xmax, length = 200)
        yy <- exppower(as.numeric(xx), m)
        dev.hold()
        plot(xx, yy, type = "l", xlim = c(xmin, xmax),xlab=expression(italic(x)),ylab=expression(paste(italic(f),"(",italic(x),")", sep = "")))
        dev.flush()        
                }
    tclServiceMode(TRUE)
    m <- tktoplevel()
    tkwm.title(m, "Visualizing Exponential Power Distributions")
    tkwm.geometry(m, "+0+0")
    tkpack(tklabel(m,text="      Visualizing Exponential Power Distributions      "))
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "f(x) =", font=c("Helvetica","9","italic")), side = "left")
    tkpack(tklabel(fr, text = "exp(", font=c("Helvetica","9","normal")),side="left")
    tkpack(tklabel(fr, text ='-|x| \u036b', font=c("Helvetica","9","italic")), side = "left", anchor ="w")
    tkpack(tklabel(fr, text = " )  ", font=c("Helvetica","9","normal")),side="left", anchor ="w")
    tkpack(tklabel(m,text=""), side = "top")
    tkpack(fr1 <- tkframe(m), side = "top")
    tkpack(tklabel(fr1, text = "m  ", font=c("Helvetica","10","italic")),side="left", anchor = "s")
    tkpack(sc <- tkscale(fr1, command = norm.refresh, from = -10, 
        to = 10, orient = "horiz", resolution = .1, showvalue = TRUE), 
        side = "left", anchor="n")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = m), envir = slider.env)
    
    tkpack(tklabel(m,text=""), side = "top")
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "Xmin:", width = 6), side = "left")
    tkpack(e <- tkentry(fr, width = 8), side = "left")
    assign("e", e, envir = slider.env)
    evalq(tkconfigure(e, textvariable = xmin), envir = slider.env)
    tkpack(tklabel(fr, text = "Xmax:", width = 6), side = "left")
    tkpack(e <- tkentry(fr, width = 8), side = "left")
    assign("e", e, envir = slider.env)
    evalq(tkconfigure(e, textvariable = xmax), envir = slider.env)
    tkpack(tklabel(m,text=""), side = "top")
    tkpack(tkbutton(m, text = "Refresh", command = norm.refresh), 
        side = "left")
    tkpack(tkbutton(m, text = "Exit", command = function() tkdestroy(m)), 
        side = "right")
 }


