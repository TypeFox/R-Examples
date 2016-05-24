see.unif.tck<-function () 
{

    if (!exists("slider.env")) 
        slider.env <- NULL; suppressWarnings(rm(slider.env)); slider.env <<- new.env()# Dummy to trick R CMD check 
    minX <- 2.5
    assign("minX", tclVar(minX), envir = slider.env)
    maxX <- 3
    assign("maxX", tclVar(maxX), envir = slider.env)
    xmin <- 0
    assign("xmin", tclVar(xmin), envir = slider.env)
    xmax <- 5.5
    assign("xmax", tclVar(xmax), envir = slider.env)
    
       
   norm.refresh <- function(...) {
        minX <- as.numeric(evalq(tclvalue(minX), envir = slider.env))
        maxX <- as.numeric(evalq(tclvalue(maxX), envir = slider.env))
        xmin <- as.numeric(evalq(tclvalue(xmin), envir = slider.env))
        xmax <- as.numeric(evalq(tclvalue(xmax), envir = slider.env))
        xx <- seq(xmin, xmax, length = 500)
        yy <- dunif(xx, minX, maxX)
        d <- dunif(minX, minX, maxX)
        umean<-(minX+maxX)/2
        dev.hold()
        plot(xx, yy, type = "n", xlim = c(xmin, xmax), xlab=expression(italic(x)),ylab=expression(paste(italic(f),"(",italic(x),")", sep = "")))
        segments(minX,d,maxX,d) 
        segments(minX,0,minX,d)
        segments(maxX,0,maxX,d)
         mtext(bquote(paste(italic(X), " ~ ", italic(UNIF), "(", .(minX), ", ", .(maxX),")", sep = "")), line = 1, side = 3)
        dev.flush()           
                    }
    tclServiceMode(TRUE)
    m <- tktoplevel()
    tkwm.title(m, "UNIF(a, b)")
    tkpack(tklabel(m,text="      Visualizing the Uniform Distribution      "))
    tkwm.geometry(m, "+0+0")
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "a", font=c("Helvetica","9","italic"),width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 0, 
        to = 2.9, orient = "horiz", resolution = 0.1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = minX), envir = slider.env)
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "b", font=c("Helvetica","9","italic"),width = "20"), 
        side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 3, 
        to = 5.5, orient = "horiz", resolution = 0.1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = maxX), envir = slider.env)
    
    
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



see.unifcdf.tck<-function () 
{

    if (!exists("slider.env")) 
        slider.env <- NULL; suppressWarnings(rm(slider.env)); slider.env <<- new.env()# Dummy to trick R CMD check 
    minX <- 2.5
    assign("minX", tclVar(minX), envir = slider.env)
    maxX <- 3
    assign("maxX", tclVar(maxX), envir = slider.env)
    xmin <- 0
    assign("xmin", tclVar(xmin), envir = slider.env)
    xmax <- 5.5
    assign("xmax", tclVar(xmax), envir = slider.env)
     
    
    dev.new(height=4,width=8);par(mar=c(4.4,4.5,1,0.5),cex=.85, oma = c(0,0,1.5,0)); layout(matrix(c(1,2), 1, 2, byrow = TRUE))
    norm.refresh <- function(...) {
        minX <- as.numeric(evalq(tclvalue(minX), envir = slider.env))
        maxX <- as.numeric(evalq(tclvalue(maxX), envir = slider.env))
        xmin <- as.numeric(evalq(tclvalue(xmin), envir = slider.env))
        xmax <- as.numeric(evalq(tclvalue(xmax), envir = slider.env))
        xx <- seq(xmin, xmax, length = 500)
        yy <- dunif(xx, minX, maxX)
        y1 <- punif(xx, minX, maxX)
        d <- dunif(minX, minX, maxX)
        dev.hold()
        plot(xx, yy, type = "n", xlim = c(xmin, xmax), xlab=expression(italic(x)),ylab=expression(paste(italic(f),"(",italic(x),")", sep = "")))
        segments(minX,d,maxX,d) 
        segments(minX,0,minX,d)
        segments(maxX,0,maxX,d)
        plot(xx, y1, type = "l", xlim = c(xmin, xmax), xlab=expression(italic(x)),ylab=expression(paste(italic(F),"(",italic(x),")", sep = "")))
        mtext(bquote(paste(italic(X), " ~ ", italic(UNIF), "(", .(minX), ", ", .(maxX),")", sep = "")), outer = TRUE, side = 3, cex = .9)
        dev.flush()
            }
    tclServiceMode(TRUE)
    m <- tktoplevel()
     tkwm.title(m, "UNIF(a, b)")
    tkpack(tklabel(m,text="      Visualizing the Uniform Distribution      "))
    tkwm.geometry(m, "+0+0")
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "a", font=c("Helvetica","9","italic"),width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 0, 
        to = 2.9, orient = "horiz", resolution = 0.1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = minX), envir = slider.env)
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "b", font=c("Helvetica","9","italic"),width = "20"), 
        side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 3, 
        to = 5.5, orient = "horiz", resolution = 0.1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = maxX), envir = slider.env)
    
    
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
