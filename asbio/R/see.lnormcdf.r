see.lnorm.tck<-function () 
{

    if (!exists("slider.env")) 
        slider.env <- NULL; suppressWarnings(rm(slider.env)); slider.env <<- new.env()# Dummy to trick R CMD check 
    mu <- 0
    sigma<-1
    assign("mu", tclVar(mu), envir= slider.env)
    assign("sigma", tclVar(sigma), envir= slider.env)
    xmin <- 0
    assign("xmin", tclVar(xmin), envir= slider.env)
    xmax <- 15
    assign("xmax", tclVar(xmax), envir= slider.env)
           
   norm.refresh <- function(...) {
        mu <- as.numeric(evalq(tclvalue(mu), envir= slider.env))
        sigma <- as.numeric(evalq(tclvalue(sigma), envir= slider.env))
        xmin <- as.numeric(evalq(tclvalue(xmin), envir= slider.env))
        xmax <- as.numeric(evalq(tclvalue(xmax), envir= slider.env))
        xx <- seq(xmin, xmax, length = 500)
        yy <- dlnorm(xx, mu, sigma)
        dev.hold()
        plot(xx, yy, type = "l", xlim = c(xmin, xmax), xlab=expression(italic(x)),ylab=expression(paste(italic(f),"(",italic(x),")", sep = "")))
        mtext(bquote(paste(italic(X), " ~ ", italic(LOGN), "(", .(mu), ", ", .(sigma),")", sep = "")), line = 1, side = 3)
        dev.flush()          
                    }
    tclServiceMode(TRUE)
    m <- tktoplevel()
    tkwm.title(m, "LNORM(\u03bc, \u03c3\u00B2)")
    tkpack(tklabel(m,text="      Visualizing the Log-normal Distribution      "))
    tkwm.geometry(m, "+0+0")
    
    
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = '\u03bc',font=c("Helvetica","9","italic"), width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 0, 
        to = 3, orient = "horiz", resolution = 0.1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir= slider.env)
    evalq(tkconfigure(sc, variable = mu), envir= slider.env)
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = '\u03c3',font=c("Helvetica","9","italic"), width = "20"), 
        side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 0.5, 
        to = 3, orient = "horiz", resolution = 0.1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir= slider.env)
    evalq(tkconfigure(sc, variable = sigma), envir= slider.env)
    
     
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "Xmin:", width = 6), side = "left")
    tkpack(e <- tkentry(fr, width = 8), side = "left")
    assign("e", e, envir= slider.env)
    evalq(tkconfigure(e, textvariable = xmin), envir= slider.env)
    tkpack(tklabel(fr, text = "Xmax:", width = 6), side = "left")
    tkpack(e <- tkentry(fr, width = 8), side = "left")
    assign("e", e, envir= slider.env)
    evalq(tkconfigure(e, textvariable = xmax), envir= slider.env)
    tkpack(tkbutton(m, text = "Refresh", command = norm.refresh), 
        side = "left")
    tkpack(tkbutton(m, text = "Exit", command = function() tkdestroy(m)), 
        side = "right")
}


see.lnormcdf.tck<-function (){ 

    if (!exists("slider.env")) 
       slider.env <- NULL; suppressWarnings(rm(slider.env)); slider.env <<- new.env()# Dummy to trick R CMD check 
    mu <- 0
    sigma<-1
    assign("mu", tclVar(mu), envir= slider.env)
    assign("sigma", tclVar(sigma), envir= slider.env)
    xmin <- 0
    assign("xmin", tclVar(xmin), envir= slider.env)
    xmax <- 15
    assign("xmax", tclVar(xmax), envir= slider.env)
           
   dev.new(height=4,width=8);par(mar=c(4.4,4.5,1,0.5),cex=.85, oma = c(0,0,1.5,0)); layout(matrix(c(1,2), 1, 2, byrow = TRUE))
   norm.refresh <- function(...) {
        mu <- as.numeric(evalq(tclvalue(mu), envir= slider.env))
        sigma <- as.numeric(evalq(tclvalue(sigma), envir= slider.env))
        xmin <- as.numeric(evalq(tclvalue(xmin), envir= slider.env))
        xmax <- as.numeric(evalq(tclvalue(xmax), envir= slider.env))
        xx <- seq(xmin, xmax, length = 500)
        yy <- dlnorm(xx, mu, sigma)
        y1<- plnorm(xx, mu, sigma)
        dev.hold()
        plot(xx, yy, type = "l", xlim = c(xmin, xmax), xlab=expression(italic(x)),ylab=expression(paste(italic(f),"(",italic(x),")", sep = "")))
        plot(xx, y1, type = "l", xlim = c(xmin, xmax), xlab=expression(italic(x)),ylab=expression(paste(italic(F),"(",italic(x),")", sep = "")))
        mtext(bquote(paste(italic(X), " ~ ", italic(LOGN), "(", .(mu), ", ", .(sigma),")", sep = "")), outer = TRUE, side = 3, cex = .9)
        dev.flush()          
                    }
    tclServiceMode(TRUE)
    m <- tktoplevel()
    tkwm.title(m, "LNORM(\u03bc, \u03c3\u00B2)")
    tkpack(tklabel(m,text="      Visualizing the Log-normal Distribution      "))
    tkwm.geometry(m, "+0+0")
    
     tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = '\u03bc',font=c("Helvetica","9","italic"), width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 0, 
        to = 3, orient = "horiz", resolution = 0.1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir= slider.env)
    evalq(tkconfigure(sc, variable = mu), envir= slider.env)
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = '\u03c3',font=c("Helvetica","9","italic"), width = "20"), 
        side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 0.5, 
        to = 3, orient = "horiz", resolution = 0.1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir= slider.env)
    evalq(tkconfigure(sc, variable = sigma), envir= slider.env)
    
   
    
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "Xmin:", width = 6), side = "left")
    tkpack(e <- tkentry(fr, width = 8), side = "left")
    assign("e", e, envir= slider.env)
    evalq(tkconfigure(e, textvariable = xmin), envir= slider.env)
    tkpack(tklabel(fr, text = "Xmax:", width = 6), side = "left")
    tkpack(e <- tkentry(fr, width = 8), side = "left")
    assign("e", e, envir= slider.env)
    evalq(tkconfigure(e, textvariable = xmax), envir= slider.env)
    tkpack(tkbutton(m, text = "Refresh", command = norm.refresh), 
        side = "left")
    tkpack(tkbutton(m, text = "Exit", command = function() tkdestroy(m)), 
        side = "right")
  }  
                            