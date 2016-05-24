see.tcdf.tck<-function () 
{

    if (!exists("slider.env")) 
        slider.env <- NULL; suppressWarnings(rm(slider.env)); slider.env <<- new.env()# Dummy to trick R CMD check 
    nu <- 5
    assign("nu", tclVar(nu), envir = slider.env)
    xmin <- -7
    assign("xmin", tclVar(xmin), envir = slider.env)
    xmax <- 7
    assign("xmax", tclVar(xmax), envir = slider.env)
    
    show.norm<-tclVar(0)    
    
    dev.new(height=4,width=8);par(mar=c(4.4,4.5,1,0.5),cex=.85, oma = c(0,0,1.5,0)); layout(matrix(c(1,2), 1, 2, byrow = TRUE))
    
    prefunc<-function(xx,yy,sny,y1,sncy,xlim,nu,show.norm=FALSE){
        dev.hold()
        plot(xx, yy, type = "l", xlim = c(xmin, xmax), xlab=expression(italic(x)),ylab=expression(paste(italic(f),"(",italic(x),")", sep = "")))
        if(show.norm==TRUE) points(xx,sny,type = "l", col =2)
        plot(xx, y1, type = "l", xlim = c(xmin, xmax), xlab=expression(italic(x)),ylab=expression(paste(italic(F),"(",italic(x),")", sep = "")))
        if(show.norm==TRUE) points(xx,sncy,type = "l", col =2)
        mtext(bquote(paste(italic(X), " ~ ", italic(t), "(", .(nu), ")", sep = "")), outer = TRUE, side = 3, cex = .9)
        dev.flush()
        }
    
      
    norm.refresh <- function(...) {
        nu <- as.numeric(evalq(tclvalue(nu), envir = slider.env))
        xmin <- as.numeric(evalq(tclvalue(xmin), envir = slider.env))
        xmax <- as.numeric(evalq(tclvalue(xmax), envir = slider.env))
        xx <- seq(xmin, xmax, length = 500)
        yy <- dt(xx, df=nu)
        y1 <- pt(xx, df=nu)
        sny <- dnorm(xx)
        sncy <- pnorm(xx)
        show.norm <- as.logical(tclObj(show.norm))
        prefunc(xx,yy,sny,y1,sncy,nu=nu,show.norm=show.norm)
              }
              
    tclServiceMode(TRUE)
    m <- tktoplevel()
    tkwm.title(m, "t(\u03bd)")
    tkpack(tklabel(m,text="      Visualizing the t Distribution      "))
    tkwm.geometry(m, "+0+0")
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = '\u03bd',font=c("Helvetica","9","italic"), width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 1, 
        to = 15, orient = "horiz", resolution = 1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = nu), envir = slider.env)
    
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr,text="      "))  
    tkpack(fr,tkcheckbutton(m, text="Show standard normal pdf",variable=show.norm),anchor="w")
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr,text="      "))  
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


see.t.tck<-function () 
{

    if (!exists("slider.env")) 
        slider.env <- NULL; suppressWarnings(rm(slider.env)); slider.env <<- new.env()# Dummy to trick R CMD check 
    nu <- 5
    assign("nu", tclVar(nu), envir = slider.env)
    xmin <- -7
    assign("xmin", tclVar(xmin), envir = slider.env)
    xmax <- 7
    assign("xmax", tclVar(xmax), envir = slider.env)
        
     show.norm<-tclVar(0)    
    
    prefunc<-function(xx,yy,sny,y1,sncy,xlim,nu,show.norm=FALSE){
        dev.hold()
        plot(xx, yy, type = "l", xlim = c(xmin, xmax), xlab=expression(italic(x)),ylab=expression(paste(italic(f),"(",italic(x),")", sep = "")))
        if(show.norm==TRUE) points(xx,sny,type = "l", col =2)
        mtext(bquote(paste(italic(X), " ~ ", italic(t), "(", .(nu), ")", sep = "")), line = 1, side = 3)
        dev.flush()         
                 }
    
    norm.refresh <- function(...) {
        nu <- as.numeric(evalq(tclvalue(nu), envir = slider.env))
        xmin <- as.numeric(evalq(tclvalue(xmin), envir = slider.env))
        xmax <- as.numeric(evalq(tclvalue(xmax), envir = slider.env))
        xx <- seq(xmin, xmax, length = 500)
        yy <- dt(xx, nu)
        sny <- dnorm(xx)
        show.norm <- as.logical(tclObj(show.norm))
        prefunc(xx,yy,sny,nu=nu,show.norm=show.norm)
              }
    
    tclServiceMode(TRUE)
    m <- tktoplevel()
    tkwm.title(m, "t(\u03bd)")
    tkpack(tklabel(m,text="      Visualizing the t Distribution      "))
    tkwm.geometry(m, "+0+0")
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = '\u03b7',font=c("Helvetica","9","italic"), width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 1, 
        to = 15, orient = "horiz", resolution = 1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = nu), envir = slider.env)
        
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr,text="      "))  
    tkpack(fr,tkcheckbutton(m, text="Show standard normal pdf",variable=show.norm),anchor="w")
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr,text="      "))  
    
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