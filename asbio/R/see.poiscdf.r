see.pois.tck<-function () 
{

    if (!exists("slider.env")) 
        slider.env <- NULL; suppressWarnings(rm(slider.env)); slider.env <<- new.env()# Dummy to trick R CMD check 
    lambda <- 5
    assign("lambda", tclVar(lambda), envir  = slider.env)
    xmin <- 0
    assign("xmin", tclVar(xmin), envir  = slider.env)
    xmax <- 20
    assign("xmax", tclVar(xmax), envir  = slider.env)
       
    show.norm<-tclVar(0) 
      prefunc<-function(xx,yy,vy,muy,xmin,xmax,lambda,show.norm=FALSE){
        dev.hold()
        plot(xx, yy, type = "h", xlim = c(xmin, xmax), xlab=expression(italic(x)),ylab=expression(paste(italic(f),"(",italic(x),")", sep = "")))
        points(xx, yy, pch =19)
        x <- NULL; rm(x); # Dummy to trick R CMD check 
        if(show.norm==TRUE) curve(dnorm(x,muy,vy),0,xmax, col =2,add=TRUE)
        mtext(bquote(paste(italic(X), " ~ ", italic(POI), "(", .(lambda), ")", sep = "")), line = 1, side = 3)
        dev.flush()
                  } 
        
        norm.refresh <- function(...) {
        lambda <- as.numeric(evalq(tclvalue(lambda), envir  = slider.env))
        xmin <- as.numeric(evalq(tclvalue(xmin), envir  = slider.env))
        xmax <- as.numeric(evalq(tclvalue(xmax), envir  = slider.env))
        xx <- seq(xmin, xmax, length = (xmax-xmin)+1)
        yy <- dpois(xx, lambda)
        show.norm <- as.logical(tclObj(show.norm))
        vy<-sqrt(lambda)
        muy<-lambda 
        prefunc(xx,yy,vy,muy,xmin,xmax,lambda=lambda,show.norm=show.norm)
                    }
    tclServiceMode(TRUE)
    m <- tktoplevel()
    tkwm.title(m, "POI(\u03bb)")
    tkpack(tklabel(m,text="      Visualizing the Poisson Distribution      "))
    tkwm.geometry(m, "+0+0")
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = '\u03bb', font=c("Helvetica","9","italic"), width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 0.1, 
        to = 10, orient = "horiz", resolution = 0.1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir  = slider.env)
    evalq(tkconfigure(sc, variable = lambda), envir  = slider.env)
   
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "Xmin:", width = 6), side = "left")
    tkpack(e <- tkentry(fr, width = 8), side = "left")
    assign("e", e, envir  = slider.env)
    evalq(tkconfigure(e, textvariable = xmin), envir  = slider.env)
    tkpack(tklabel(fr, text = "Xmax:", width = 6), side = "left")
    tkpack(e <- tkentry(fr, width = 8), side = "left")
    assign("e", e, envir  = slider.env)
    evalq(tkconfigure(e, textvariable = xmax), envir  = slider.env)
    tkpack(fr,tkcheckbutton(m, text="Show normal approx.",variable=show.norm),anchor="w")
    tkpack(tkbutton(m, text = "Refresh", command = norm.refresh), 
        side = "left")
    tkpack(tkbutton(m, text = "Exit", command = function() tkdestroy(m)), 
        side = "right")
}


see.poiscdf.tck<-function () 
{

    if (!exists("slider.env")) 
        slider.env <- NULL; suppressWarnings(rm(slider.env)); slider.env <<- new.env()# Dummy to trick R CMD check 
    lambda <- 5
    assign("lambda", tclVar(lambda), envir  = slider.env)
    xmin <- 0
    assign("xmin", tclVar(xmin), envir  = slider.env)
    xmax <- 20
    assign("xmax", tclVar(xmax), envir  = slider.env)
   
   show.norm<-tclVar(0)  
   dev.new(height=4,width=8);par(mar=c(4.4,4.5,1,0.5),cex=.85, oma = c(0,0,1.5,0)); layout(matrix(c(1,2), 1, 2, byrow = TRUE))
   prefunc<-function(xx,yy,vy,muy,y1,xmin, xmax,lambda, show.norm=FALSE){
        dev.hold()
        plot(xx, yy, type = "h", xlim = c(xmin, xmax),xlab=expression(italic(x)),ylab=expression(paste(italic(f),"(",italic(x),")", sep = "")))
        points(xx, yy, pch =19)
        x <- NULL; rm(x); # Dummy to trick R CMD check 
        if(show.norm==TRUE) curve(dnorm(x,muy,vy),0,xmax, col =2,add=TRUE)
        plot(xx, y1, type = "n", xlim = c(xmin, xmax), xlab=expression(italic(x)),ylab=expression(paste(italic(F),"(",italic(x),")", sep = "")))
        segments(xx, y1, xx+1, y1, pch=19)
        points(xx,y1,pch=19)
        points(xx+1,y1,pch=1)
        if(show.norm==TRUE) curve(pnorm(x,muy,vy),0,xmax, col =2,add=TRUE)
        mtext(bquote(paste(italic(X), " ~ ", italic(POI), "(", .(lambda), ")", sep = "")), outer = TRUE, side = 3, cex = .9)
        dev.flush()
                  } 
   
   
   norm.refresh <- function(...) {
        lambda <- as.numeric(evalq(tclvalue(lambda), envir  = slider.env))
        xmin <- as.numeric(evalq(tclvalue(xmin), envir  = slider.env))
        xmax <- as.numeric(evalq(tclvalue(xmax), envir  = slider.env))
        xx <- seq(xmin, xmax, length = (xmax-xmin)+1)
        yy <- dpois(xx, lambda)
        y1<- ppois(xx, lambda)
        show.norm <- as.logical(tclObj(show.norm))
        vy<-sqrt(lambda)
        muy<-lambda 
        prefunc(xx,yy,vy,muy,y1,xmin,xmax,lambda,show.norm=show.norm)
                            }
    tclServiceMode(TRUE)
    m <- tktoplevel()
    tkwm.title(m, "POI(\u03bb)")
    tkpack(tklabel(m,text="      Visualizing the Poisson Distribution      "))
    tkwm.geometry(m, "+0+0")
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = '\u03bb', font=c("Helvetica","9","italic"), width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 0.1, 
        to = 10, orient = "horiz", resolution = 0.1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir  = slider.env)
    evalq(tkconfigure(sc, variable = lambda), envir  = slider.env)
   
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "Xmin:", width = 6), side = "left")
    tkpack(e <- tkentry(fr, width = 8), side = "left")
    assign("e", e, envir  = slider.env)
    evalq(tkconfigure(e, textvariable = xmin), envir  = slider.env)
    tkpack(tklabel(fr, text = "Xmax:", width = 6), side = "left")
    tkpack(e <- tkentry(fr, width = 8), side = "left")
    assign("e", e, envir  = slider.env)
    evalq(tkconfigure(e, textvariable = xmax), envir  = slider.env)
    tkpack(fr,tkcheckbutton(m, text="Show normal approx.",variable=show.norm),anchor="w")
    tkpack(tkbutton(m, text = "Refresh", command = norm.refresh), 
        side = "left")
    tkpack(tkbutton(m, text = "Exit", command = function() tkdestroy(m)), 
        side = "right")
}


