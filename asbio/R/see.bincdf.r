see.bin.tck<-function () 
{

    if (!exists("slider.env")) 
    slider.env <- NULL; suppressWarnings(rm(slider.env)); slider.env <<- new.env()# Dummy to trick R CMD check 
        slider.env <<- new.env()
    n <- 10
    p <- 0.5
    assign("n", tclVar(n), envir = slider.env)
    assign("p", tclVar(p), envir = slider.env)           
   
   show.norm<-tclVar(0) 
      prefunc<-function(xx,yy,vy,muy,n,p,show.norm=FALSE){
        dev.hold()
        plot(xx, yy, type = "h", xlab=expression(italic(x)),ylab=expression(paste(italic(f),"(",italic(x),")", sep = "")))
        points(xx, yy, pch =19)
        x <- NULL; rm(x); # Dummy to trick R CMD check
        if(show.norm==TRUE) curve(dnorm(x,muy,vy),0,n, col =2,add=TRUE)
        mtext(bquote(paste(italic(X), " ~ ", italic(BIN), "(", .(n), ", ", .(p),")", sep = "")), line = 1, side = 3)
        dev.flush()
                  }   
   norm.refresh <- function(...) {
        n <- as.numeric(evalq(tclvalue(n), envir = slider.env))
        p <- as.numeric(evalq(tclvalue(p), envir = slider.env))
        xx <- seq(0, n, length = n+1)
        yy <- dbinom(xx,n,p)
        show.norm <- as.logical(tclObj(show.norm))
        vy<-sqrt(p*n*(1-p))
        muy<-n*p 
        prefunc(xx,yy,vy,muy,n,p=p,show.norm=show.norm)
                    }
    tclServiceMode(TRUE)
    m <- tktoplevel()
    tkwm.title(m, "BIN(n, \u03C0)")
    tkpack(tklabel(m,text="      Visualizing the Binomial Distribution      "))
    tkwm.geometry(m, "+0+0")
    tkpack(fr <- tkframe(m), side = "top")
    
    tkpack(tklabel(fr, text = "n", font=c("Helvetica","9","italic"),width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 1, 
        to = 30, orient = "horiz", resolution = 1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = n), envir = slider.env)
    
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = '\u03C0', font=c("Helvetica","9","italic"),width = "20"), 
        side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 0, 
        to = 1, orient = "horiz", resolution = 0.01, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = p), envir = slider.env)  
    tkpack(fr,tkcheckbutton(m, text="Show normal approx.",variable=show.norm),anchor="w")
     tkpack(tkbutton(m, text = "Refresh", command = norm.refresh), 
        side = "left")
    tkpack(tkbutton(m, text = "Exit", command = function() tkdestroy(m)), 
        side = "right")
}

 
see.bincdf.tck<-function () 
{

    if (!exists("slider.env")) 
        slider.env <- NULL; suppressWarnings(rm(slider.env)); slider.env <<- new.env()# Dummy to trick R CMD check
    n <- 10
    p <- 0.5
    assign("n", tclVar(n), envir = slider.env)
    assign("p", tclVar(p), envir = slider.env)           
    
    show.norm<-tclVar(0)
    
    dev.new(height=4,width=8);par(mar=c(4.4,4.5,1,0.5),cex=.85, oma = c(0,0,1.5,0)); layout(matrix(c(1,2), 1, 2, byrow = TRUE)) 
   
    prefunc<-function(xx,yy,y1,vy,muy,n,p,show.norm=FALSE){
        dev.hold()
        plot(xx, yy, type = "h", xlab=expression(italic(x)),ylab=expression(paste(italic(f),"(",italic(x),")", sep = "")))
        points(xx, yy, pch =19)
        x <- NULL; rm(x); # Dummy to trick R CMD check
        if(show.norm==TRUE) curve(dnorm(x,muy,vy),0,n, col =2, add = TRUE)
        plot(xx, y1, type = "n", xlab=expression(italic(x)),ylab=expression(paste(italic(F),"(",italic(x),")", sep = "")))
        points(xx, y1, pch =19)
        segments(xx, y1,xx+1,y1)
        points(xx+1, y1, pch =1)
        if(show.norm==TRUE) curve(pnorm(x,muy,vy),0,n, col =2, add=TRUE)
        mtext(bquote(paste(italic(X), " ~ ", italic(BIN), "(", .(n), ", ", .(p),")", sep = "")), outer = TRUE, side = 3, cex = .9)
        dev.flush()
                  }   
   
   norm.refresh <- function(...) {
        n <- as.numeric(evalq(tclvalue(n), envir = slider.env))
        p <- as.numeric(evalq(tclvalue(p), envir = slider.env))
        xx <- seq(0, n, length = n+1)
        yy <- dbinom(xx,n,p)
        y1 <- pbinom(xx,n,p)
        show.norm <- as.logical(tclObj(show.norm))
        vy<-sqrt(p*n*(1-p))
        muy<-n*p 
        prefunc(xx, yy, y1, vy, muy, n, p, show.norm = show.norm)
                   }
    tclServiceMode(TRUE)
    m <- tktoplevel()
    tkwm.title(m, "BIN(n, \u03C0)")
    tkpack(tklabel(m,text="      Visualizing the Binomial Distribution      "))
    tkwm.geometry(m, "+0+0")
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "n", font=c("Helvetica","9","italic"),width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 1, 
        to = 30, orient = "horiz", resolution = 1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = n), envir = slider.env)
    
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = '\u03C0', font=c("Helvetica","9","italic"),width = "20"), 
        side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 0, 
        to = 1, orient = "horiz", resolution = 0.01, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = p), envir = slider.env)  
    tkpack(fr,tkcheckbutton(m, text="Show normal approx.",variable=show.norm),anchor="w")
     tkpack(tkbutton(m, text = "Refresh", command = norm.refresh), 
        side = "left")
    tkpack(tkbutton(m, text = "Exit", command = function() tkdestroy(m)), 
        side = "right")
}
 