see.hyper.tck<-function () 
{

    if (!exists("slider.env")) 
    slider.env <- NULL; suppressWarnings(rm(slider.env)); slider.env <<- new.env()# Dummy to trick R CMD check 
    n <- 1
    M <- 1
    N <- 20
    assign("n", tclVar(n), envir = slider.env)
    assign("M", tclVar(M), envir = slider.env)           
    assign("N", tclVar(N), envir = slider.env)
    norm.refresh <- function(...) {
        n <- as.numeric(evalq(tclvalue(n), envir = slider.env))
        M <- as.numeric(evalq(tclvalue(M), envir= slider.env))
        N <- as.numeric(evalq(tclvalue(N), envir= slider.env))
        xx <- seq(0, n, length = n+1)
        yy <- dhyper(xx,M,N-M,n)
        dev.hold()
        plot(xx, yy, type = "h", xlab=expression(italic(x)),ylab=expression(paste(italic(f),"(",italic(x),")", sep = "")))
        points(xx, yy, pch =19)
        mtext(bquote(paste(italic(X), " ~ ", italic(HYP), "(", .(n), ", ", .(M),", ", .(N),")", sep = "")), line = 1, side = 3)
        dev.flush()
                    }
    tclServiceMode(TRUE)
    m <- tktoplevel()
   tkwm.title(m, "HYP(n, M, n)")
    tkpack(tklabel(m,text="      Visualizing the Hypergeometric Distribution      "))
    tkwm.geometry(m, "+0+0")
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "n", font=c("Helvetica","9","italic"),width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 1, 
        to = 20, orient = "horiz", resolution = 1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir= slider.env)
    evalq(tkconfigure(sc, variable = n), envir= slider.env)
    
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "M", font=c("Helvetica","9","italic"),width = "20"), 
        side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 1, 
        to = 20, orient = "horiz", resolution = 1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir= slider.env)
    evalq(tkconfigure(sc, variable = M), envir= slider.env)  

tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "N", font=c("Helvetica","9","italic"),width = "20"), 
        side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 20, 
        to = 40, orient = "horiz", resolution = 1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir= slider.env)
    evalq(tkconfigure(sc, variable = N), envir= slider.env)  
}

 
see.hypercdf.tck<-function () 
{

    if (!exists("slider.env")) 
    slider.env <- NULL; suppressWarnings(rm(slider.env)); slider.env <<- new.env()# Dummy to trick R CMD check 
    n <- 1
    M <- 1
    N <- 20
    assign("n", tclVar(n), envir= slider.env)
    assign("M", tclVar(M), envir= slider.env)           
    assign("N", tclVar(N), envir= slider.env)
    dev.new(height=4,width=8);par(mar=c(4.4,4.5,1,0.5),cex=.85, oma = c(0,0,1.5,0)); layout(matrix(c(1,2), 1, 2, byrow = TRUE))
   norm.refresh <- function(...) {
        n <- as.numeric(evalq(tclvalue(n), envir= slider.env))
        M <- as.numeric(evalq(tclvalue(M), envir= slider.env))
        N <- as.numeric(evalq(tclvalue(N), envir= slider.env))
        xx <- seq(0, n, length = n+1)
        yy <- dhyper(xx,M,N-M,n)
        y1 <- phyper(xx,M,N-M,n)
        dev.hold()
        plot(xx, yy, type = "h", xlab=expression(italic(x)),ylab=expression(paste(italic(f),"(",italic(x),")", sep = "")))
        points(xx, yy, pch =19)
        plot(xx, y1, type = "n", xlab=expression(italic(x)),ylab=expression(paste(italic(F),"(",italic(x),")", sep = "")))
        points(xx, y1, pch =19)
        segments(xx, y1,xx+1,y1)
        points(xx+1, y1, pch =1)
        mtext(bquote(paste(italic(X), " ~ ", italic(HYP), "(", .(n), ", ", .(M),", ", .(N),")", sep = "")), outer = TRUE, side = 3, cex = .9)
        dev.flush()
                    }
    tclServiceMode(TRUE)
    m <- tktoplevel()
    tkwm.title(m, "HYP(n, M, n)")
    tkpack(tklabel(m,text="      Visualizing the Hypergeometric Distribution      "))
    tkwm.geometry(m, "+0+0")
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "n", font=c("Helvetica","9","italic"),width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 1, 
        to = 20, orient = "horiz", resolution = 1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir= slider.env)
    evalq(tkconfigure(sc, variable = n), envir= slider.env)
    
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "M", font=c("Helvetica","9","italic"),width = "20"), 
        side = "right")
     
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 1, 
        to = 20, orient = "horiz", resolution = 1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir= slider.env)
    evalq(tkconfigure(sc, variable = M), envir= slider.env)  

tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "N", font=c("Helvetica","9","italic"),width = "20"), 
        side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 20, 
        to = 40, orient = "horiz", resolution = 1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir= slider.env)
    evalq(tkconfigure(sc, variable = N), envir= slider.env)  
}     