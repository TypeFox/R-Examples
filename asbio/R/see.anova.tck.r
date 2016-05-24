see.anova.tck <- function () 
{

    if (!exists("slider.env")) 
        slider.env <- NULL; suppressWarnings(rm(slider.env)); slider.env <<- new.env()# Dummy to trick R CMD check
    mu1 <- -1.5
    assign("mu1", tclVar(mu1), envir = slider.env)
    mu2 <- 0
    assign("mu2", tclVar(mu2), envir = slider.env)
    mu3 <- 1
    assign("mu3", tclVar(mu3), envir = slider.env)
    sigma <- 1
    assign("sigma", tclVar(sigma), envir = slider.env)
    n <- 10
    assign("n", tclVar(n), envir = slider.env)
    color<-tclVar(1)
    
    ylim <- c(0, dnorm(0, 0, 0.5))
    norm.refresh <- function(...) {
        dev.hold()
        mu1 <- as.numeric(evalq(tclvalue(mu1), envir = slider.env))
        mu2 <- as.numeric(evalq(tclvalue(mu2), envir = slider.env))
        mu3 <- as.numeric(evalq(tclvalue(mu3), envir = slider.env))
        sigma <- as.numeric(evalq(tclvalue(sigma), envir = slider.env))
        n <- as.numeric(evalq(tclvalue(n), envir = slider.env))
        color <- as.logical(tclObj(color))
        
        mu <- (mu1+mu2+mu3)/3
        E.MSA <- sigma^2 + (n*((mu-mu1)^2)/2 + n*((mu-mu2)^2)/2 + n*((mu-mu3)^2)/2)
        
        
        y1 <- rnorm(n,mu1,sigma)
        y2 <- rnorm(n,mu2,sigma)
        y3 <- rnorm(n,mu3,sigma)
        x <- rep(c("1","2","3"),each=n)
        
        y <- c(y1,y2,y3)
        a <- anova(lm(y~factor(x)))
        MSA <- a$"Mean Sq"[1]
        MSE <- a$"Mean Sq"[2]
        F.star <- MSA/MSE
        P <- pf(F.star,2,n-3,lower.tail=FALSE) 
               
        cols <- c(rgb(red=0.4,blue=0.4,green=0.5,alpha=.3), rgb(red=0.4,blue=0.8,green=0.5,alpha=.2),rgb(red=0.8,blue=0.4,green=0.5,alpha=.3))
        rcols <- c(rep(rgb(red=0.4,blue=0.4,green=0.5),n), rep(rgb(red=0.4,blue=0.8,green=0.5),n),rep(rgb(red=0.8,blue=0.4,green=0.5),n))
        
        curve(dnorm(x,mu1,sigma),from=-5,to=5,ylim=ylim,xlab=expression(italic(y)),ylab=expression(paste(italic(f),"(",italic(y),")", sep = "")))
        mtext(side = 3, "ANOVA from random sample of populations:", line = 2.5) 
        mtext(side = 3, line = 1, bquote(paste(italic(MSTR)," = ", .(round(MSA, 2)), "  ",italic(MSE)," = ", 
        .(round(MSE, 2)), "  ",italic(F),  "* = ", .(round(F.star, 2)),"   ", italic(P), "-value = ", .(round(P, 5)))))
        
        xx <- seq(-10, 10, length = 500)
        yy1 <- dnorm(xx, mu1, sigma)
        yy2 <- dnorm(xx, mu2, sigma)
        yy3 <- dnorm(xx, mu3, sigma)
    
    if(color == TRUE){
        polygon(c(xx[xx <= 10], 10), c(yy1[xx <= 10], yy1[xx == 
            -10]), col = cols[1])
        polygon(c(xx[xx <= 10], 10), c(yy2[xx <= 10], yy2[xx == 
            -10]), col = cols[2])
        polygon(c(xx[xx <= 10], 10), c(yy3[xx <= 10], yy3[xx == 
            -10]), col = cols[3])
        mtext(side = 1, at = c(y1,y2,y3), line = -1, c(rep(1,n),rep(2,n),rep(3,n)), adj = c(0.5,0), cex = .9, col = rcols)
        legend("topright", pch = 22, pt.cex = 1.7, pt.bg = cols,  
            legend = c(expression(italic(X)[1]), expression(italic(X)[2]), expression(italic(X)[3])), 
            bty = "n")}
        
    if(color == FALSE){
        polygon(c(xx[xx <= 10], 10), c(yy1[xx <= 10], yy1[xx == 
            -10]), angle = 45, density = 20)
        polygon(c(xx[xx <= 10], 10), c(yy2[xx <= 10], yy2[xx == 
            -10]), angle = 0, density = 20)
        polygon(c(xx[xx <= 10], 10), c(yy3[xx <= 10], yy3[xx == 
            -10]), angle = 135, density = 20)
        mtext(side = 1, at = c(y1,y2,y3), line = -1, c(rep(1,n),rep(2,n),rep(3,n)), adj = c(0.5,0), cex = .9, col = 1)
        legend("topright", angle = c(45, 0, 135), density = 20,  
            legend = c(expression(italic(X)[1]), expression(italic(X)[2]), expression(italic(X)[3])), 
            bty = "n")}
            
        legend("topleft", legend = c(expression(paste(italic(alpha)[1], 
            "  =")), expression(paste(italic(alpha)[2], "  =")), expression(paste(italic(alpha)[3], 
            "  = ")), "", expression(paste(italic(E), "(", italic(MSTR), 
            ") = ")), expression(paste(italic(E), "(", italic(MSE), 
            ") = "))), bty = "n")
        legend(-4.7, 0.833, legend = c(round(mu1 - mu, 2), round(mu2 - 
            mu, 2), round(mu3 - mu, 2),"", paste("              ", 
            round(E.MSA, 2)), paste("            ", round(sigma^2, 
            2))), bty = "n")
        
        abline(v=mu,lty=2)
        text(mu + 0.2, max(ylim),expression(mu))
        dev.flush()
    }
    tclServiceMode(TRUE)
    m <- tktoplevel()
    tkwm.title(m, "Visualizing ANOVA")
    tkpack(tklabel(m, text = "      Visualizing ANOVA      "))
    tkwm.geometry(m, "+0+0")
    
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "\u03BC\u2081", font = c("Helvetica", "11", 
        "normal"), width = "10"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = -4, 
        to = 4, orient = "horiz", resolution = 0.1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = mu1), envir = slider.env)
    
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "\u03BC\u2082", font = c("Helvetica", "11", 
        "normal"), width = "10"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = -4, 
        to = 4, orient = "horiz", resolution = 0.1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = mu2), envir = slider.env)
    
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "\u03BC\u2083", font = c("Helvetica", "11", 
        "normal"), width = "10"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = -4, 
        to = 4, orient = "horiz", resolution = 0.1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = mu3), envir = slider.env)
    
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "\u03C3", font = c("Helvetica", "10", 
        "normal"), width = "10"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 0.1, 
        to = 2, orient = "horiz", resolution = 0.1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = sigma), envir = slider.env)
  
   tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "n", font = c("Helvetica", "10", 
        "normal"), width = "10"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 2, 
        to = 20, orient = "horiz", resolution = 1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = n), envir = slider.env) 
    tkpack(tkcheckbutton(m, text="Color", variable=color))
    
    tkpack(tkbutton(m, text = "Refresh", command = norm.refresh), 
        side = "left")
    tkpack(tkbutton(m, text = "Exit", command = function() tkdestroy(m)), 
        side = "right")
}
