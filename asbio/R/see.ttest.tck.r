see.ttest.tck <- function () 
{

        tclRequire("BWidget")
    if (!exists("slider.env")) 
        slider.env <- NULL; suppressWarnings(rm(slider.env)); slider.env <<- new.env()# Dummy to trick R CMD check
    mu1 <- -1
    assign("mu1", tclVar(mu1), envir = slider.env)
    mu2 <- 0
    assign("mu2", tclVar(mu2), envir = slider.env)
    sigma <- 1
    assign("sigma", tclVar(sigma), envir = slider.env)
    n <- 10
    assign("n", tclVar(n), envir = slider.env)
    
alt <-tclVar("two.sided")
var.equal<-tclVar(0)
    
    ylim <- c(0, dnorm(0, 0, 0.5))
    norm.refresh <- function(...) {
        dev.hold()
        par(mfrow=c(2, 1))
        mu1 <- as.numeric(evalq(tclvalue(mu1), envir = slider.env))
        mu2 <- as.numeric(evalq(tclvalue(mu2), envir = slider.env))
        sigma <- as.numeric(evalq(tclvalue(sigma), envir = slider.env))
        n <- as.numeric(evalq(tclvalue(n), envir = slider.env))
        alt <- tclvalue(alt)       
        ve <- as.logical(tclObj(var.equal))    
        y1 <- rnorm(n,mu1,sigma)
        y2 <- rnorm(n,mu2,sigma)
        
        if(alt == "two.sided")test <- "two"
        if(alt == "less")test <- "lower" 
        if(alt == "greater")test <- "upper"
        x <- rep(c("1","2","3"),each=n)
        
              
        a <- t.test(y1, y2, alternative = alt, var.equal = ve)
        t <- a$statistic
        P <- a$p.value 
        nu <- as.numeric(a$parameter)
               
        cols <- c(rgb(red=0.5,blue=0.6,green=0.5,alpha=.3), rgb(red=0.5,blue=0.5,green=0.6,alpha=.3))
        rcols <- c(rep(rgb(red=0.5,blue=0.6,green=0.5),n), rep(rgb(red=0.5,blue=0.5,green=0.6),n))
        
        curve(dnorm(x,mu1,sigma),from=-5,to=5,ylim=ylim,xlab=expression(italic(x)),ylab=expression(paste(italic(f),"(",italic(x),")", sep = "")))
        
        
        mtext(side = 3, expression(paste(italic(t), "-test from random sample of populations:", sep = "")), line = -1.5) 
        mtext(side = 3, line = -2.5, bquote(paste("alternative = ", .(alt), "  ", italic(t)," = ", .(round(t, 2)), "  ",italic(P)," = ",.(round(P, 5))))) 
        
        
        xx <- seq(-10, 10, length = 500)
        yy1 <- dnorm(xx, mu1, sigma)
        yy2 <- dnorm(xx, mu2, sigma)
        
        polygon(c(xx[xx <= 10], 10), c(yy1[xx <= 10], yy1[xx == 
            -10]), col = cols[1])
        polygon(c(xx[xx <= 10], 10), c(yy2[xx <= 10], yy2[xx == 
            -10]), col = cols[2])
        text(c(y1,y2),rep(0.008,n*2), c(rep(1,n),rep(2,n)), adj = c(0.5,0), cex = .9, col = rcols)
        legend("topright", pch = 22, pt.cex = 1.7, pt.bg = cols, 
            legend = c(expression(italic(X)[1]), expression(italic(X)[2])), 
            bty = "n")
        par(mar=c(4,4,0.5,2))
        shade.t(t, nu=nu, tail= test)
        dev.flush()
    }
    tclServiceMode(TRUE)
    m <- tktoplevel()
    tkwm.title(m, "Visualizing t-tests")
    tkpack(tklabel(m, text = "      Visualizing t-tests      "))
    tkwm.geometry(m, "+0+0")
    
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "\u03BC\u2081", font = c("Helvetica", "10", 
        "normal"), width = "10"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = -4, 
        to = 4, orient = "horiz", resolution = 0.1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = mu1), envir = slider.env)
    
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "\u03BC\u2082", font = c("Helvetica", "10", 
        "normal"), width = "10"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = -4, 
        to = 4, orient = "horiz", resolution = 0.1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = mu2), envir = slider.env)
    
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
    
    alternative <- c("two.sided", "less", "greater")
    comboBox1 <- tkwidget(m,"ComboBox", editable=FALSE, values=alternative, textvariable = alt, width = 10)
      tkpack(tklabel(m, text = "Alternative")) 
      tkpack(comboBox1) 
    v.but <- tkcheckbutton(m, text="Variance equal", variable=var.equal) 
      tkpack(v.but) 
    tkpack(tkbutton(m, text = "Refresh", command = norm.refresh), 
        side = "left")
    tkpack(tkbutton(m, text = "Exit", command = function() tkdestroy(m)), 
        side = "right")
}
