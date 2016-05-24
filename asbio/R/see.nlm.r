see.nlm<-function(){
if(any(dev.list()>1)) graphics.off()
options()$device(xpos=20)
par(mar=c(.1,.1,.1,.1))
plot(seq(1,10),seq(1,10.5,length.out=10),type="n",xaxt="n",yaxt="n",xlab="",ylab="")
text(5.8,10.2,"Important non-linear models",cex=1.5)
rect(2.25, 1.15, 9.75, 9.75,lwd=1.5)
segments(2.25,9.1,9.75,9.1,lwd=1.5)
segments(5,9.75,5,1.15,lwd=1.5)
segments(8,9.75,8,1.15,lwd=1.5)
segments(2.25,8.45,9.75,8.45, col="gray")
segments(2.25,7.45,9.75,7.45, col="gray")
segments(2.25,6.45,9.75,6.45, col="gray")
segments(2.25,5.8,9.75,5.8, col="gray")
segments(2.25,4.8,9.75,4.8, col="gray")
segments(2.25,3.8,9.75,3.8, col="gray")
segments(2.25,2.8,9.75,2.8, col="gray")
segments(2.25,2.15,9.75,2.15, col="gray")
segments(2.25,1.15,9.75,1.15)
text(3.1,8.6,"Asymptotic",font=2)
text(3.05,5.95,"Sigmoidal",font=2)
text(3.25,2.3,"Hump-shaped",font=2)
text(3.5,9.25,"Name",font=2)
text(6.5,9.25,"Equation",font=2)
text(8.9,9.25,"Parameters",font=2)
text(3.31,7.6,"Michaelis-Menten",cex=.9)
text(3.65,6.6,"2-parameter exponential",cex=.9)
text(3.4,4.95,"2-parameter logistic",cex=.9)
text(3.4,3.95,"3-parameter logistic",cex=.9)
text(2.85,2.95,"Gompertz",cex=.9)
text(2.71,1.3,"Ricker",cex=.9)
text(8.8,7.6,"a, b",font=3, cex=.9)
text(8.8,6.6,"a, b",font=3, cex=.9)
text(8.8,4.95,"a, b",font=3,cex=.9)
text(8.85,3.95,"a, b, c",font=3,cex=.9)
text(8.85,2.95,"a, b, c",font=3,cex=.9)
text(8.8,1.3,"a, b",font=3,cex=.9)
text(6.3,7.9,expression(paste(italic("f(x) = "),italic(frac(ax,b+x)))),font=3)
text(6.3,6.9,expression(paste(italic("f(x) = a","(1-"),italic(e^(-bx)))),font=3)
text(6.3,5.25,expression(paste(italic("f(x) = "),italic(frac(e^(a+bx),1+e^(a+bx))))),font=3)
text(6.3,4.25,expression(paste(italic("f(x) = "),italic(frac(a,1+be^(-cx))))),font=3)
text(6.3,3.25,expression(paste(italic("f(x) = a"),italic(e^(-be)^(-cx)))),font=3)
text(6.3,1.6,expression(paste(italic("f(x) = ax"),italic(e^(-bx)))),font=3)
x<-rep(1.5,5,6)
y<-c(7.95,6.95,5.3,4.3,3.3,1.65)
points(x,y,cex=1.5)

fp<-function(){
ans <- identify(x, y, n = 1, plot = FALSE)
yw <- y[ans]
if(yw==7.95)
{points(1.5,7.95,pch=21,bg="red",cex=1.5);dev.new();com="see.MM"}
if(yw==6.95)
{points(1.5,6.95,pch=21,bg="red",cex=1.5);dev.new();com="see.2PE"}
if(yw==5.3)
{points(1.5,5.3,pch=21,bg="red",cex=1.5);dev.new();com="see.2PL"}
if(yw==4.3)
{points(1.5,4.3,pch=21,bg="red",cex=1.5);dev.new();com="see.3PL"}
if(yw==3.3)
{points(1.5,3.3,pch=21,bg="red",cex=1.5);dev.new();com="see.G"}
if(yw==1.65)
{points(1.5,1.65,pch=21,bg="red",cex=1.5);dev.new();com="see.R"}
com
}

com <- fp()

    if (!exists("slider.env")) 
        slider.env <- NULL; suppressWarnings(rm(slider.env)); slider.env <<- new.env()# Dummy to trick R CMD check
    a <- 1
    b <- 1
    c <- 1
    assign("a", tclVar(a), envir= slider.env)
    assign("b", tclVar(b), envir= slider.env)
    assign("c", tclVar(c), envir= slider.env)
    xmin <- 0
    assign("xmin", tclVar(xmin), envir= slider.env)
    xmax <- 10
    assign("xmax", tclVar(xmax), envir= slider.env)
        
   
    
    norm.refresh <- function(...) {
        a <- as.numeric(evalq(tclvalue(a), envir= slider.env))
        b <- as.numeric(evalq(tclvalue(b), envir= slider.env))
        c <- as.numeric(evalq(tclvalue(c), envir= slider.env))
        xmin <- as.numeric(evalq(tclvalue(xmin), envir= slider.env))
        xmax <- as.numeric(evalq(tclvalue(xmax), envir= slider.env))
        xx <- seq(xmin, xmax, length = 500)
        par(mar=c(5, 4.3, 4, 2), cex=1.2)
        
        if(com=="see.MM"){yy<-a*xx/(b+xx);main=bquote(paste("Michaelis-Menten Model   ",italic(f),"(",italic(x),") = ",.(a),italic(x),"/(",.(b)," + ",italic(x),")"))}
        if(com=="see.2PE"){yy<-a*exp(-b*xx);main="2 Parameter Exponential"}
        if(com=="see.2PL"){yy<-exp(a+b*xx)/(1+exp(a+b*xx));main="2 Parameter Logistic"}
        if(com=="see.3PL"){yy<-a/(1+b*exp(-c*xx));main="3 Parameter Logistic"}
        if(com=="see.G"){yy<-a*exp(-b*exp(-c*xx));main="Gompertz"}
        if(com=="see.R"){yy<-a*xx*exp(-b*xx);main="Ricker"}
        dev.hold()
        par(cex=1.3)
        plot(xx, yy, type = "l", xlim = c(xmin, xmax), xlab=expression(italic(x)),ylab=expression(paste(italic(f),"(",italic(x),")", sep = "")),main=main, cex.main=1.1)
        dev.flush()    
    }
    
    tw <- function(){
    tkdestroy(m)
    see.nlm()
    }
    
    
    tclServiceMode(TRUE)
    m <- tktoplevel()
    tkwm.geometry(m, "+600+4")
    tkwm.title(m, "Visualizing Non-linear Models")
    tkpack(tklabel(m, text = "      Visualizing Non-linear Models      "))
    tkwm.geometry(m, "+0+0")
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "a", font = c("Helvetica", 
        "9", "italic"), width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 0, 
        to = 50, orient = "horiz", resolution = 0.1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir= slider.env)
    evalq(tkconfigure(sc, variable = a), envir= slider.env)
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "b", font = c("Helvetica", 
        "9", "italic"), width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 0, 
        to = 20, orient = "horiz", resolution = 0.1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir= slider.env)
    evalq(tkconfigure(sc, variable = b), envir= slider.env)
     tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "c", font = c("Helvetica", 
        "9", "italic"), width = "20"), side = "right")
    tkpack(sc <- tkscale(fr, command = norm.refresh, from = 0, 
        to = 10, orient = "horiz", resolution = 0.1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir= slider.env)
    evalq(tkconfigure(sc, variable = c), envir= slider.env)
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "Xmin:", width = 6), side = "left")
    tkpack(e <- tkentry(fr, width = 8), side = "left")
    assign("e", e, envir= slider.env)
    evalq(tkconfigure(e, textvariable = xmin), envir= slider.env)
    tkpack(tklabel(fr, text = "Xmax:", width = 6), side = "left")
    tkpack(e <- tkentry(fr, width = 8), side = "left")
    assign("e", e, envir= slider.env)
    evalq(tkconfigure(e, textvariable = xmax), envir= slider.env)
    tkpack(tkbutton(m, text = "New model", command = function() tw()))   
    tkpack(tkbutton(m, text = "Refresh", command = norm.refresh), 
        side = "left")
    tkpack(tkbutton(m, text = "Exit", command = function() tkdestroy(m)), 
        side = "right")
   
}

