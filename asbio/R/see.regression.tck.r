see.regression.tck <- function () 
{

    if (!exists("slider.env")) 
        slider.env <- NULL; suppressWarnings(rm(slider.env)); slider.env <<- new.env()# Dummy to trick R CMD check
    beta0 <- 0
    assign("beta0", tclVar(beta0), envir = slider.env)
    beta1 <- 0.6
    assign("beta1", tclVar(beta1), envir = slider.env)
    sigma <- 2
    assign("sigma", tclVar(sigma), envir = slider.env)
    n <- 5
    assign("n", tclVar(n), envir = slider.env)
    

xpts<- c(5,15)
dev.new(height = 5, width = 10)

plot_it<-function(...){

dev.hold()
        beta0 <- as.numeric(evalq(tclvalue(beta0), envir = slider.env))
        beta1 <- as.numeric(evalq(tclvalue(beta1), envir = slider.env))
        sigma <- as.numeric(evalq(tclvalue(sigma), envir = slider.env))
        n <- as.numeric(evalq(tclvalue(n), envir = slider.env))

X<-Y<-0:20; Z<-seq(0,.6,.03)
pts <- runif(n, 0, 20) 


layout(matrix(c(1,1,1,1,1,1,2,2,2,2), 1, 10, byrow = TRUE))
e<-scatterplot3d(X,Y,Z,scale.y=1.75,type="n",angle=50,xlab=expression(italic(x)),ylab=expression(italic(y)),
zlab=expression(paste(italic(f),"(",italic(y),")", sep = "")),lty.hide=0,box=FALSE, mar = c(3, 3, 0, 2))
e$points3d(seq(0,20),(beta1*seq(0,20))+beta0, rep(0,21),type="l",lwd=2) #true reg. line

for(i in 1:length(xpts)){
R1 <- seq(max(((beta1*xpts[i])+beta0)-3*sigma, -5),min(((beta1*xpts[i])+beta0)+3*sigma, 25), .1)
L1 <- length(R1)
e$points3d(rep(xpts[i],L1),R1,dnorm(R1,beta1*xpts[i]+beta0,sigma),type="l")
e$points3d(rep(xpts[i],L1),R1,rep(0,L1),type="l")
e$points3d(rep(xpts[i],L1),rep(beta1*xpts[i]+beta0,L1),seq(0,round(dnorm(beta1*xpts[i]+beta0,beta1*xpts[i]+beta0,sigma),3),length=L1),type="l")
}

epsilon <- matrix(nrow = n, ncol = 1)
loc <- matrix(nrow = n, ncol = 1)

for(i in 1:n){
loc[i] <- rnorm(1, beta1*pts[i]+beta0, sigma)
epsilon[i] <- (loc[i] - beta1*pts[i]+beta0)^2
e$points3d(pts[i],loc[i],0, pch=16)
e$points3d(rep(pts[i], 20), seq(loc[i], beta1*pts[i]+beta0, length = 20), rep(0,20),type="l",col="red",lty=2)
}

text(3, 10.75, "Population regression line", cex=1.6)
text(3, 9.5, bquote(paste(italic(E),"(",italic(Y[i]),") = ", .(beta0), " + ", .(beta1),italic(X[i]),sep="")), cex=1.5)
text(3, 8.75, bquote(paste(italic(epsilon[i]), " ~ ", italic(N), "(0, ", .(sigma^2),")", sep = "")), cex=1.5)
text(3, 8, bquote(paste(italic(E), "(", italic(MSR), ") = ", .(sigma^2 + beta1^2 * sum((pts-mean(pts))^2)), ",   ", italic(E), "(", italic(MSE), ") = ", .(sigma^2),sep = "")), cex=1.5)

reg <- lm(loc ~ pts)
par(mar = c(5,4,8,3))
p <- plotCI.reg(pts, loc, CI = FALSE, PI = FALSE, resid = TRUE, resid.col=2, resid.lty=2, pch = 16, cex=1.2, cex.lab=1.4, cex.axis=1.4, xlab = expression(italic(x)), ylab = expression(italic(y)))
b <- round(coef(p), 2)
MS <- round(anova(p)$"Mean Sq", 2)
P <- round(anova(p)$"Pr(>F)"[1], 3)

mtext(side = 3, "Sample regression line", cex=1, line = 6)
mtext(side = 3, bquote(paste(italic(Y[i])," = ", .(b[1]), " + ", .(b[2]),italic(X[i]),sep="")), cex=.9, line = 3)
mtext(side = 3, bquote(paste(italic(MSR), " = ", .(MS[1]), ",   ", italic(MSE), " = ", .(MS[2]),",   ", italic(P),"-value = ", .(P), sep = "")), cex=.9, line = 1)
dev.flush()
}


tclServiceMode(TRUE)
    m <- tktoplevel()
    tkwm.title(m, "Visualizing regression")
    tkpack(tklabel(m, text = "      Visualizing regression      "))
    tkwm.geometry(m, "+0+0")
    
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "\u03B2\u2080", font = c("Helvetica", "10", 
        "normal"), width = "10"), side = "right")
    tkpack(sc <- tkscale(fr, command = plot_it, from = -3, 
        to = 14, orient = "horiz", resolution = 0.1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = beta0), envir = slider.env)
    
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "\u03B2\u2081", font = c("Helvetica", "10", 
        "normal"), width = "10"), side = "right")
    tkpack(sc <- tkscale(fr, command = plot_it, from = -0.5, 
        to = 1, orient = "horiz", resolution = 0.1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = beta1), envir = slider.env)
    
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "\u03C3", font = c("Helvetica", "10", 
        "normal"), width = "10"), side = "right")
    tkpack(sc <- tkscale(fr, command = plot_it, from = 0.1, 
        to = 4, orient = "horiz", resolution = 0.1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = sigma), envir = slider.env)
  
   tkpack(fr <- tkframe(m), side = "top")
    tkpack(tklabel(fr, text = "n", font = c("Helvetica", "10", 
        "normal"), width = "10"), side = "right")
    tkpack(sc <- tkscale(fr, command = plot_it, from = 2, 
        to = 20, orient = "horiz", resolution = 1, showvalue = TRUE), 
        side = "left")
    assign("sc", sc, envir = slider.env)
    evalq(tkconfigure(sc, variable = n), envir = slider.env) 
    
    
    tkpack(tkbutton(m, text = "Refresh", command = plot_it), 
        side = "left")
    tkpack(tkbutton(m, text = "Exit", command = function() tkdestroy(m)), 
        side = "right")

        
}
