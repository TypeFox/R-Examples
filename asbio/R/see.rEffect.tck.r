r.eff <- function(){
varA = 1
dev.hold()
layout(matrix(c(1,1,1,1,2,3,3,3,3,4,5,5,5,5,6),3,5,byrow=TRUE))
par(mar=c(2.8,2,1,1),oma=c(0.5,0,0,0))
plot(seq(-6,6,.001),dnorm(seq(-6,6,.001),mean=0,sd=sqrt(varA)),ylab="",cex=.1,xaxt="n",yaxt="n")
abline(v=0,lty=2)
mtext("Density", side = 2, line =.6,cex=.9)
mtext(expression(paste(italic(E),"(",italic(Y),"|",italic(A[i]),")")), side = 1, line = 1.2,cex=.9)
A1 <- rnorm(1,0,sqrt(varA))
A2 <- rnorm(1,0,sqrt(varA))
abline(v=A1,lty=2,col="gray")
abline(v=A2,lty=2,col="gray")
text(-.25,.5*dnorm(0,0,sd=sqrt(varA)), expression(mu))


plot(seq(1,20),xaxt="n",yaxt="n",col="white",bty="n")
text(8,12.5, expression(paste("Mean = ",mu)),cex=1.5)
text(8,8.5, expression(paste("Var = ",sigma[A]^2)) ,cex=1.5)


plot(seq(-6,6,.001),dnorm(seq(-6,6,.001),mean=A1,sd=sqrt(varA*4)),ylab="",cex=.1,xaxt="n",yaxt="n")
abline(v=A1,lty=2,col="gray")
mtext("Density", side = 2, line = .6,cex=.9)
mtext(expression(paste("Yield given fertilizer 1; i.e.", italic(Y),"|",italic(A)[1])), side = 1, line = 1.2,cex=.9)
text(A1-.25,.5*dnorm(A1,A1,sd=sqrt(varA*4)), expression(mu[1]))
s<-rnorm(3,A1, sd=sqrt(varA*2))
points(s[1],.001)
text(s[1],.02, expression(Y[11]))
points(s[2],.001)
text(s[2],.02, expression(Y[12]))
points(s[3],.001)
text(s[3],.02, expression(Y[13]))

plot(seq(1,20),xaxt="n",yaxt="n",col="white",bty="n")
text(8-.25,12.5, expression(paste("Mean = ",mu[1])),cex=1.5)
text(8,8.5, expression(paste("Var = ",sigma^2)),cex=1.5)


plot(seq(-6,6,.001),dnorm(seq(-6,6,.001),mean=A2,sd=sqrt(varA*4)),ylab="",cex=.1,xaxt="n",yaxt="n")
text(A2-.25,.5*dnorm(A2,A2,sd=sqrt(varA*4)), expression(mu[2]))
abline(v=A2,lty=2,col="gray")
mtext("Density", side = 2, line = .6,cex=.9)
mtext(expression(paste("Yield given fertilizer 2; i.e.", italic(Y),"|",italic(A)[2])), side = 1,  line = 1.2,cex=.9)
s<-rnorm(3,A2, sd=sqrt(varA*2))
points(s[1],.001)
text(s[1],.02, expression(Y[11]))
points(s[2],.001)
text(s[2],.02, expression(Y[12]))
points(s[3],.001)
text(s[3],.02, expression(Y[13]))

plot(seq(1,20),xaxt="n",yaxt="n",col="white",bty="n")
text(8,12.5, expression(paste("Mean = ",mu[2])),cex=1.5)
text(8,8.5, expression(paste("Var = ",sigma^2)),cex=1.5)
dev.flush()
}



see.rEffect.tck <- function () 
{

    if (!exists("slider.env")) 
    slider.env <- NULL; suppressWarnings(rm(slider.env)); slider.env <<- new.env()# Dummy to trick R CMD check
    sigmaA <- 1
    assign("sigmaA", tclVar(sigmaA), envir = slider.env)
    
    norm.refresh <- function(...) {
          r.eff()
    }
    
    tclServiceMode(TRUE)
    m <- tktoplevel()
    tkwm.title(m, "Visualizing random effects")
    tkpack(tklabel(m, text = " Visualizing random effects "))
    tkwm.geometry(m, "+0+0")
   
    tkpack(fr <- tkframe(m), side = "top")
    tkpack(tkbutton(fr, text = "Sample", command = norm.refresh), 
        side = "left")
    tkpack(tkbutton(m, text = "Exit", command = function() tkdestroy(m)), 
        side = "right")
}
