see.HW<-function(parg){
p<-seq(0,1,length=1000);pp<-p^2;q<-1-p;qq<-q^2;twoqp<-2*q*p
pn<-parg;qn<-round(1-pn,3);ppn<-round(pn^2,3);qqn<-round(qn^2,3);twoqpn<-signif(2*qn*pn,3)
dev.hold()
layout(matrix(c(1,2,rep(3,8)),10,1,byrow=TRUE))
par(mar=c(0,0,0,0))
plot(seq(0,1),seq(0,1),xaxt="n",yaxt="n",xlab="",ylab="",bty="n",type="n")
legend("center",legend="Hardy-Weinberg genotypic proportions  ",cex=1.6,bty="n")
plot(seq(0,1),seq(0,1),xaxt="n",yaxt="n",xlab="",ylab="",bty="n",type="n")
legend("center",pch=22,pt.bg=c("coral","slateblue1","salmon4"),pt.cex=3,legend=c(paste("pp = ", ppn),paste("pq = ",twoqpn), paste("qq = ",qqn)),ncol=3,cex=1.5,bty="n",yjust=1)
par(mar=c(4.5,4,0,.6))
fig<-plot(seq(0,1,length=1000),seq(0,1,length=1000),xlab=expression(italic(p)), ylab = "Genotype proportions",type="n")
grid(fig)
points(seq(0,1,length=1000),pp,col="coral",type="l",lwd=2)
points(seq(0,1,length=1000),twoqp,col="slateblue1","l",lwd=2)
points(seq(0,1,length=1000),qq,col="salmon4","l",lwd=2)
abline(v=parg,lwd=2)
dev.flush()
}

see.HW.tck<-function (){

if (!exists("slider.envir")) 
slider.env <- NULL; suppressWarnings(rm(slider.env)); slider.env <<- new.env()# Dummy to trick R CMD check 
p <- 0.5
assign("p", tclVar(p), envir = slider.env)

refresh <- function(...) {
p <- as.numeric(evalq(tclvalue(p), envir = slider.env))
see.HW(p)
}
tclServiceMode(TRUE)
m <- tktoplevel()
tkwm.title(m, "the Hardy Weinberg Equilibrium")
tkwm.geometry(m, "+0+0")
tkpack(tklabel(m,text="      Visualizing the Hardy Weinberg Equilibrium      "))
tkpack(fr <- tkframe(m), side = "top")
tkpack(tklabel(fr, text = "p  ", font=c("Helvetica","10","italic")),side="left", anchor = "s")
tkpack(sc <- tkscale(fr, command = refresh, from = 0, 
    to = 1, orient = "horiz", resolution = .01, showvalue = TRUE), 
    side = "left", anchor="n")
assign("sc", sc, envir = slider.env)
evalq(tkconfigure(sc, variable = p), envir = slider.env)
tkpack(tkbutton(m, text = "Exit", command = function() tkdestroy(m)))
 }