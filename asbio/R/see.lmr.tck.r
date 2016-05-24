pm1<-function(Y,X,sz=1,showXY=TRUE){
X<-matrix(nrow=length(Y),data=X)
Y<-as.matrix(Y)
n<-length(Y)
p<-ncol(X)
b<-solve(t(X)%*%X)%*%t(X)%*%Y
b1<-round(b,8)
MSR<-round(((t(b)%*%t(X)%*%Y)-(1/n)*t(Y)%*%matrix(1,nrow=n,ncol=n)%*%Y)/(p-1),10)
MSE<-round((t(Y)%*%(Y)-t(b)%*%t(X)%*%Y)/(n-p),10)
Yhat<-round(X%*%b,8)

par(mar=c(0,0,0,0))
layout(matrix(c(rep(1,1),rep(2,2),rep(3,1),rep(4,1),rep(5,1)),2,3, byrow = TRUE))
plot(seq(1,10),seq(1,10),type="n",xlab="",ylab="",xaxt="n",yaxt="n",bty="n")
if(showXY==TRUE){
legend("center",ncol=1,legend=Y,bty="n",title="Y",cex=1.3*sz)}
plot(seq(1,10),seq(1,10),type="n",xlab="",ylab="",xaxt="n",yaxt="n",bty="n")
if(showXY==TRUE)
legend("center",ncol=ncol(X),legend=c(as.vector(X)),bty="n",title="X",cex=1.3*sz)

plot(seq(1,10),seq(1,10),type="n",xlab="",ylab="",xaxt="n",yaxt="n")
text(5,9.5,"Coefficients",cex=1.6*sz)
text(1.8,8.1,expression(paste(hat(beta)," = ")),cex=1.4*sz)
text(4.1,8,paste("(X'X","\u0029","\u02c9","\u00b9","X'Y = ",sep=""),cex=1.4*sz)
legend(5,9,b1,cex=1.3*sz,bty="n")
     
plot(seq(1,10),seq(1,10),type="n",xlab="",ylab="",xaxt="n",yaxt="n")
text(5,9.5,"Mean squares",cex=1.6*sz)
text(1.8,8.1,expression(paste("MSR = ", hat(beta),"'","X'Y - ",frac(1,n)," Y'1Y/(p - 1)")),cex=1.4*sz,adj=0)
text(1.8,6.8,bquote(paste("       = ",.(MSR))),cex=1.4*sz,adj=0)
text(1.8,5.5,expression(paste("MSE = Y'",hat(beta),"'X'Y/(n - p)")),cex=1.4*sz,adj=0)
text(1.8,4.1,paste("       = ",bquote(.(MSE))),cex=1.4*sz,adj = 0)

plot(seq(1,10),seq(1,10),type="n",xlab="",ylab="",xaxt="n",yaxt="n")
text(5,9.5,"Fitted values",cex=1.6*sz)
text(1.8,8.1,expression(paste(hat(Y)," = X", hat(beta)," = ")),cex=1.4*sz,adj=0)
legend(5,9,Yhat,cex=1.3*sz,bty="n")
}

see.lmr.tck<-function (){

local({
have_ttk <- as.character(tcl("info", "tclversion")) >= "8.5"
if (have_ttk) {
    tkbutton <- ttkbutton
    tkcheckbutton <- ttkcheckbutton
    tkentry <- ttkentry
    tkframe <- ttkframe
    tklabel <- ttklabel
    tkradiobutton <- ttkradiobutton
}
tclServiceMode(FALSE) 
dialog.sd <- function() {
  tt <- tktoplevel()
  tkwm.title(tt,"Linear models")
  y.entry <- tkentry(tt, textvariable=Y, width =45)
  x0.entry <- tkentry(tt, textvariable=X0, width =45) 
  x1.entry <- tkentry(tt, textvariable=X1, width =45) 
  x2.entry <- tkentry(tt, textvariable=X2, width =45) 
  x3.entry <- tkentry(tt, textvariable=X3, width =45) 
  x4.entry <- tkentry(tt, textvariable=X4, width =45) 
  nvars <- tkentry(tt, textvariable=NV, width =5) 
showXY<-tclVar(1) 
  done <- tclVar(0)
 
reset<-function(){
Y<-"c(20,30,10,15,5,45,60,55,45)"
X0<-"c(rep(1,9))"
X1<-"c(13,20,10,11,2,25,30,25,23)"
X2<-"c(1.2,2,1.5,1,0.3,2,3,2.7,2.5)"
X3<-"c(15,14,16,12,10,18,25,24,20)"
X4<-"(45,120,100,56,5,20,5,15,15)"
NV<-"4"
  tclvalue(showXY)<-"1"
}

reset.but <- tkbutton(tt, text = "Reset", command = reset)
submit.but <- tkbutton(tt, text = "Submit", command = function() tclvalue(done) <- 1)



build <- function() {
  showXY <- as.logical(tclObj(showXY))
  Y <-parse(text=tclvalue(Y))[[1]]
  NV<-parse(text=tclvalue(NV))[[1]]
  
if(NV == "2"){
X0 <- parse(text=tclvalue(X0))[[1]]
X1 <- parse(text=tclvalue(X1))[[1]]
X<-substitute(cbind(as.numeric(X0),as.numeric(X1)))}  
if(NV == "3"){
X0 <- parse(text=tclvalue(X0))[[1]]
X1 <- parse(text=tclvalue(X1))[[1]]
X2 <- parse(text=tclvalue(X2))[[1]]
X<-substitute(cbind(as.numeric(X0),as.numeric(X1),as.numeric(X2)))}
if(NV =="4"){
X0 <- parse(text=tclvalue(X0))[[1]]
X1 <- parse(text=tclvalue(X1))[[1]]
X2 <- parse(text=tclvalue(X2))[[1]]
X3 <- parse(text=tclvalue(X3))[[1]]
X<-substitute(cbind(as.numeric(X0),as.numeric(X1),as.numeric(X2),as.numeric(X3)))}
if(NV=="5"){
X0 <- parse(text=tclvalue(X0))[[1]]
X1 <- parse(text=tclvalue(X1))[[1]]
X2 <- parse(text=tclvalue(X2))[[1]]
X3 <- parse(text=tclvalue(X3))[[1]]
X4 <- parse(text=tclvalue(X4))[[1]]
X<-substitute(cbind(as.numeric(X0),as.numeric(X1),as.numeric(X2),as.numeric(X3),as.numeric(X4)))}
substitute(pm1(as.numeric(Y),as.numeric(X),showXY=showXY))
}                

nc.cbut <- tkcheckbutton(tt, text="Show XY", variable=showXY)
  tkgrid(tklabel(tt, text = "Linear model (Regression)"), 
      columnspan = 2)
  tkgrid(tklabel(tt, text = ""))
  tkgrid(tklabel(tt, text = 'Y'), y.entry)
  tkgrid(tklabel(tt, text = 'X0'), x0.entry)
  tkgrid(tklabel(tt, text = 'X1'), x1.entry)
  tkgrid(tklabel(tt, text = "X2"), x2.entry)
  tkgrid(tklabel(tt, text = "X3"), x3.entry)
  tkgrid(tklabel(tt, text = "X4"), x4.entry)
  tkgrid(tklabel(tt, text = "Number of X vars."), nvars)
  tkgrid(tklabel(tt, text = ""))
  tkgrid(nc.cbut)  
  tkgrid(tklabel(tt, text = ""))
  tkgrid(submit.but, reset.but, sticky ="e")
  
  tkbind(tt, "<Destroy>", function() tclvalue(done) <- 2)
  tkwait.variable(done)
  if (tclvalue(done) == "2") 
      stop("aborted")
  tkdestroy(tt)
  cmd <- build()
  eval.parent(cmd)
}
Y<-tclVar("c(20,30,10,15,5,45,60,55,45)")
X0<-tclVar("c(rep(1,9))")
X1<-tclVar("c(13,20,10,11,2,25,30,25,23)")
X2<-tclVar("c(1.2,2,1.5,1,0.3,2,3,2.7,2.5)")
X3<-tclVar("c(15,14,16,12,10,18,25,24,20)")
X4<-tclVar("c(45,120,100,56,5,20,5,15,15)")
NV<-tclVar("5")
dialog.sd()
})
}    
                                                            