pm<-function(Y, X1, X2, X1X2, change.order=FALSE, delete = 0){

if(delete != 0){
Y <- Y[-delete]
X1 <- X1[-delete]
X2 <- X2[-delete]
X1X2 <- X1X2[-delete]
}

n <- tapply(Y,X1,length)
status <- as.character(ifelse(n[1]==n[2], "Ballanced", "Unballanced"))

par(mar=c(0,0,0,0))
layout(matrix(c(rep(1,1),rep(2,2),rep(3,1),rep(4,1),rep(5,1)),2,3, byrow = TRUE))
plot(seq(1,10),seq(1,10),type="n",xlab="",ylab="",xaxt="n",yaxt="n",bty="n")
text(2.5,5,"Y = ",font=2,cex=2)
segments(c(4.5,7),c(2,2),c(4.5,7),c(9,9),lwd=2)
segments(c(4.5,7),c(2,9),c(4.7,6.8),c(2,9),lwd=2)
segments(c(4.5,7),c(9,2),c(4.7,6.8),c(9,2),lwd=2)

legend("center",ncol=1,legend=Y,bty="n",cex=.9)
plot(seq(1,10),seq(1,10),type="n",xlab="",ylab="",xaxt="n",yaxt="n",bty="n")
text(4,5, "X = ", font=2,cex=2)
segments(c(4.5,6.6),c(2,2),c(4.5,6.6),c(9,9),lwd=2)
segments(c(4.5,6.6),c(2,9),c(4.6,6.5),c(2,9),lwd=2)
segments(c(4.5,6.6),c(9,2),c(4.6,6.5),c(9,2),lwd=2)
legend("center",ncol=3,legend=c(X1,X2,X1X2),bty="n",cex=.9)
if(status == "Unballanced") legend("center", yjust = 3, legend = status, text.col = 2, bty="n", cex = 1.2)
if(delete != 0 & change.order == FALSE) legend("topright", legend = paste("Obs. no. ", delete, " deleted"), bty="n", cex = .9)
if(delete != 0 & change.order == TRUE) legend("topright", legend = c(paste("Obs. no. ", delete, " deleted"), "Predictor order reversed"), bty="n", cex = .9)
if(delete == 0 & change.order == TRUE) legend("topright", legend = "Predictor order reversed", bty="n", cex = .9)
#SS

 
 SSE.mu<-anova(lm(Y~1))$"Sum Sq"[1]
 SSE.X1<-anova(lm(Y~X1))$"Sum Sq"[2]
 SSE.X2<-anova(lm(Y~X2))$"Sum Sq"[2]
 SSE.noX1<-anova(lm(Y~X2+X1X2))$"Sum Sq"[3]
 SSE.noX2<-anova(lm(Y~X1+X1X2))$"Sum Sq"[3]
 SSE.X1pX2<-anova(lm(Y~X1+X2))$"Sum Sq"[3]
 SSE.X1iX2<-anova(lm(Y~X1+X2+X1X2))$"Sum Sq"[4]

 if(change.order==FALSE)anovaI<-anova(lm(Y~X1+X2+X1X2))
 if(change.order==TRUE)anovaI<-anova(lm(Y~X2+X1+X1X2))
 Dfs <- anovaI$Df
 Dfi <- Dfs[3]
 Dfe <- Dfs[4]
 MSE <- tail(anovaI$"Mean Sq", 1)

 ##type I
 
 plot(seq(1,10),seq(1,10),type="n",xlab="",ylab="",xaxt="n",yaxt="n",bty="n")
 rect(1,1,10,10,col="wheat")
legend("top",inset=.2,title=ifelse(change.order==FALSE,"Type I Sum of Squares\n\nModel: Y = X1 + X2 + X1:X2","Type I Sum of Squares\n\nModel: Y = X2 + X1 + X1:X2"),legend=c(
 ifelse(change.order==FALSE,paste("SSE.X1 = ",round(SSE.mu-SSE.X1,3)),paste("SSE.X2 = ",round(SSE.mu-SSE.X2,3))), ifelse(change.order==FALSE,paste("SSE.X2 = ",round(SSE.X1-SSE.X1pX2,3)),paste("SSE.X1 = ",round(SSE.X2-SSE.X1pX2,3))),paste("SSE.X1:X2 = ",round(SSE.X1pX2-SSE.X1iX2,3)),paste("SSE = ",round(SSE.X1iX2,3))),bty="n")

legend("bottom",title="P-value",legend=c(ifelse(change.order==FALSE,paste("X1 = ",round(anovaI$"Pr(>F)"[1],6)),paste("X2 = ",round(anovaI$"Pr(>F)"[1],6))),ifelse(change.order==FALSE,(paste("X2 = ",round(anovaI$"Pr(>F)"[2],6))),paste("X1 = ",round(anovaI$"Pr(>F)"[2],6))),paste("X1:X2 = ",round(anovaI$"Pr(>F)"[3],6))),bty="n",inset=.2)

 
#type II 
SSE.Xn1 <- SSE.X2 - SSE.X1pX2
SSE.Xn2 <- SSE.X1  -SSE.X1pX2
SSE.X1X2 <- SSE.X1pX2-SSE.X1iX2
SSE <- SSE.X1iX2

plot(seq(1,10),seq(1,10),type="n",xlab="",ylab="",xaxt="n",yaxt="n",bty="n")
rect(1,1,10,10,col="wheat3") 

legend("top",inset=.2,title=ifelse(change.order==FALSE,"Type II Sum of Squares\n\nModel: Y = X1 + X2 + X1:X2","Type II Sum of Squares\n\nModel: Y = X2 + X1 + X1:X2"),legend=c(
ifelse(change.order==FALSE,paste("SSE.X1 = ",round(SSE.Xn1,3)),paste("SSE.X2 = ",round(SSE.Xn2, 3))), ifelse(change.order==FALSE,paste("SSE.X2 = ",round(SSE.Xn2,3)),
paste("SSE.X1 = ",round(SSE.Xn1,3))),paste("SSE.X1:X2 = ",round(SSE.X1X2,3)),paste("SSE = ",round(SSE.X1iX2,3))),bty="n")

if(change.order==FALSE)P.val <- c(pf(SSE.Xn1/MSE, 1, Dfe, lower.tail = FALSE), pf(SSE.Xn2/MSE, 1, Dfe, lower.tail = FALSE), pf(SSE.X1X2/MSE, Dfi, Dfe, lower.tail = FALSE))   
if(change.order== TRUE)P.val <- c(pf(SSE.Xn2/MSE, 1, Dfe, lower.tail = FALSE), pf(SSE.Xn1/MSE, 1, Dfe, lower.tail = FALSE), pf(SSE.X1X2/MSE, Dfi, Dfe, lower.tail = FALSE))   

legend("bottom", title="P-value", legend=c(ifelse(change.order == FALSE, paste("X1 = ",round(P.val[1],6)), paste("X2 = ",round(P.val[1],6))), 
ifelse(change.order == FALSE, paste("X2 = ",round(P.val[2],6)), paste("X1 = ",round(P.val[2],6))),
paste("X1:X2 = ",round(P.val[3],6))),bty="n",inset=.2)

#type III
SSE.X1 <- SSE.noX1-SSE.X1iX2
SSE.X2 <- SSE.noX2-SSE.X1iX2
SSE.X1X2 <- SSE.X1pX2-SSE.X1iX2
SSE <- SSE.X1iX2

plot(seq(1,10),seq(1,10),type="n",xlab="",ylab="",xaxt="n",yaxt="n",bty="n")
rect(1,1,10,10,col="wheat4") 
legend("top",inset=.2,title=ifelse(change.order==FALSE,"Type III Sum of Squares\n\nModel: Y = X1 + X2 + X1:X2","Type III Sum of Squares\n\nModel: Y = X2 + X1 + X1:X2"),legend=c(
 ifelse(change.order==FALSE,paste("SSE.X1 = ",round(SSE.X1,3)),paste("SSE.X2 = ",round(SSE.X2,3))), ifelse(change.order==FALSE,paste("SSE.X2 = ",round(SSE.X2,3)),paste("SSE.X1 = ",round(SSE.X1,3))),paste("SSE.X1:X2 = ",round(SSE.X1pX2-SSE.X1iX2,3)),paste("SSE = ",round(SSE.X1iX2,3))),bty="n")

if(change.order==FALSE)P.val <- c(pf(SSE.X1/MSE, 1, Dfe, lower.tail = FALSE), pf(SSE.X2/MSE, 1, Dfe, lower.tail = FALSE), pf(SSE.X1X2/MSE, Dfi, Dfe, lower.tail = FALSE))   
if(change.order== TRUE)P.val <- c(pf(SSE.X2/MSE, 1, Dfe, lower.tail = FALSE), pf(SSE.X1/MSE, 1, Dfe, lower.tail = FALSE), pf(SSE.X1X2/MSE, Dfi, Dfe, lower.tail = FALSE))   

legend("bottom", title="P-value", legend=c(ifelse(change.order == FALSE, paste("X1 = ",round(P.val[1],6)), paste("X2 = ",round(P.val[1],6))), ifelse(change.order == FALSE, paste("X2 = ",round(P.val[2],6)), paste("X1 = ",round(P.val[2],6))),
paste("X1:X2 = ",round(P.val[3],6))),bty="n",inset=.2)


}






see.lmu.tck<-function (){

tclRequire("BWidget")
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
  tkwm.title(tt,"Type I, II, and III SS")
  y.entry <- tkentry(tt, textvariable=Y, width =45)
  x1.entry <- tkentry(tt, textvariable=X1, width =45) 
  x2.entry <- tkentry(tt, textvariable=X2, width =45) 
  x1x2.entry <- tkentry(tt, textvariable=X1X2, width =45) 

Delete <-tclVar("0")
 
  
 
  r.order<-tclVar(0)
  done <- tclVar(0)

submit.but <- tkbutton(tt, text = "Submit", command = function() tclvalue(done) <- 1)

build <- function() {
  Y <-parse(text=tclvalue(Y))[[1]]
  X1 <- parse(text=tclvalue(X1))[[1]]
  X2 <- parse(text=tclvalue(X2))[[1]]
  X1X2 <- parse(text=tclvalue(X1X2))[[1]]
  r.order <- as.logical(tclObj(r.order))
  delete <- tclvalue(Delete)
X<-substitute(cbind(as.numeric(X1),as.numeric(X2),as.numeric(X1X2)))
substitute(pm(as.numeric(Y),as.numeric(X1),as.numeric(X2),as.numeric(X1X2),change.order = r.order, delete = as.numeric(delete)))
}                

 
  tkgrid(tklabel(tt, text = " ANOVA linear model "), 
      columnspan = 2)
  tkgrid(tklabel(tt, text = ""))
  tkgrid(tklabel(tt, text = 'Y'), y.entry)
  tkgrid(tklabel(tt, text = 'X1'), x1.entry)
  tkgrid(tklabel(tt, text = "X2"), x2.entry)
  tkgrid(tklabel(tt, text = "X1X2"), x1x2.entry)
  tkgrid(tklabel(tt, text = ""))

  val <- seq(0, length(eval(parse(text=tclvalue(Y))[[1]])))
  val <- as.character(val)
  comboBox <- tkwidget(tt,"ComboBox", editable=FALSE, values=val,textvariable = Delete, width = 3)
  nc.cbut <- tkcheckbutton(tt, text="Reverse X1 X2 order", variable=r.order)
  tkgrid(nc.cbut)
  tkgrid(tklabel(tt, text = "Delete obs.")) 
  tkgrid(comboBox) 
  tkgrid(tklabel(tt, text = ""), submit.but)
  
  tkbind(tt, "<Destroy>", function() tclvalue(done) <- 2)
  tkwait.variable(done)
  if (tclvalue(done) == "2") 
      stop("aborted")
  tkdestroy(tt)
  cmd <- build()
  eval.parent(cmd)
}
Y<-tclVar("c(16.5,18.4,12.7,14.0,12.8,14.5,11.0,10.8,14.3,10.0,39.1,26.2,21.3,35.8,40.2,32.0,23.8,28.8,25.0,29.3)")
X1<-tclVar("c(1,1,1,1,1,-1,-1,-1,-1,-1,1,1,1,1,1,-1,-1,-1,-1,-1)")
X2<-tclVar("c(1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1)")
X1X2<-tclVar("c(1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1)")
dialog.sd()
})
}                         
