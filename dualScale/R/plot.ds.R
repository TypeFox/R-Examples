plot.ds <-
function(x,dim1=1,dim2=2, type="Asy1", ...)
{
	# BASIC INFORMATION
	#
	labels<-c(x$ItONa,x$SubNa)
	labelsS<-c(x$SubNa)
	labelsO<-c(x$ItONa)
	cat("\nPrior Analysis: ") 
	print(x$Call)
	
if(x$tipo=="O"){ #Ordinary MC Dual Scaling
  StanOp<-as.matrix(x$Norm.Op_O) #to adapt Antonio improvements...
  AdjstOp<-as.matrix(x$Proj.Op_O)
  StanS<-as.matrix(x$Norm.Su_O)
  AdjstS<-as.matrix(x$Proj.Su_O)
	results<-x$Out_O
	color<-rep(1:x$N.Item,times=x$NoOpI)
	delta1<-round(x$Out_O[dim1,5],1)
	delta2<-round(x$Out_O[dim2,5],1)
	xlab = paste("Component ", dim1, " (delta: ", delta1,"%)", sep="")
	ylab = paste("Component ", dim2, " (delta: ", delta2,"%)", sep="")  
	}
#
if(x$tipo=="A"){ #Force Classification DS
  StanOp<-x$Norm.Op_A
  AdjstOp<-x$Proj.Op_A
  StanS<-x$Norm.Su_A
  AdjstS<-x$Proj.Su_A #### ojo
	results<-x$Out_A
	color<-rep(1:x$N.Item,times=x$NoOpI)
	delta1<-round(x$Out_A[dim1,4],1)
	delta2<-round(x$Out_A[dim2,4],1)
	xlab = paste("Component ", dim1, " (delta: ", delta1,"%)", sep="")
	ylab = paste("Component ", dim2, " (delta: ", delta2,"%)", sep="")  
	}

	if(type=="Sym"){Type<-c("Symmetric Plot")}
	if(type=="Asy1"){Type<-c("Asymmetric Plot I")}
	if(type=="Asy2"){Type<-c("Asymmetric Plot II")}
	if(type=="Sub"){Type<-c("Only Subjects")}
	if(type=="Ite"){Type<-c("Only Items Options")}
	
cat("Type of Graph:", Type, "of dimensions", print(c(dim1,dim2)));
cat("\nCumulative Delta: ", delta1+delta2,"%\n",sep="")
#
####control
#
if(ncol(StanOp)<max(dim1,dim2))
{stop('dsHELP: It is impossible to plot dim ', max(dim1, dim2))}
#
##
switch(type, 
Sym={
	X<-rbind(StanOp,StanS)
	limits<-c(floor(min(X[,dim1],X[,dim2])),ceiling(max(X[,dim1],X[,dim2])))
	aco1<-round(acos(results$SingValue[dim1])*360/(2*pi),2)
	aco2<-round(acos(results$SingValue[dim2])*360/(2*pi),2)
	plot(StanOp[,dim1],StanOp[,dim2], type = "p", col=color, pch=c(19), xlim = limits, ylim = limits, log="", main = "Symmetric Plot:\nproj.op weights and norm.su scores", sub = "\nprojected options weights and normed subject scores", 
	xlab = paste("Component ", dim1, " (row-column discrepancy: ",aco1,"o)", sep=""), 
	ylab = paste("Component ", dim2, " (row-column discrepancy: ",aco2,"o)", sep=""), 
	ann = par("ann"))
	points(StanS[,dim1],StanS[,dim2], type = "p", pch=c(3))
	text(X[,dim1]+0.01,X[,dim2]+0.1, labels,cex=0.6)},
#
#
Asy1={#### this is the default plot####
	X<-rbind(AdjstOp,StanS)
	limits<-c(floor(min(X[,dim1],X[,dim2])),ceiling(max(X[,dim1],X[,dim2])))
	#
	plot(AdjstOp[,dim1], AdjstOp[,dim2], type = "p", col=color, pch=c(19), xlim = limits, ylim = limits, log="", main = "Asymmetric Plot I:\nproj.op weights and norm.su scores", sub = NULL, 
	xlab, ylab, 
	ann = par("ann"))
	points(StanS[,dim1],StanS[,dim2], type = "p", pch=c(3))
	text(X[,dim1]+0.01,X[,dim2]+0.1, labels,cex=0.6)},
#
#
Asy2={
	X<-rbind(StanOp,AdjstS)
	limits<-c(floor(min(X[,dim1],X[,dim2])),ceiling(max(X[,dim1],X[,dim2])))
	#
	plot(StanOp[,dim1], StanOp[,dim2], type = "p", col=color, pch=c(19), xlim = limits, ylim = limits, log="", main = "Asymmetric Plot II\nproj.op weights and proj.su scores", sub = NULL, 
	#xlab = paste("Component ", dim1), 
	xlab, ylab, 
	ann = par("ann"))
	points(AdjstS[,dim1],AdjstS[,dim2], type = "p", pch=c(3))
	text(X[,dim1]+0.01,X[,dim2]+0.1, labels,cex=0.6)},
#
#
Sub={
	X<-AdjstS
  limits<-c(floor(min(X[,dim1],X[,dim2])),ceiling(max(X[,dim1],X[,dim2])))
 plot(AdjstS[,dim1], AdjstS[,dim2], type = "p", pch=c(3), xlim = limits, ylim = limits, 
       log="", main = "Subjects Plot\nProjected Scores", sub = NULL, 
	xlab, ylab, 
	ann = par("ann"))
	text(X[,dim1]+0.01,X[,dim2]+0.1, labelsS,cex=0.6)},
#
#
Ite={
	X<-rbind(AdjstOp)
	limits<-c(floor(min(X[,dim1],X[,dim2])),ceiling(max(X[,dim1],X[,dim2])))
	#
	plot(AdjstOp[,dim1], AdjstOp[,dim2], type = "p", col=color, pch=c(19), xlim = limits, ylim = limits, log="", main = "Options Plot: \nProjected Weights", sub = NULL, 
	xlab, ylab, 
	ann = par("ann"))
	text(X[,dim1]+0.01,X[,dim2]+0.1, labelsO,cex=0.6)},
#	
)
#
#
}
