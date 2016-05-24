plot.criterionRkh <-
function(x,choice.H,choice.K,...){
if (!inherits(x, "criterionRkh")) stop("use only with \"criterionRkh\" objects")
method <- x$method
H <- x$H
K <- x$K
p <- dim(x$Rkhbootmean)[2]
n <- x$n
Met <- x$method

# plot for different values of h
if(missing(choice.H)){cat("none value of h for ploting criterion for several values of h","\n")}
else{if(sum(!(choice.H %in% H))>0){
cat("the values of h must be included in : ",x$H)
cat("\n")}
else{
for (i in choice.H){
dev.new()
par(mfrow=c(2,1),mar=c(5,5,4,2))
laby<-expression(hat(R)[paste("k,h")]^paste("(b)"))
ind<-which(H==i)

titre1<-paste(Met,": h=",H[ind],", n=",n)
boxplot(x$Rkhboot[[ind]],ylab=laby,xlab="k",ylim=c(0,1),main=titre1,axes=F)
axis(1,K)
axis(2,seq(0,1,by=0.2))

laby<-expression(hat(R)[paste("k,h")])
plot(K,x$Rkhbootmean[ind,],type="o",ylim=c(0,1),xlab="k",ylab=laby,main=titre1,axes=F)
axis(1,K)
axis(2,seq(0,1,by=0.2))}
}
}

# plot for different values of k
if(missing(choice.K)){cat("none values of k for ploting criterion for several value of k","\n")}
else{if(sum(!(choice.K %in% K))>0){
cat("the values of k must be included in : ",x$K)
cat("\n")}
else{
for (i in choice.K){
dev.new()
par(mfrow=c(2,1),mar=c(5,5,4,2))
laby<-expression(hat(R)[paste("k,h")]^paste("(b)"))

crit<-NULL
titre1<-paste(Met,": k=",i,", n=",n)
for(j in 1:length(H)){
crit<-cbind(crit,x$Rkhboot[[j]][,i])
}

boxplot(crit,ylab=laby,xlab="h",ylim=c(0,1),main=titre1,axes=F)
axis(1,1:length(H),H)
axis(2,seq(0,1,by=0.2))

laby<-expression(hat(R)[paste("k,h")])
titre1<-paste(Met,": k=",i,", n=",n)
plot(1:length(H),x$Rkhbootmean[,i],type="o",ylim=c(0,1),xlab="h",ylab=laby,main=titre1,axes=F)
axis(1,1:length(H),H)
axis(2,seq(0,1,by=0.2))}
}
}


open3d()
persp3d(x$K,x$H,t(x$Rkhbootmean),theta=10,phi=20,zlim=c(0,1),xlab="k",ylab="h",zlab="R(k,h)",box=T,front="lines",col = "red",axes=F)
bbox3d(xat=x$K, yat=x$H, zat=c(0.2,0.4,0.8,1),color=c("white","black"),emission="#333377",specular="#3333FF")  
title3d(Met,line=4,cex=1.5)
}

