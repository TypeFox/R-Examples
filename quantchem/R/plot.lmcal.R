"plot.lmcal" <-
function (x,type=c("curve","residuals","chronologic","qqplot","cook","optimization"),
trend=TRUE,confidence=TRUE,prediction=TRUE,lines=c(1,2,3,4),colors=c(1,1,1,1,1),xlab=NULL,ylab=NULL,...) 
{

obj = x

plotcurve <- function(modelname,plotname) {
plot(obj$x,obj$y,xlab=xlab,ylab=ylab,main=plotname,col=colors[1],...);
lines(obj$graph$grid,eval(parse(text=paste("obj$graph$fit$",modelname))),lty=lines[1],col=colors[1]);
if (prediction) { lines(obj$graph$grid,eval(parse(text=paste("obj$graph$upperp$",modelname))),lty=lines[2],col=colors[2]);
lines(obj$graph$grid,eval(parse(text=paste("obj$graph$lowerp$",modelname))),lty=lines[2],col=colors[2]); }
if (confidence) { lines(obj$graph$grid,eval(parse(text=paste("obj$graph$upperc$",modelname))),lty=lines[3],col=colors[3]);
lines(obj$graph$grid,eval(parse(text=paste("obj$graph$lowerc$",modelname))),lty=lines[3],col=colors[3]); }
}

plotresiduals <- function(modelname,plotname) {
res = eval(parse(text=paste("residuals(obj$models$",modelname,")")));
xx = obj$x;
plot(xx,res,col=colors[1],main=plotname,xlab=xlab,ylab=ylab); abline(h=0,lty=lines[1],col=colors[2],...);
if (trend) {
fitr = loess(res ~ xx,span=2); fitrp = predict(fitr,newdata=data.frame(xx=obj$graph$grid));
lines(obj$graph$grid,fitrp,lty=lines[2],col=colors[3]);
}
}

plotchronologic <- function(modelname,plotname) {
res = eval(parse(text=paste("residuals(obj$models$",modelname,")")));
xx = 1:length(obj$x);
plot(xx,res,col=colors[1],main=plotname,xlab=xlab,ylab=ylab); abline(h=0,lty=lines[1],col=colors[2],...);
if (trend) {
fitr = loess(res ~ xx,span=2); 
lines(xx,fitted(fitr),lty=lines[2],col=colors[3]);
}
}

plotqqplot <- function(modelname,plotname) {
res = eval(parse(text=paste("residuals(obj$models$",modelname,")")));
qqnorm(res,xlab=xlab,ylab=ylab,main=plotname,...); qqline(res,lty=lines[1],col=colors[1]);
}

plotcook <- function(modelname,plotname) {
res = eval(parse(text=paste("cooks.distance(obj$models$",modelname,")")));
xx = 1:length(obj$x);
plot(xx,res,col=colors[1],main=plotname,xlab=xlab,ylab=ylab,type="h",...);

}


    s <- match.arg(type)
    pt <- switch(s, curve = 0, residuals = 1, chronologic = 2, qqplot = 3, cook = 4, optimization = 5);

if (pt == 0) {
if (is.null(xlab)) xlab="Amount"; 
if (is.null(ylab)) ylab="Response";

plotcurve("p1","Linear");
plotcurve("p2","Quadratic");
plotcurve("p3","Cubical");
plotcurve("p4","4th Order");
plotcurve("P1","Linear weighted");
plotcurve("P2","Quadratic weighted");
plotcurve("P3","Cubical weighted");
plotcurve("P4","4th Order weighted");
plotcurve("l1","Log-Log Linear");
plotcurve("l1","Log-Log Quadratic");
plotcurve("bx","Box-Cox (x)");
plotcurve("by","Box-Cox (y)");
plotcurve("r1","Robust linear");
plotcurve("r2","Robust quadratic");
plotcurve("r3","Robust cubical");
plotcurve("r4","Robust 4th order");
plotcurve("R1","Robust linear, weighted");
plotcurve("R2","Robust quadratic, weighted");
plotcurve("R3","Robust cubical, weighted");
plotcurve("R4","Robust 4th order, weighted");

}
else if (pt == 1) {

if (is.null(xlab)) xlab="Amount"; 
if (is.null(ylab)) ylab="Residuals";

plotresiduals("p1","Linear");
plotresiduals("p2","Quadratic");
plotresiduals("p3","Cubical");
plotresiduals("p4","4th Order");
plotresiduals("P1","Linear weighted");
plotresiduals("P2","Quadratic weighted");
plotresiduals("P3","Cubical weighted");
plotresiduals("P4","4th Order weighted");
plotresiduals("l1","Log-Log Linear");
plotresiduals("l1","Log-Log Quadratic");
plotresiduals("bx","Box-Cox (x)");
plotresiduals("by","Box-Cox (y)");
plotresiduals("r1","Robust linear");
plotresiduals("r2","Robust quadratic");
plotresiduals("r3","Robust cubical");
plotresiduals("r4","Robust 4th order");
plotresiduals("R1","Robust linear, weighted");
plotresiduals("R2","Robust quadratic, weighted");
plotresiduals("R3","Robust cubical, weighted");
plotresiduals("R4","Robust 4th order, weighted");


}
else if (pt == 2) {

if (is.null(xlab)) xlab="Sample"; 
if (is.null(ylab)) ylab="Residuals";

plotchronologic("p1","Linear");
plotchronologic("p2","Quadratic");
plotchronologic("p3","Cubical");
plotchronologic("p4","4th Order");
plotchronologic("P1","Linear weighted");
plotchronologic("P2","Quadratic weighted");
plotchronologic("P3","Cubical weighted");
plotchronologic("P4","4th Order weighted");
plotchronologic("l1","Log-Log Linear");
plotchronologic("l1","Log-Log Quadratic");
plotchronologic("bx","Box-Cox (x)");
plotchronologic("by","Box-Cox (y)");
plotchronologic("r1","Robust linear");
plotchronologic("r2","Robust quadratic");
plotchronologic("r3","Robust cubical");
plotchronologic("r4","Robust 4th order");
plotchronologic("R1","Robust linear, weighted");
plotchronologic("R2","Robust quadratic, weighted");
plotchronologic("R3","Robust cubical, weighted");
plotchronologic("R4","Robust 4th order, weighted");


}
else if (pt == 3) {

if (is.null(xlab)) xlab="Theoretical Quantiles"; 
if (is.null(ylab)) ylab="Residuals Quantiles";

plotqqplot("p1","Linear");
plotqqplot("p2","Quadratic");
plotqqplot("p3","Cubical");
plotqqplot("p4","4th Order");
plotqqplot("P1","Linear weighted");
plotqqplot("P2","Quadratic weighted");
plotqqplot("P3","Cubical weighted");
plotqqplot("P4","4th Order weighted");
plotqqplot("l1","Log-Log Linear");
plotqqplot("l1","Log-Log Quadratic");
plotqqplot("bx","Box-Cox (x)");
plotqqplot("by","Box-Cox (y)");
plotqqplot("r1","Robust linear");
plotqqplot("r2","Robust quadratic");
plotqqplot("r3","Robust cubical");
plotqqplot("r4","Robust 4th order");
plotqqplot("R1","Robust linear, weighted");
plotqqplot("R2","Robust quadratic, weighted");
plotqqplot("R3","Robust cubical, weighted");
plotqqplot("R4","Robust 4th order, weighted");


}

else if (pt == 4) {

if (is.null(xlab)) xlab="Sample"; 
if (is.null(ylab)) ylab="Cook's Distance";

plotcook("p1","Linear");
plotcook("p2","Quadratic");
plotcook("p3","Cubical");
plotcook("p4","4th Order");
plotcook("P1","Linear weighted");
plotcook("P2","Quadratic weighted");
plotcook("P3","Cubical weighted");
plotcook("P4","4th Order weighted");
plotcook("l1","Log-Log Linear");
plotcook("l1","Log-Log Quadratic");
plotcook("bx","Box-Cox (x)");
plotcook("by","Box-Cox (y)");
}

else if (pt == 5) {

ylab="Mean Relative Error"; 
xlab="Gamma";
plot(mrex~gamma,obj$weigh,main="Weighing on x",xlab=xlab,ylab=ylab,type="l",col=colors[1],lty=lines[1]);
abline(v=obj$wx,lty=lines[2],col=colors[2]);
mtext(as.character(obj$wx),at=obj$wx);
plot(mrey~gamma,obj$weigh,main="Weighing on y",xlab=xlab,ylab=ylab,type="l",col=colors[1],lty=lines[1]);
abline(v=obj$wy,lty=lines[2],col=colors[2]);
mtext(as.character(obj$wy),at=obj$wy);
xlab="Lambda"; 
ylab="Log - Likelihood";
plot(x$bx$x,x$bx$y,main="Box-Cox on x",xlab=xlab,ylab=ylab,type="l",col=colors[1],lty=lines[1]);
abline(v=obj$px,lty=lines[2],col=colors[2]);
mtext(as.character(obj$px),at=obj$px);
plot(x$by$x,x$by$y,main="Box-Cox on y",xlab=xlab,ylab=ylab,type="l",col=colors[1],lty=lines[1]);
abline(v=obj$py,lty=lines[2],col=colors[2]);
mtext(as.character(obj$py),at=obj$py);

}


}

