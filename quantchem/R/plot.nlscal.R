"plot.nlscal" <-
function (x,type=c("curve","residuals","chronologic","qqplot"),
trend=TRUE,lines=c(1,2,3,4),colors=c(1,1,1,1,1),xlab=NULL,ylab=NULL,...) 
{

obj=x

plotcurve <- function(modelname,plotname) {
plot(obj$x,obj$y,xlab=xlab,ylab=ylab,main=plotname,col=colors[1],...);
lines(obj$graph$grid,eval(parse(text=paste("obj$graph$fit$",modelname))),lty=lines[1],col=colors[1]);
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

try(plotcurve("a1","Asymptotic"));
try(plotcurve("a2","Asymptotic (origin)"));
try(plotcurve("g1","Logistic"));
try(plotcurve("g2","4PL"));
try(plotcurve("m1","Michaelis-Menten"));
try(plotcurve("s1","Spline"));
}
else if (pt == 1) {

if (is.null(xlab)) xlab="Amount"; 
if (is.null(ylab)) ylab="Residuals";

try(plotresiduals("a1","Asymptotic"));
try(plotresiduals("a2","Asymptotic (origin)"));
try(plotresiduals("g1","Logistic"));
try(plotresiduals("g2","4PL"));
try(plotresiduals("m1","Michaelis-Menten"));
try(plotresiduals("s1","Spline"));

}
else if (pt == 2) {

if (is.null(xlab)) xlab="Sample"; 
if (is.null(ylab)) ylab="Residuals";

try(plotchronologic("a1","Asymptotic"));
try(plotchronologic("a2","Asymptotic (origin)"));
try(plotchronologic("g1","Logistic"));
try(plotchronologic("g2","4PL"));
try(plotchronologic("m1","Michaelis-Menten"));
try(plotchronologic("s1","Spline"));

}
else if (pt == 3) {

if (is.null(xlab)) xlab="Theoretical Quantiles"; 
if (is.null(ylab)) ylab="Residuals Quantiles";


try(plotqqplot("a1","Asymptotic"));
try(plotqqplot("a2","Asymptotic (origin)"));
try(plotqqplot("g1","Logistic"));
try(plotqqplot("g2","4PL"));
try(plotqqplot("m1","Michaelis-Menten"));
try(plotqqplot("s1","Spline"));
}
}

