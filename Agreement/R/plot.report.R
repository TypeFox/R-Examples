plot.report <-
function(x,...){

refx <- function(y){
	x=y;
}

X <- x$Data$x;
Y <- x$Data$y;
min <- x$Data$min;
max <- x$Data$max;
range <- max-min;

if (x$Data$error=="prop"){
    Y <- exp(Y);
    X <- exp(X);
}

if (x$Data$error=="prop"){
    plot(refx,xlim=c(min,max),ylim=c(min,max),xaxt="n",yaxt="n",log="xy",xlab=x$Data$xlab,ylab=x$Data$ylab);
    axis(1,at=round(c(min,x$Data$by,max),x$Data$dec));
    axis(2,at=round(c(min,x$Data$by,max),x$Data$dec));
}else{
    minn <- min-range/20;
    maxx <- max+range/20;
    plot(refx,xlim=c(minn,maxx),ylim=c(minn,maxx),xaxt="n",yaxt="n",xlab=x$Data$xlab,ylab=x$Data$ylab);
    axis(1,at=seq(min,max,x$Data$by));
    axis(2,at=seq(min,max,x$Data$by));
}
points(X,Y,cex=2,...);
}

