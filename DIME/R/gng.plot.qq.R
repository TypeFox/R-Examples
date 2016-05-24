gng.plot.qq <-
function (data, obj, resolution=10,xlab=NULL,ylab=NULL,main=NULL,pch=NULL,...)
{
x <- unlist(data);
if(is.null(resolution)) resolution <- 10;
    r <- runif(length(x)*resolution)
    y <- apply(matrix(r), 1, gng.qq.plot.internal, obj);
   len <- length(x);
   if(is.null(xlab)) xlab = "GNG";
   if(is.null(ylab)) ylab = "Observations";
   if(is.null(main)) main = "QQ-plot";
   if(is.null(pch)) pch = "*"; 
   qqplot(y, x, xlab = xlab, ylab = ylab, main = main, pch = pch,...);
    abline(0, 1);

}

