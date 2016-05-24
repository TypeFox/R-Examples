boxqq <-function(formula,...){
    UseMethod("boxqq")
}

boxqq.formula<-function(formula, data = NULL, ...)
### Box plots and Normal Quantile Quantile plots
{
  par("mar") -> mar

  if (missing(formula) || (length(formula) != 3))
      stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if (is.matrix(eval(m$data, parent.frame())))
      m$data <- as.data.frame(data)

  m[[1]] <- as.name("model.frame")
  mf <- eval(m, parent.frame())
  data<-mf

  par(mfrow=c(1,2))
  par(mar=c(5.1,4.1,4.1,0))
  boxplot(formula, data = data, xlab=names(data)[2],ylab=names(data)[1],main=paste(names(data)[1],"vs.",
                                                      names(data)[2]),cex=.75,col="yellow",...)
  names <- unique(data[,2])
  nv <- length(names)
  y <- data[,1]
  x <- NULL
  y <- NULL
  for (i in 1:nv) {
    foo <- qqnorm(data[data[,2]==names[i],1],plot.it=FALSE)
    x <- c(x,foo$x)
    y <- c(y,foo$y)
  }
  xl <- min(x)
  xm <- max(x)
  yl <- min(y)
  ym <- max(y)
  foo <- qqnorm(data[data[,2]==names[1],1],plot.it=FALSE)
  par(mar=c(5.1,0,4.1,2.1))
  plot(sort(foo$x),sort(foo$y),pch=1,xlim=c(xl,xm),ylim=c(yl,ym),xlab="Normal Quantiles",yaxt="n",type="b"
       ,main="Normal QQ Plots")
  leg.txt <- as.character(names)
  legend(x=-1.9,y=max(y)-.1*sd(y),legend=leg.txt,pch=c(1:nv),col=c(1:nv),bg="antiquewhite2")
  abline(mean(data[data[,2]==names[1],1]),sd(data[data[,2]==names[1],1]))
  for (j in 2:nv) {
    goo <- qqnorm(data[data[,2]==names[j],1],plot.it=FALSE)
    points(sort(goo$x),sort(goo$y),pch=j,col=j,type="b")
    abline(mean(data[data[,2]==names[j],1]),sd(data[data[,2]==names[j],1]),col=j)
  }
  par(mar=mar)
  par(mfrow=c(1,1))
}

