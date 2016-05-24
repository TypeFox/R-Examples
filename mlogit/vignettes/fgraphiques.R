## graphique illustrant les séquences de Halton
ghalton <- function(n){
  require(mlogit)
  nf <- layout(matrix(c(1,2,3),3,1))
  layout.show(nf)
  opar <- par(bg="lightyellow")

  nb <- halton(3,27,0)
  namnb <- c("0","1/3","2/3",
             "0+1/9","1/3+1/9","2/3+1/9",
             "0+2/9","1/3+2/9","2/3+2/9",
             "0+1/27","1/3+1/27","2/3+1/27",
             "0+1/9+1/27","1/3+1/9+1/27","2/3+1/9+1/27",
             "0+2/9+1/27","1/3+2/9+1/27","2/3+2/9+1/27",
             "0+2/27","1/3+2/27","2/3+2/27",
             "0+1/9+2/27","1/3+1/9+2/27","2/3+1/9+2/27",
             "0+2/9+2/27","1/3+2/9+2/27","2/3+2/9+2/27")
  
  nbtext <- function(ind,dist,col){
    for (i in ind){
      lines(c(nb[i],nb[i]),c(-.1,.1),col=col,lwd=2)
      text(nb[i],-dist,namnb[i],cex=1)
      text(nb[i],.2,i,cex=1)
    }
  }
  
  yrange <- c(-.6,.6)
  opar <- par(mar=c(0,0,0,0))
  
  plot(0,xlim=c(-.2,1.2),ylim=yrange,type="n",lwd=2,axes=F,ann=F)
  rect(-0.2,-.6,1.2,.6,col="lightblue",lty=0)
  lines(c(0,1),c(0,0),lwd=2)
  nbtext(1:3,.2,"red")
  
  if (n>1){
    opar <- par(mar=c(0,0,0,0))
    plot(0,xlim=c(-.2,1.2),ylim=yrange,type="n",lwd=2,axes=F,ann=F)
    rect(-0.2,-.6,1.2,.6,col="lightblue",lty=0)
    lines(c(0,1),c(0,0),lwd=2)
    nbtext(1:3,.2,"red")
    nbtext(4:6,.3,"blue")
    nbtext(7:9,.3,"blue")
  }
  if (n>2){
    opar <- par(mar=c(0,0,0,0))
    plot(0,xlim=c(-.2,1.2),ylim=yrange,type="n",lwd=2,axes=F,ann=F)
    rect(-0.2,-.6,1.2,.6,col="lightblue",lty=0)
    lines(c(0,1),c(0,0),lwd=2)
    nbtext(1:3,.2,"red")
    nbtext(4:6,.3,"blue")
    nbtext(7:9,.3,"blue")
    nbtext(10:12,.4,"green")
    nbtext(13:15,.4,"green")
    nbtext(16:18,.4,"green")
  }

}

## Graphique illustrant la méthode de Newton-Ralphson


## nf <- layout(matrix(c(1,2,3, 4),2,2))
## layout.show(nf)
## goptimisation(1)
## goptimisation(2)
## goptimisation(3)
## goptimisation(4)



goptimisation <- function(n){
  a <- 1
  b <- 2
  f <- function(x) a*log(x)-b*x^2
  fp <- function(x) a/x-2*b*x
  fs <- function(x) -a/x^2-2*b
  fa <- function(x,xb) f(xb)+fp(xb)*(x-xb)+0.5*fs(xb)*(x-xb)^2
  xinf <- .1
  xsup <- 2.1
  xlim <- c(xinf,xsup)
  yinf <- -5
  ysup <- -0.5
  ylim <- c(yinf,ysup)
  maxa <- function(xb) xb-fp(xb)/fs(xb)
  xi <- 2
  opar <- par(bg="lightyellow", mar=rep(1, 4))
  curve(f,ann=F,axes=F,type="n",xaxs="i",yaxs="i",ylim=c(yinf,ysup),xlim=c(xinf,xsup), add=FALSE)
  rect(xinf,yinf,xsup,ysup,lty=0,col="lightblue")
  curve(f,add=T)
  #axis(side=1,tic=F,labels=F)
  #axis(side=2,c(seq(0,0.4,.1)),seq(0,0.4,.1),las=1)
  #abline(v=0)
  #title(main="Numerical maximisation",xlab="",
  #ylab="",cex.sub=.5,col.lab="blue")
  par(opar)
  for (i in 1:n){
    curve(fa(x,xb=xi),add=T,lty="dotted",col=palette()[i+1],lwd=2)
    xio <- xi
    xi <- maxa(xio)
    segments(xi,yinf,xi,max(f(xi),fa(xi,xio)),lty="dashed")
    segments(xinf,f(xi),xi,f(xi),lty="dashed")
    text(xi,yinf-.1,paste("x",i,sep=""))
  }
}

##  Graphique illustrant la génération de nombre de Gumbell
gunifgumbel <- function(n = 20){
  qgumbel <- function(x) -log(-log(x))
  qfunc <- function(x) qnorm(x)
  qfunc <- function(x) qgumbel(x)
  def.par <- par(no.readonly = TRUE , bg="lightyellow")
  x <- runif(n)
  x <- halton(3,n,10)
  y <- qfunc(x)
  xrange <- c(0,1)
  yrange <- c(-3,3.5)
  #yrange <- c(-2,2)
  xhist <- hist(x, breaks=seq(xrange[1],xrange[2],0.2), plot=FALSE)
  yhist <- hist(y[y>yrange[1] & y<yrange[2]],
                breaks=seq(yrange[1],yrange[2],.5), plot=FALSE)
  top <- max(c(xhist$counts, yhist$counts))
  nf <- layout(matrix(c(2,0,1,3),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)
  layout.show(nf)
  par(mar=c(3,3,1,1))
  #curve(qnorm,xlim=c(0,1),ylim=c(-3,3),las=1,xaxs="i",yaxs="i")
  curve(qfunc,xlim=xrange,ylim=yrange,las=1,xaxs="i",yaxs="i",type="n")
  rect(xrange[1],yrange[1],xrange[2],yrange[2],col="lightblue",lty=1)
  curve(qfunc,add=T)
  points(x,rep(yrange[2],length(x)),pch=15)
  points(rep(1,length(x)),qfunc(x),pch=15)
  abline(h=0)
  ## for (i in x){
  ##   lines(c(i,i),c(yrange[2],qfunc(i)),lty="dotted")
  ##   lines(c(i,1),c(qfunc(i),qfunc(i)),lty="dotted")
  ## }
  ##plot(x, y, xlim=xrange, ylim=yrange, xlab="", ylab="")
  par(mar=c(0,3,1,1))
  barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
  par(mar=c(3,0,1,1))
  barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)
  par(def.par)#- reset to default
}

## Graphique illustrant la corrélation entre variables aléatoires
gcorrelation <- function(n = 200){
  opar <- par(bg="lightyellow",las=1,xaxs="i",yaxs="i")
  nf <- layout(matrix(c(1,2),1,2))
  layout.show(nf)
  opar <- par(mar=c(3,3,3,3))
  x <- rnorm(n)
  y <- rnorm(n)
  v <- matrix(c(.5,.8,.8,2),2,2)
  xrange <- c(-3,3)*sqrt(v[1,1])
  yrange <- c(-3,3)*sqrt(v[2,2])
  ch <- chol(v)
  tr <- cbind(x,y)%*%ch
  x2 <- tr[,1]
  y2 <- tr[,2]
  plot(y~x,xlim=xrange,ylim=yrange,type="n")
  rect(xrange[1],yrange[1],xrange[2],yrange[2],col="lightblue",lty=1)
  points(x,y)
  plot(y2~x2,xlim=xrange,ylim=yrange,type="n")
  rect(xrange[1],yrange[1],xrange[2],yrange[2],col="lightblue",lty=1)
  points(x2,y2)
}

haltonvsrandom <- function(n = 100){
  opar <- par(bg="lightyellow",las=1,xaxs="i",yaxs="i")
  nf <- layout(matrix(c(1,2),1,2))
  layout.show(nf)
  xr <- runif(n)
  yr <- runif(n)
  xh <- halton(2,n,100)
  yh <- halton(3,n,100)
  opar <- par(mar=c(3,3,3,3))
  plot(yr~xr,ann=F,xlim=c(0,1),ylim=c(0,1),type="n")
  rect(0,0,1,1,col="lightblue",lty=1)
  points(xr,yr)
  opar <- par(mar=c(3,3,3,3))
  plot(yh~xh,ann=F,xlim=c(0,1),ylim=c(0,1),type="n")
  rect(0,0,1,1,col="lightblue",lty=1)
  points(xh,yh)
}
