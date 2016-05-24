"triplot" <-
function(x, y=NULL, z=NULL,
         labels=dimnames(x)[[2]],
         txt=dimnames(x)[[1]], legend=NULL,
         legend.split=NULL,
         inner=TRUE, inner.col=c('lightblue','pink'),
         inner.lty=c(2,3), add=FALSE, main="", ...){

  old.par <- par(xpd=TRUE)
  on.exit(par(old.par))

  if( is.data.frame(x) ) x <- as.matrix(x)
  x <- cbind(x,y,z)
  if( ncol(x) < 2 || ncol(x) > 3 ){
    stop("need 2 or 3 columns")
  }
  if( ncol(x)==3 ){
    x <- sweep(x,1,FUN="/",apply(x,1,sum))
  }
  if( ncol(x)==2 ){
    x <- cbind(x, 1-x[,1]-x[,2])
  }

  if(dev.cur()==1){
    dev.new()
    add <- FALSE
  }

  if( !add ){

    pin <- par("pin")
    xstar <- (pin[1]/pin[2]*sqrt(3)-2)/2

    plot( c(0,1,2,0), c(0,sqrt(3),0,0), type="l",
         lwd=3, xlim=c(-xstar,2+xstar),
         xlab="",ylab="",axes=FALSE, main=main)

    if(inner){
      lines( c(1,1.5,0.5,1), c(0,sqrt(3)/2,sqrt(3)/2,0),
            lwd=.5, col=inner.col[1], lty=inner.lty[1])
      lines( c(1.25, 1, .75, 1.25),
            c(sqrt(3)/4, sqrt(3)/2, sqrt(3)/4, sqrt(3)/4),
            lwd=0.25, col=inner.col[2],lty=inner.lty[2])
    }

    if(length(labels)==0){
      labels <- c("X","Y","Z")
    }

    ystar <- par("cxy")[2] * 1.1
    text( c(0,2,1), c(-ystar,-ystar,sqrt(3)+ystar),
         labels, cex=1.5 )
  }

  newy <- x[,3] * sqrt(3)
  newx <- 2-2*x[,1]-x[,3]

  if(length(txt)==length(newx)){
    text(newx,newy,txt,...)
  } else {
    points(newx,newy,...)
  }

  if(length(legend)==length(newx)){
    labpos <- function(y){
      strh <- par("cxy")[2]*1.15
      y2 <- sort(y)
      df <- y2[-1] - y2[-length(y2)]
      i <- 1
      while(any (df < strh)){
        y2[c(df < strh, FALSE)] <- y2[ c(df < strh,FALSE)] - strh/10
        y2[c(FALSE, df < strh)] <- y2[ c(FALSE,df < strh)] + strh/10
        if(min(y2)<0){y2 <- y2 - min(y2)}
        y2 <- sort(y2)
        df <- y2[-1] - y2[ -length(y2)]
        i <- i+1
        if(i>100){break}
      }
      y2
    }

    if(length(legend.split)==1){
      tmp.x <- quantile(newx, legend.split)
      y1 <- newy[newx <= tmp.x]
      y1 <- labpos(y1)[order(order(y1))]
      text(rep(-0.01,length(y1)), y1, legend[newx<=tmp.x], adj=1)
      segments(rep(0,length(y1)), y1, newx[newx<=tmp.x], newy[newx<=tmp.x])

      y2 <- newy[newx>tmp.x]
      y2 <- labpos(y2)[order(order(y2))]
      text(rep(2.01,length(y2)), y2, legend[newx>tmp.x], adj=0)
      segments(rep(2,length(y2)), y2, newx[newx>tmp.x], newy[newx>tmp.x])
    } else {
      if(any(newx <= 1)){
        y1 <- newy[newx<=1]
        y1 <- labpos(y1)[order(order(y1))]
        text(rep(-0.01,length(y1)), y1, legend[newx<=1],adj=1)
        segments(rep(0,length(y1)), y1, newx[newx<=1], newy[newx<=1])
      }
      if(any(newx > 1)){
        y2 <- newy[newx>1]
        y2 <- labpos(y2)[order(order(y2))]
        text(rep(2.01,length(y2)), y2, legend[newx>1],adj=0)
        segments(rep(2,length(y2)), y2, newx[newx>1], newy[newx>1])
      }
    }
  }

  invisible(cbind(x=newx,y=newy))
}

