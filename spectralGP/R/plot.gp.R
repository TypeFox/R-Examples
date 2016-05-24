"plot.gp" <-
function(x,type="l",col=terrain.colors(32),...){
  # plot an image of the process
  # require('fields')  # no longer needed as I added image.plot code to spectralGP
  m1=x$gridsize[1]
  m2=x$gridsize[2]
  if(x$d==2){
    gr=getgrid(x) 
    pred=predict(x)
    if(min(pred)==max(pred)){      
      image_plot(gr$x1,gr$x2,pred,xlab=expression(x[1]),ylab=expression(x[2]),col=col,zlim=c(min(pred)-1,max(pred)+1),...)
    } else{
      image_plot(gr$x1,gr$x2,pred,xlab=expression(x[1]),ylab=expression(x[2]),col=col,...)
    }
  } else{
    plot(getgrid(x),predict(x),xlab="x",ylab="f",type,...)
  }
}
