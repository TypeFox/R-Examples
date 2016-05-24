sym.scatterplot.ggplot <-
function(sym.var.x,sym.var.y,labels=FALSE,...) {
  x<-0
  y<-0
  xmin<-0
  xmax<-0
  ymin<-0
  ymax<-0
  if(((sym.var.x$var.type!='$C')||(sym.var.y$var.type!='$C'))&&
       ((sym.var.x$var.type!='$I')||(sym.var.y$var.type!='$I')))
    stop("Impossible to plot this type of variable")
  if((sym.var.x$var.type=='$C')&&(sym.var.y$var.type=='$C')) {
    df <- data.frame(sym.var.x$var.data.vector,sym.var.y$var.data.vector)
    names(df) <- c("x","y")
    p <- ggplot(df,aes(x,y)) + labs(x=sym.var.x$var.name,y=sym.var.y$var.name)
    if(labels==FALSE){
      p <- p + geom_point()
    }
    else {
      ltext<-sym.var.x$obj.names
      p <- p + geom_text(label=ltext)
    }
  }
  if((sym.var.x$var.type=='$I')&&(sym.var.y$var.type=='$I')) {    
    df <- data.frame(sym.var.x$var.data.vector,sym.var.y$var.data.vector)
    names(df) <- c("xmin","xmax","ymin","ymax")
    p <- ggplot(df) + labs(x=sym.var.x$var.name,y=sym.var.y$var.name) +
      geom_rect(data=df,aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),fill=alpha(1:sym.var.x$N,2/3))
    
    if(labels==TRUE) {
      ltext<-sym.var.x$obj.names
      p <- p + geom_text(aes((xmin+xmax)/2,(ymin+ymax)/2),label=ltext)
    }
  }
  print(p)
  invisible(p)
}
