draw.xy <-
function(x,y,xx,yy,xlim=NULL,ylim=NULL,width=1,height=0.5,bg=NULL,border=1,type='p',col=1,silent=TRUE,...){
  n <- length(x)
  if(length(y)!=n | length(xx)!=n | length(yy)!=n) stop("Arguments 'x', 'y', 'xx' and 'yy' should all be numeric vectors of the same length")
  if(length(col)>1 & length(col)!=n) warning("Argument 'col' should be single colour or a vector of colours with the same length as 'x'")
  if(is.null(xlim)) xlim <- range(xx,na.rm=T)
  if(is.null(ylim)) ylim <- range(yy,na.rm=T)  
  marg <- 0.8
  xxs <- x+marg*width*(xx-xlim[1])/(xlim[2]-xlim[1])-0.5*width*marg
  yys <- y+marg*height*(yy-ylim[1])/(ylim[2]-ylim[1])-0.5*height*marg
  col <- rep(col,length.out=n)
  df <- data.frame(x,y,xxs,yys,col)
  xy <- unique(paste(x,y))
  if(!silent)message('Converting data to list...')
  lst <- lapply(xy,function(xy) subset(df,paste(x,y)==xy))
  if(!silent) message('...data converted.')
  if(!silent) message('Plotting data:')  
  nxy <- length(lst)
  pm <- setProgressMsg(1,nxy)
  for(i in 1:nxy){
    lsti <- lst[[i]]
    xi <- lsti$x[1]
    yi <- lsti$y[1]
    xxi <-lsti$xxs
    yyi <- lsti$yys
    yy0 <- yi-0.5*height*marg
    coli <- as.character(lsti$col)
    rect(xi-0.5*width,yi-0.5*height,xi+0.5*width,yi+0.5*height,col=bg,border=border)
    if(type=='h') segments(xxi,yyi,xxi,yy0,lend=2,coli,...) else points(xxi,yyi,type,col=coli,...)
    if(!silent) pm <- progressMsg(pm,i)
    }
}

