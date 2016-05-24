

#2-d plot with high regret colored outputs:

#object is a dataset, with output dimensions thresholded
#xdim and ydim are the column numbers referring to the x and y axis
#outdim is colnum for output dimension
#lowcol gives the color to fill in points below threshold
#hicol gives the color to fill in for points above threshold
#... is in case of passing on xlim and ylim

colptplot <- function(object,xdim,ydim,outdim=ncol(object),lowcol="transparent",
     hicol="black",xname=colnames(object)[xdim],yname=colnames(object)[ydim],xlims=NULL,ylims=NULL,...){
  
  x <- object[,xdim]
  y <- object[,ydim]
  
  collogi <- object[,outdim]==1         #logical vector of high/low
  colvect <- rep(lowcol, nrow(object))  #vector of all lowcolors
  
  colvect[collogi] <- hicol             #replace high values with hicol
  
  
  
  plot(x,y,col=hicol,pch=(collogi*15+1), xlab=xname, ylab=yname,xlim=xlims,ylim=ylims,...)
  #This format has all points the color of high color, but only those that are 
  #equal to 1 are actually filled in (controlled with pch - 1 or 16
  
}
  
