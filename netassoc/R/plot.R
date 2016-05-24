
# plot network
plot_netassoc_network <- function(network, layout = layout.fruchterman.reingold(network), 
                                  vertex.label = V(network)$name, 
                                  vertex.color = NA, 
                                  vertex.shape = "none",
                                  vertex.label.color = "black", 
                                  vertex.label.family = "sans",
                                  edge.width = NULL, 
                                  edge.color = NULL, 
                                  edge.arrow.size = 0.2, 
                                  vertex.label.cex = 0.5, 
                                  legend=TRUE,
                                  ...)
{    
  if(is.null(edge.width))
  {
    if(length(E(network)$weight)==0)
    {
      edge.width=1
    }
    else
    {
      edge.width=abs(E(network)$weight)
    }
  }
  
  if(is.null(edge.color))
  {
    if(length(E(network)$weight)==0)
    {
      edge.color <- 'black'  
      zlmin <- -1
      zlmax <- 1
    }
    else
    {
      edge.color <- ifelse(E(network)$weight > 0, rgb(0,0,1),rgb(1,0,0))
      
      if (length(E(network)$weight)>0)
      {
        zlmax <- max(abs(E(network)$weight),na.rm=T)
        zlmin = -1*zlmax
      }
      else
      {
        zlmin <- -1
        zlmax <- 1
      }
    }
  }
  
  plot(network,
       layout=layout,
       vertex.label=vertex.label,
       edge.color=edge.color,
       edge.width=edge.width,
       vertex.color=vertex.color,
       vertex.label.color=vertex.label.color,
       vertex.shape=vertex.shape,
       edge.arrow.size=edge.arrow.size,
       vertex.label.cex=vertex.label.cex,
       vertex.label.family=vertex.label.family,
       ...)
  
  colors=colorRampPalette(c('red','white','blue'))(51)
  
  if (legend==TRUE)
  {
    legend('topleft',adj=c(0,0),legend=format(c(zlmin,zlmin/2+zlmax/2,zlmax),digits=2),fill=c(colors[1],colors[ceiling(length(colors)/2)],colors[length(colors)]),bg='white',cex=0.5)
  }
}




plot_netassoc_matrix <- function(data, colors, onesided=FALSE, main="", legend=TRUE, axis=TRUE, title=TRUE, cex.axis=0.5)
{
  
  if (length(na.omit(as.numeric(data))) > 0)
  {
    zlmax <- max(abs(as.numeric(data)),na.rm=T)
    if (is.infinite(zlmax))
    {
      zlmax <- 1
    }
  }
  else
  {
    zlmax <- 1
  }
  if (onesided==TRUE)
  {
    zlmin = 0
  }
  else
  {
    zlmin = -1*zlmax
  }
  
  image(t(data),col=colors,axes=F,zlim=c(zlmin, zlmax),main=ifelse(title==TRUE,main,""))
  
  if (axis==TRUE)
  {
    axis(side=1,dimnames(data)[[2]],at=seq(0,1,length.out=length(dimnames(data)[[2]])),cex.axis=cex.axis,las=2)
    axis(side=2,dimnames(data)[[1]],at=seq(0,1,length.out=length(dimnames(data)[[1]])),cex.axis=cex.axis,las=1)
  }
  
  if (legend==TRUE)
  {
    legend('topleft',adj=c(0,0),legend=format(c(zlmin,zlmin/2+zlmax/2,zlmax),digits=2),fill=c(colors[1],colors[ceiling(length(colors)/2)],colors[length(colors)]),bg='white',cex=0.6)
  }
  
  box()
}


