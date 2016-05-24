draw.grid <-
function(grd,breaks=NULL,col=NULL){
  
  # use default settintgs if no breaks are provided
  if(missing(breaks)) breaks <- breaks.grid(grd)
  ncol <- length(breaks)-1
#  if(missing(col)) col <- c('white',heat.colors(ncol-1)[(ncol-1):1])
  if(missing(col)) col <- colorRampPalette(c("lightyellow","yellow","orange","red", "brown4"))(ncol)
  if(length(col)!=ncol) 
    stop('The number of breakpoints should be one more than the number of colours')
  
  # replace values outside breakpoins so they are not blank
  grd <- ifelse(grd>max(breaks),max(breaks),grd)
  grd <- ifelse(grd<min(breaks),min(breaks),grd)
  
  # add heatmap to existing plot
  x <- as.numeric(rownames(grd))
  y <- as.numeric(colnames(grd))
  image(x,y,grd,col=col,breaks=breaks,add=T)
  box()
  }

