
showHeatmaps <- function(result,order_by=0,show_all=TRUE){
#result is an object of class JIVE_result
#order_by specifies how to order the rows and columns of the heatmap.  If order_by=0, orderings are determined by joint structure.  Otherwise, order_by gives the number of the individual structure dataset to determine the ordering.  In all cases orderings are determined by complete-linkage hiearchichal clustering of Euclidian distances.
#show_all specifies whether to show the full decomposition of the data, JIVE estimates, and noise.  If show_all=FALSE, only the matrix (or matrices) that determined the column ordering is shown.   

  old.par <- par(no.readonly = TRUE) # all par settings which could be changed
  on.exit(par(old.par))
  
  l <- length(result$data)
  ####Get row/column orderings
  Mat_ColOrder <- do.call(rbind,result$joint)
  if(order_by>0) Mat_ColOrder = result$individual[[order_by]]
  col.dist<-dist(t(Mat_ColOrder))
  rm(Mat_ColOrder)
  col.order<-hclust(col.dist)$order
  row.orders = list()
  for(i in 1:l){
    if(order_by==0) row.dist <- dist(result$joint[[i]])
    else row.dist <- dist(result$individual[[i]])
    row.orders[[i]] <- hclust(row.dist)$order
  }
  
  Image_Joint = list()
  for(i in 1:l){ Image_Joint[[i]] = as.matrix(result$joint[[i]][row.orders[[i]],col.order]) }
  
  Image_Data = list()
  for(i in 1:l){ Image_Data[[i]] = as.matrix(result$data[[i]][row.orders[[i]],col.order]) }
  
  Image_Indiv = list()
  for(i in 1:l){ Image_Indiv[[i]] = as.matrix(result$individual[[i]][row.orders[[i]],col.order]) }
  
  Image_Noise = list()
  for(i in 1:l){ Image_Noise[[i]] = Image_Data[[i]]-Image_Joint[[i]]-Image_Indiv[[i]] }
  
  par(mar=c(1,1,1,1))
  ### Make heatmap of all estimates###
  if(show_all){
    Heights=c(1,rep(3,l))
    Widths=c(1,5,1,5,1,5,1,5)
    layoutVals = c((8+7*l):(8+8*l),2:(2+l),(5+4*l):(5+5*l),(3+l):(3+2*l),(6+5*l):(6+6*l),1,(4+2*l):(3+3*l),(7+6*l):(7+7*l),(4+3*l):(4+4*l))
    layout(matrix(layoutVals,1+l,8), heights = Heights,widths = Widths)
    layout.show(8+8*l)
    CEX <- textplot("Individual") ###Plot this first, to get cex size
    textplot("Data", cex =CEX)
    for(i in 1:l) show.image(Image_Data[[i]],ylab=names(result$data)[i])  
    textplot("Joint",cex = CEX)
    for(i in 1:l) show.image(Image_Joint[[i]])
    #textplot("Individual") 
    for(i in 1:l) show.image(Image_Indiv[[i]])
    textplot("Noise", cex = CEX)
    for(i in 1:l) show.image(Image_Noise[[i]])
    textplot(' ')
    for(i in 1:l) textplot('=')
    textplot(' ')
    for(i in 1:l) textplot('+')
    textplot(' ')
    for(i in 1:l) textplot('+')
    textplot(' ')
    par(mar=c(0,0,0,0))
    for(i in 1:l){ plot.new();text(1/2,1/2,names(result$data)[[i]],srt=90,cex=CEX*0.75)}
    }
  else{ 
    if(order_by==0){
      layout(matrix(c((l+2):(2*l+2),1:(l+1)),l+1,2),heights=c(1,rep(3,l)),widths=c(1,5))
      layout.show(2+2*l)
      CEX=textplot("Joint")
      for(i in 1:l) show.image(Image_Joint[[i]])
      textplot(' ')
      par(mar=c(0,0,0,0))
      for(i in 1:l){ plot.new();text(1/2,1/2,names(result$data)[[i]],srt=90,cex=CEX*0.5)}
    }
    if(order_by>0){
      layout(matrix(c(3,4,1,2),2,2),heights=c(1,5),widths=c(1,5))
      layout.show(2+l)
      CEX=textplot("Individual")
      show.image(Image_Indiv[[order_by]])
      textplot(' ')
      plot.new();text(1/2,1/2,names(result$data)[[order_by]],srt=90,cex=CEX*0.75)
    }
}
}


show.image = function(Image,ylab=''){
  lower = mean(Image)-3*sd(Image)
  upper = mean(Image)+3*sd(Image)
  Image[Image<lower] = lower
  Image[Image>upper] = upper
  image(x=1:dim(Image)[2], y=1:dim(Image)[1], z=t(Image), zlim = c(lower,upper),axes=FALSE,col=bluered(100),xlab="",ylab=ylab)
} 



