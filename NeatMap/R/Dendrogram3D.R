draw.dendrogram3d<-function(cluster,positions,direction=c(0,0,-1),scale=NULL,heights=NULL,labels=NULL,label.colors=NULL,label.size=3)
{
  if(class(cluster)!="hclust")
  {
    stop("cluster must be of type hclust");
  }
  merge<-cluster$merge;
  n.points<-dim(merge)[1]+1;
  positions<-as.matrix(positions);
  if(dim(positions)[1]!=n.points)  
  {
    stop("dimensionalities of cluster and positions do not match");
  }
  if(dim(positions)[2]==2)
  {
    positions<-cbind(positions,rep(0,n.points));
  }
  if(is.null(heights))
  {
    heights<-cluster$height/(max(cluster$height));
  }
  else
  {
        heights<-as.vector(heights);
        if(length(heights)!=n.points-1) stop("heights has incorrect dimensionality");
  }
 if(is.null(scale))
  {
    scale=((max(positions[,1])+max(positions[,2])+max(positions[,3]))-(min(positions[,1])+min(positions[,2])+min(positions[,3])))/6;
  }
 
  
  
  leaf.pos<-matrix(nrow=n.points,ncol=3);
  node.pos<-matrix(nrow=n.points-1,ncol=3);
  dendro.points<-matrix(nrow=6*(n.points-1),ncol=3);  
  for(i in 1:n.points)
  {      
      leaf.pos[i,]<-positions[i,];
  }
  
  for(i in 1:(n.points-1))
  {
      
      if(cluster$merge[i,1]<0)
      {
	  pos1<-as.vector(leaf.pos[-merge[i,1],]);      
	  height.old<-0;
      }
      else
      {
	  pos1<-as.vector(node.pos[merge[i,1],]);
	  height.old<-heights[merge[i,1]];
      }
      pos1a<-pos1+scale*(heights[i]-height.old)*direction;
      if(cluster$merge[i,2]<0)
      {
	  pos2=as.vector(leaf.pos[-merge[i,2],]);      
	  height.old<-0;
      }
      else
      {
	  pos2=as.vector(node.pos[merge[i,2],]);
	  height.old<-heights[merge[i,2]];
      }
      pos2a<-pos2+scale*(heights[i]-height.old)*direction;
      
      dendro.points[6*(i-1)+1,]<-pos1;    
      dendro.points[6*(i-1)+2,]<-pos1a;    
      dendro.points[6*(i-1)+3,]<-pos1a;    
      dendro.points[6*(i-1)+4,]<-pos2a;    
      dendro.points[6*(i-1)+5,]<-pos2a;    
      dendro.points[6*(i-1)+6,]<-pos2;    

      node.pos[i,]<-(pos1a+pos2a)/2.0;
    }
    if(!is.null(labels))
    {
       labels<-as.vector(labels);
       if(length(labels)!=n.points) stop("incorrect number of labels");
      if(is.null(label.colors))
      {
	text3d(positions,texts=labels,cex=label.size);
      }
      else
      {
        labels.colors<-as.vector(labels.colors);
         if(length(label.colors)!=n.points) stop("incorrect number of label.colors");
	text3d(positions,texts=labels,color=label.colors,cex=label.size);
      }
    }
    segments3d(dendro.points,color='black');  
    
}
