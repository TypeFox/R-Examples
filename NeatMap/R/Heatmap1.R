data.reduction<-function(profiles,method="nMDS",metric="pearson",random.seed=NULL)
{
  METHODS=c("nMDS","PCA","average.linkage","complete.linkage");
  meth<-pmatch(method,METHODS);
  if(is.na(meth)) stop("Invalid Data Reduction Method");
  if(meth == -1) stop("Ambiguous Data Reduction Method");
  METRICS<-c("pearson","euclidean");
  metr<-pmatch(metric,METRICS);
  if(is.na(metr)) stop("Invalid Metric");
  if(metr == -1) stop("Ambiguous Metric");
  switch(meth,if(metr==1) nMDS(profiles,random.seed=random.seed) else nMDS(profiles,metric="euclidean",random.seed=random.seed),prcomp(profiles),if(metr==1) hclust(as.dist(1-cor(t(profiles),use="pairwise.complete.obs")),method="average") else hclust(dist(profiles),method="average"),if(metr==1) hclust(as.dist(1-cor(t(profiles),use="pairwise.complete.obs")),method="complete") else hclust(dist(profiles),method="complete"))
}

find.order<-function(reduction.result)
{
        METHODS=c("nMDS","prcomp","hclust");
        meth<-pmatch(class(reduction.result),METHODS);
        if(is.na(meth)) stop("Invalid Data Reduction Method");
        if(meth == -1) stop("Ambiguous Data Reduction Method");
        switch(meth,order(RadialCoords(reduction.result$x)[,2]),order(RadialCoords(reduction.result$x[,1:2])[,2]),reduction.result$order);
}

make.heatmap1<-function(profiles,row.method="nMDS",column.method="none",row.metric="pearson",column.metric="pearson",row.cluster.method="average",column.cluster.method="average",column.labels=NULL,row.labels=NULL,row.label.size=3,column.label.size=3,row.normalize=F, row.random.seed=NULL,column.random.seed=NULL)
{
  
  if(data.class(profiles)!="matrix")
  {
      profiles<-as.matrix(profiles);
  }
  METHODS=c("none","nMDS","PCA","average.linkage","complete.linkage");
  rmeth<-pmatch(row.method,METHODS);
  if(is.na(rmeth)) stop("Invalid Row Method");
  if(rmeth == -1) stop("Ambiguous Row Method");  
  if(rmeth ==1) row.order<-1:dim(profiles)[1] else
  row.order<-find.order(data.reduction(profiles,method=row.method,metric=row.metric,random.seed=row.random.seed));

  cmeth<-pmatch(column.method,METHODS);
  if(is.na(cmeth)) stop("Invalid Column Method");
  if(cmeth == -1) stop("Ambiguous Column Method");  
  if(cmeth ==1) column.order<-1:dim(profiles)[2] else
  column.order<-find.order(data.reduction(t(profiles),method=column.method,metric=column.metric,random.seed=column.random.seed));

  CMETHODS=c("none","average.linkage","complete.linkage");
  r.clust<-pmatch(row.cluster.method,CMETHODS);
  if(is.na(r.clust)) stop("Invalid Row Cluster Method");
  if(r.clust == -1) stop("Ambiguous Row Cluster Method");
  if(r.clust==1) row.cluster = NULL else row.cluster<-data.reduction(profiles, method=row.cluster.method,metric=row.metric);

  c.clust<-pmatch(column.cluster.method,CMETHODS);
  if(is.na(c.clust)) stop("Invalid Column Cluster Method");
  if(c.clust == -1) stop("Ambiguous Column Cluster Method");
  if(c.clust==1) column.cluster=NULL else column.cluster<-data.reduction(t(profiles), method=column.cluster.method,metric=column.metric);
  
  heatmap1(profiles,row.order=row.order,column.order=column.order,row.cluster=row.cluster,column.cluster=column.cluster,row.labels=row.labels,column.labels=column.labels,row.label.size=row.label.size,column.label.size=column.label.size,row.normalize=row.normalize);
  
}


heatmap1<-function(profiles,row.order=NULL,column.order=NULL,row.cluster=NULL,column.cluster=NULL,column.labels=NULL,row.labels=NULL,column.label.size=3,row.label.size=3,row.normalize=F)
{
  profiles<-as.matrix(profiles);
  n.points=dim(profiles)[1];
  profile.length=dim(profiles)[2];
  row.ranks<-vector(length=n.points);
  column.ranks<-vector(length=profile.length);
  if(is.null(row.order))
  {
    row.order=1:n.points;
  }
  else
  {
       if(!setequal(row.order,1:n.points))
       {
        stop("row.order not a permutation of 1:(number_of_rows)");
       }
  }
  if(is.null(column.order))
  {
    column.order=1:profile.length;
  }
  else
  {
       if(!setequal(column.order,1:profile.length))
       {
        stop("column.order not a permutation of 1:(number_of_columns)");
       }
  }
  for(i in 1:n.points)
  {
    row.ranks[row.order[i]]=i;
  }
  for(i in 1:profile.length)
  {
    column.ranks[column.order[i]]=i;
  }
  
 profiles1<-profiles;
 if(row.normalize)
 {
        for(i in 1:n.points)
        {
                profiles1[i,]<-(profiles[i,]-mean(profiles[i,],na.rm=T))/sd(profiles[i,],na.rm=T);
        }
 }

  myplot<-ggplot(data.frame(expand.grid(y=row.ranks,x=column.ranks),values=as.vector(profiles1)));
  if(!is.null(row.cluster))
  {
    if(class(row.cluster)!="hclust") stop("row.cluster must be of class hclust");
    myplot<-myplot+draw.dendrogram(row.cluster,scale=profile.length/10,leaf.order=row.order);  
  }
  if(!is.null(column.cluster))
  {
    if(class(column.cluster)!="hclust") stop("column.cluster must be of class hclust");
    myplot<-myplot+draw.dendrogram(column.cluster,scale=n.points/10,leaf.order=column.order,dendro.dir="up",order.dir="right",origin=as.vector(c(0,n.points)));    
  }
 if(!is.null(column.labels))
 {
        column.labels<-as.vector(column.labels);
        if(length(column.labels)!=profile.length) stop("column.labels must have number_of_columns elements");
        myplot<-myplot+draw.labels(column.labels,column.order,origin=c(0,0),order.dir="right",angle=270,size=column.label.size);
 }
 if(!is.null(row.labels))
 {
        row.labels<-as.vector(row.labels);
        if(length(row.labels)!=n.points) stop("row.labels must have number_of_rows elements");
        myplot<-myplot+draw.labels(row.labels,row.order,origin=c(profile.length+1,0),size=row.label.size)  ;
 }
#fill.data=data.frame(expand.grid(row.ranks,column.ranks),values=as.vector(profiles1));
 #myplot+geom_tile(data=data.frame(expand.grid(y=row.ranks,x=column.ranks),values=as.vector(profiles1)),aes_string(x=x,y=y,fill=values))+scale_fill_gradient2(low="green",high="red",mid="black",midpoint=mean(profiles1,na.rm=T))
 myplot+geom_tile(aes_string(x="x",y="y",fill="values"))+scale_fill_gradient2(low="green",high="red",mid="black",midpoint=mean(profiles1,na.rm=T));
}


primary.direction<-function(direction)
{
  DIRECTIONS=c("up","left","down","right");
  dir.no<-pmatch(direction,DIRECTIONS);
  if(!is.na(dir.no))
  {
    dirvec<-round(as.vector(c(cos(pi*dir.no/2),sin(pi*dir.no/2))));
  }
  else stop("ambiguous direction");
  dirvec
}


draw.labels<-function(labels,label.order=NULL,order.dir="up",label.colors=NULL,origin=c(0,0),angle=0,size=3)
{
        n.points=length(labels);
        if(is.null(label.order))
        {
                label.order=1:n.points;
        }
        label.dir<-primary.direction(order.dir);
        label.pos<-matrix(nrow=n.points,ncol=2);
        label.reorder=label.order;
	for( i in 1:n.points)
	{
	  label.reorder[label.order[i]]=i;
	}
        for(i in 1:n.points)
        {
                label.pos[i,]<-origin+label.dir*label.reorder[i];        
        }
        label.data<-data.frame(x=label.pos[,1],y=label.pos[,2],names=labels);
        geom_text(data=label.data,aes_string(x="x",y="y",label="names"),hjust=0,angle=angle,size=size);
}

draw.dendrogram<-function(cluster,leaf.order=NULL,scale=10,dendro.dir="left",order.dir="up",origin=as.vector(c(0.5,0)),heights=NULL)
{
  dendro.vec<-primary.direction(dendro.dir);
  order.vec<-primary.direction(order.dir);
  if(class(cluster)!="hclust") stop("cluster must be of class hclust");
  merge<-cluster$merge;
  if(is.null(heights))
  {
    dendro.heights<-cluster$height/(max(cluster$height));
  }
  else
  {
    heights<-as.vector(heights);
    if(length(heights)<dim(merge)) stop("all merge heights not specified");
    dendro.heights<-heights;
  }
  n.points<-dim(merge)[1]+1;
  if(is.null(leaf.order))
  {
        leaf.order<-as.vector(cluster$order);
  }
  else
  {
        leaf.order<-as.vector(leaf.order);
        if(length(leaf.order)!=n.points) stop("leaf.order has incorrect dimensionality (or type)");
  }
  
  height<-0;
  leaf.pos<-matrix(nrow=n.points,ncol=2);
  node.pos<-matrix(nrow=n.points-1,ncol=2);
  dendro.points<-matrix(nrow=4*(n.points-1),ncol=2);
  dendro.group<-matrix(nrow=4*(n.points-1),ncol=1); 
  for(i in 1:n.points)
  {      
      leaf.pos[leaf.order[i],]<-origin+i*order.vec;
  }
  for(i in 1:(n.points-1))
  {
      
      height<-dendro.heights[i]; 
      
      if(cluster$merge[i,1]<0)
      {
	  pos1<-as.vector(leaf.pos[-merge[i,1],]);      
      }
      else
      {
	  pos1<-as.vector(node.pos[merge[i,1],]);
      }
      pos1a<-origin+drop(pos1%*%order.vec)*order.vec+scale*height*dendro.vec;
      if(cluster$merge[i,2]<0)
      {
	  pos2=as.vector(leaf.pos[-merge[i,2],]);      
      }
      else
      {
	  pos2=as.vector(node.pos[merge[i,2],]);
      }
      pos2a<-origin+drop(pos2%*%order.vec)*order.vec+scale*height*dendro.vec;


      dendro.group[(4*(i-1)+1:4)]<-i;
      dendro.points[4*(i-1)+1,]<-pos1;    
      dendro.points[4*(i-1)+2,]<-pos1a;    
      dendro.points[4*(i-1)+3,]<-pos2a;    
      dendro.points[4*(i-1)+4,]<-pos2;   
      
      node.pos[i,]<-(pos1a+pos2a)/2.0;
    }
    dendrogram<-data.frame(x=dendro.points[,1],y=dendro.points[,2],group=dendro.group);

    geom_path(data=dendrogram,aes_string(x="x",y="y",group="group"),size=0.1)

}

