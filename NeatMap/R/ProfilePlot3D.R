RadialCoords<-function(pos)
{
  pos<-as.matrix(pos);
  n.points=dim(pos)[1];
  radial<-matrix(nrow=n.points,ncol=2);
  xc<-mean(pos[,1]);
  yc<-mean(pos[,2]);
  for(i in 1:n.points)
  {
    radial[i,1]<-sqrt((pos[i,2]-yc)^2+(pos[i,1]-xc)^2);
    radial[i,2]<-atan2(pos[i,2]-yc,pos[i,1]-xc);    
  }
  rownames(radial)<-rownames(pos);
  colnames(radial)<-c("r","theta");
  radial
}

make.profileplot3d<-function(profiles,row.method="nMDS", normalize.rows=T, column.method="average.linkage",row.metric="pearson",column.metric="pearson",row.cluster.method="average",column.cluster.method="average", point.size=3,col=NULL,color_scaling_function=NULL, labels=NULL, label.colors=NULL,label.size=0.5, row.random.seed=NULL,column.random.seed=NULL)
{
  profiles<-as.matrix(profiles);
  ROW.METHODS=c("nMDS","PCA");
  meth<-pmatch(row.method,ROW.METHODS);
  if(is.na(meth)) stop("Invalid Row Method");
  if(meth == -1) stop("Ambiguous Row Method");
  positions<-data.reduction(profiles,row.method,row.metric,row.random.seed)$x[,1:2];

  
  COLUMN.METHODS=c("none","nMDS","PCA","average.linkage","complete.linkage");
  cmeth<-pmatch(column.method,COLUMN.METHODS);
  if(is.na(cmeth)) stop("Invalid Column Method");
  if(cmeth == -1) stop("Ambiguous Column Method");  
  if(cmeth ==1) column.order<-1:dim(profiles)[2] else column.order<-find.order(data.reduction(t(profiles),method=column.method,metric=column.metric,random.seed=column.random.seed));
 

  CMETHODS=c("none","average.linkage","complete.linkage");
  r.clust<-pmatch(row.cluster.method,CMETHODS);
  if(is.na(r.clust)) stop("Invalid Row Cluster Method");
  if(r.clust == -1) stop("Ambiguous Row Cluster Method");
  if(r.clust==1) row.cluster = NULL else row.cluster<-data.reduction(profiles, method=row.cluster.method,metric=row.metric);

  c.clust<-pmatch(column.cluster.method,CMETHODS);
  if(is.na(c.clust)) stop("Invalid Column Cluster Method");
  if(c.clust == -1) stop("Ambiguous Column Cluster Method");
  if(c.clust==1) column.cluster=NULL else column.cluster<-data.reduction(t(profiles), method=column.cluster.method,metric=column.metric);

  profileplot3d(positions,profiles,column.order=column.order,normalize.rows=normalize.rows,row.cluster=row.cluster,column.cluster=column.cluster,point.size=point.size,labels=labels,label.colors=label.colors,label.size=label.size,col=col,color_scaling_function=color_scaling_function)
}


DynamicLabels<-function(button=3, id0, positions,labels,colors=NULL)
{
        start<-list();
        starty<-0.5;
        id<-id0;
        height<-NULL;

        begin<-function(x,y)
        {
                start$viewport<-par3d("viewport");
                width<-start$viewport[3];
                height<-start$viewport[4];
                starty<-y;
        }
        update <- function(x,y) 
        {
                lambda<-0.5*(starty-y)/height;
                mysize <- clamp(par3d("cex")+lambda, 0.1, 10)
                rgl.pop(id=id)
                par3d("cex"=mysize)
                if(is.null(colors))
                {
                        id<-text3d(positions,texts=labels,cex=mysize,adj=0)
                }
                else
                {
                        id<-text3d(positions,texts=labels,cex=mysize,color=colors,adj=0)
                }
        }
        rgl.setMouseCallbacks(button, begin, update=update,end=NULL)
}


vlen <- function(a) sqrt(sum(a^2))
clamp <- function(x, min, max)  min(max(x, min), max)
profileplot3d<-function(pos,profiles,normalize.rows=T,column.order=NULL,row.cluster=NULL,column.cluster=NULL,labels=NULL,col=NULL,color_scaling_function=NULL,point.size=3,label.colors=NULL,label.size=0.5)
{
    
    if(is.null(col))
    {
      ColorFunction<-function(x)
      {
	if(is.na(x))
	{
	    rgb(150,150,150,maxColorValue=255)  
	}
	else
	{
	  rgb(colorRamp(c("Red","Black","Green"),space="rgb")(x),maxColorValue=255)      
	}
      }
    }
    else
    {
      ColorFunction<-function(x)
      {
      if(is.na(x))
	{
	    rgb(150,150,150,maxColorValue=255)  
	}
	else
	{
	  rgb(colorRamp(col,space="rgb")(x),maxColorValue=255)      
	}
      }
	
    }
   
    profiles<-as.matrix(profiles);
    pos<-as.matrix(pos);
    pos<-pos[,1:2];
    n.points<-dim(profiles)[1];
    profile.length<-dim(profiles)[2];
    if(dim(pos)[1]!=n.points) stop("pos and profiles must have the same number of rows");
    profiles1<-matrix(nrow=n.points,ncol=profile.length);
    
    scale<-(max(pos[,1],na.rm=T)+max(pos[,2],na.rm=T)-(min(pos[,1],na.rm=T)+min(pos[,2],na.rm=T)))/2.0; 
    if(is.null(column.order)) 
    {
        column.order<-1:profile.length;
    }
    else
    {
       if(!setequal(column.order,1:profile.length))
       {
        stop("column.order not a permutation of 1:(number_of_columns)");
       }
    }
    if(!is.null(column.cluster))
    {
      if(class(column.cluster)!="hclust") stop("column.cluster must be of type hclust");
      radial<-RadialCoords(pos);	  
      R<-max(radial[,1]);
      theta<-radial[which.max(radial[,1]),2];
      dirvec<-as.vector(c(cos(theta),sin(theta),0));
      column.dendrogram.pos<-matrix(nrow=profile.length,ncol=3);
      xc=mean(pos[,1]);
      yc=mean(pos[,2]);
      for(i in 1:profile.length)
      {
	column.dendrogram.pos[i,1]=xc+R*cos(theta);
	column.dendrogram.pos[i,2]=yc+R*sin(theta);
	column.dendrogram.pos[column.order[i],3]=scale*(i-1)/profile.length;
      }
      draw.dendrogram3d(column.cluster,column.dendrogram.pos,dirvec,scale/5)
      
    }
    if(!is.null(row.cluster))
    {
      if(class(row.cluster)!="hclust") stop("row.cluster must be of type hclust");
      draw.dendrogram3d(row.cluster,cbind(pos,rep(0,n.points)),c(0,0,-1),scale/3)
    }
    if(normalize.rows)
    {
      for(i in 1:n.points)
      {
	profiles1[i,]<-(profiles[i,]-mean(profiles[i,],na.rm=T))/sd(profiles[i,],na.rm=T);
      }
     
    }
     else
      {
	profiles1=profiles;
      }
    profiles.max<-max(profiles1,na.rm=T);
    profiles.min<-min(profiles1,na.rm=T);
    positions<-matrix(nrow=n.points*profile.length,ncol=3);
    colors<-matrix(nrow=n.points*profile.length,ncol=1);
    for(i in 1:(n.points))
    {
        for(j in 1:(profile.length))
        {
            positions[(i-1)*profile.length+j,1]=pos[i,1];
            positions[(i-1)*profile.length+j,2]=pos[i,2];
            positions[(i-1)*profile.length+j,3]=scale*(j-1)/profile.length;
            if(is.na(profiles[i,column.order[j]]))
            {
              colors[(i-1)*profile.length+j]=gray(0.5);
            }
            else
            {
              if(is.null(color_scaling_function))
              {
	        colors[(i-1)*profile.length+j]=ColorFunction((profiles1[i,column.order[j]]-profiles.min)/(profiles.max-profiles.min));
              }
              else
              {
                colors[(i-1)*profile.length+j]=ColorFunction(color_scaling_function((profiles1[i,column.order[j]]-profiles.min)/(profiles.max-profiles.min)));
              }
            }
        }
    }
    
    if(!is.null(labels))
    {
      labels<-as.vector(labels);
      if(length(labels)!=n.points) stop("labels has incorrect dimensions");
      textpos<-matrix(nrow=n.points,ncol=3);
      for(i in 1:n.points)
      {
	textpos[i,1]<-pos[i,1];
	textpos[i,2]<-pos[i,2];
	textpos[i,3]<-1.05*scale;
      }
      if(is.null(label.colors))
      {
	id<-text3d(textpos,texts=labels,cex=label.size,adj=0);
        DynamicLabels(id0=id,positions=textpos,labels=labels);
      }
      else
      {
        label.colors<-as.vector(label.colors);
        if(length(label.colors)!=n.points) stop("label.colors has incorrect dimensions");
	id<-text3d(textpos,texts=labels,color=label.colors,cex=label.size,adj=0);
        DynamicLabels(id0=id,positions=textpos,labels=labels,colors=label.colors);
      }
    }
    
    
    points3d(positions,size=point.size,col=colors);
}
