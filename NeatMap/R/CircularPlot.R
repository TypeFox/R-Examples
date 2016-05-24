make.circularmap<-function(profiles,method="nMDS",column.method="none",cluster.method="average.linkage",metric="pearson",column.metric="pearson",Rin=10,Rout=30,thickness=3,label.names=NULL,Rlabel=32,label.size=1.5,normalize.profiles=T, row.random.seed=NULL,column.random.seed=NULL)
{
        
  METHODS=c("nMDS","PCA");
  meth<-pmatch(method,METHODS);
  if(is.na(meth)) stop("Invalid Method");
  if(meth == -1) stop("Ambiguous Method");
  positions<-data.reduction(profiles,method,metric,random.seed=row.random.seed)$x[,1:2];
  
  METHODS1=c("none","nMDS","PCA","average.linkage","complete.linkage");
  cmeth<-pmatch(column.method,METHODS1);
  if(is.na(cmeth)) stop("Invalid Column Method");
  if(cmeth == -1) stop("Ambiguous Column Method");  
  if(cmeth ==1) column.order<-1:dim(profiles)[2] else
  column.order<-find.order(data.reduction(t(profiles),method=column.method,metric=column.metric,random.seed=column.random.seed));

  
  CMETHODS=c("none","average.linkage","complete.linkage");
  r.clust<-pmatch(cluster.method,CMETHODS);
  if(is.na(r.clust)) stop("Invalid Row Cluster Method");
  if(r.clust == -1) stop("Ambiguous Row Cluster Method");
  if(r.clust==1) row.cluster = NULL else row.cluster<-data.reduction(profiles, method=cluster.method,metric=metric);

  circularmap(positions,profiles,row.cluster,column.order=column.order,Rin=Rin,Rout=Rout,thickness=thickness,label.names=label.names,Rlabel=Rlabel,label.size=label.size,normalize.profiles=normalize.profiles);
}

normalize<-function(data)
{
        norm<-matrix(nrow=dim(data)[1],ncol=dim(data)[2]);
        for(i in 1:(dim(norm)[1]))
        {
                norm[i,]<-(data[i,]-mean(data[i,],na.rm=T))/sd(data[i,],na.rm=T);
        }
        return(norm);
}

circularmap<-function(pos,profiles,column.order=NULL,cluster.result=NULL,cluster.heights=NULL,Rin=10,Rout=30,thickness=3,label.names=NULL,Rlabel=32,label.size=1.5,normalize.profiles=T)
{
    profiles<-as.matrix(profiles);
    pos<-as.matrix(pos);
    profiles.length<-dim(profiles)[2];
    n.points<-dim(profiles)[1];
    theta<-matrix(nrow=n.points,ncol=1);

    if(dim(pos)[1]!=n.points) stop("pos and profiles must have same number of rows");
    if(dim(pos)[2]>=2)
    {
      theta<-RadialCoords(pos[,1:2])[,2];  
    }
    else
    {
      if(dim(pos)[2]==1)
      {
	theta<-as.vector(pos);
      }
    }
    if(is.null(column.order))
    {
        column.order=1:profiles.length;
    }
    else
    {
        if(!setequal(column.order,1:profiles.length))
        {
          stop("column.order not a permutation of 1:(number_of_columns)");
        }
    }
    x<-c();
    y<-c();
    color<-c();
    group<-c();
    labelx<-matrix(nrow=n.points,ncol=1);
    labely<-matrix(nrow=n.points,ncol=1);
    labelname<-matrix(nrow=n.points,ncol=1);
    labelangle<-matrix(nrow=n.points,ncol=1);
    i=0;
    origin<-data.frame(x=0,y=0);
    myplot<-qplot(x,y,data=origin);
    if(normalize.profiles)
    {
        profiles1<-normalize(profiles);
    }
    else
    {
        profiles1<-profiles;
    }
 for(j in 0:((length(theta)-0.5)%/%1000))
    {
        x0<-c();
        y0<-c();
        color0<-c();
        group0<-c();
        if(j<(length(theta)%/%1000))
        {
                lim=1000;
        }
        else
        {
                lim=length(theta)%%1000;
        }
        for(k in 1:lim)
        {
              i<-i+1;
              xstart<-Rin*cos(theta[i]);
              xstop<-Rout*cos(theta[i]);
              ystart<-Rin*sin(theta[i]);
              ystop<-Rout*sin(theta[i]);
              xvals<-seq(xstart,xstop,(xstop-xstart)/(profiles.length))
              yvals<-seq(ystart,ystop,(ystop-ystart)/(profiles.length))
              colorvals<-c(profiles1[i,column.order],0);
              x0<-c(x0,xvals); 
              y0<-c(y0,yvals); 
              color0<-c(color0,colorvals); 
              group0<-c(group0,rep(i,profiles.length+1)); 
              labelx[i]<-Rlabel*cos(theta[i]);
              labely[i]<-Rlabel*sin(theta[i]);
              labelangle[i]<-180*theta[i]/pi;
        }
        x<-c(x,x0);
        y<-c(y,y0);
        color<-c(color,color0);
        group<-c(group,group0);
        data<-data.frame(x=x0,y=y0,color=color0,group=group0,row.names=1:length(group0))
        myplot<-myplot+geom_path(data=data,aes_string(x="x",y="y",color="color",group="group"),size=thickness)+scale_colour_gradient2(low="green",high="red",mid="black",midpoint=mean(profiles1,na.rm=T));
    }
    
    origin<-data.frame(x=0,y=0);
  

    if(!is.null(label.names))
    {
      label.names<-as.vector(label.names);
      if(length(label.names)!=n.points) stop("The number of labels is not equal to the number of rows");
      labels<-data.frame(x=labelx,y=labely,angle=labelangle,label=label.names);
      myplot<-myplot+geom_text(data=labels,aes_string(x="x",y="y",angle="angle",label="label"),size=label.size,hjust=0);
    }
    
    if(!is.null(cluster.result))
    {
      if(class(cluster.result)!="hclust") stop("cluster.result must be of type hclust");
      myplot<-myplot+CircularDendrogram(pos,cluster.result,heights=cluster.heights,Rout=Rin);
    }
    myplot

}

CircularDendrogram<-function(pos,cluster,Rout=10,Rin=0,heights=NULL)
{
  n.points=dim(pos)[1];
  merge<-cluster[["merge"]];
  m.pos<-matrix(nrow=n.points-1,ncol=2);
  dendro.points<-matrix(nrow=4*(n.points-1),ncol=2);
  dendro.group<-matrix(nrow=4*(n.points-1),ncol=1); 

  for(i in 1:(n.points-1))
  {
    if(!is.null(heights))
    {
      Rnow=Rout-(Rout-Rin)*heights[i]/max(heights[i]);
    }
    else
    {
     Rnow=Rout-((Rout-Rin)*(cluster$height[i])/(max(cluster$height))); 
    }    
    if(merge[i,1]<0)
    {
      pos1=as.vector(pos[-merge[i,1],]);      
      theta<-atan2(pos1[2],pos1[1]);
      pos1<-as.vector(c(Rout*cos(theta),Rout*sin(theta)));
    }
    else
    {
      pos1=as.vector(m.pos[merge[i,1],]);
    }
    
    theta<-atan2(pos1[2],pos1[1]);
    pos1a<-as.vector(c(Rnow*cos(theta),Rnow*sin(theta)));

    if(merge[i,2]<0)
    {
      pos2=as.vector(pos[-merge[i,2],]);      
      theta<-atan2(pos2[2],pos2[1]);
      pos2<-as.vector(c(Rout*cos(theta),Rout*sin(theta)));
    }
    else
    {
      pos2=as.vector(m.pos[merge[i,2],]);
    }
    theta<-atan2(pos2[2],pos2[1]);
    pos2a<-as.vector(c(Rnow*cos(theta),Rnow*sin(theta)));
    dendro.group[(4*(i-1)+1:4)]<-i;
    dendro.points[4*(i-1)+1,]<-pos1;    
    dendro.points[4*(i-1)+2,]<-pos1a;    
    dendro.points[4*(i-1)+3,]<-pos2a;    
    dendro.points[4*(i-1)+4,]<-pos2;    
    m.pos[i,1]<-(pos1a[1]+pos2a[1])/2.0;
    m.pos[i,2]<-(pos1a[2]+pos2a[2])/2.0;
    
  }
  dendrogram<-data.frame(x=dendro.points[,1],y=dendro.points[,2],group=dendro.group);
  geom_path(data=dendrogram,aes_string(group="group",x="x",y="y"),size=0.1)
}
