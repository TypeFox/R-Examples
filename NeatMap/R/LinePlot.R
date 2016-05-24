clip.values<-function(x,max,min)
{
        val<-x;
        if(!is.na(x))
        {
                if(x>max) val=max;
                if(x<min) val=min;
        }
        val
}


lineplot<-function(pos,profiles,n.div.x=10,n.div.y=10,normalize=F,ylim=NULL,clipped=F)
{
    

    if(data.class(pos)!="matrix")
    {
        pos<-as.matrix(pos);
    }
    if(data.class(profiles)!="matrix")
    {
        profiles<-as.matrix(profiles);
    }
    if(dim(profiles)[1] != dim(pos)[1]) stop("pos and profiles must have same number of rows");
    
    n.points<-dim(profiles)[1];
    length.profile<-dim(profiles)[2];

    #normalize profiles if needed
    profiles.norm<-matrix(nrow=n.points,ncol=length.profile);
    if(normalize)
    {
      for(i in 1:n.points)
      {
	  profiles.norm[i,]<-(profiles[i,]-mean(profiles[i,],na.rm=T))/sd(profiles[i,],na.rm=T);
      }
    }
    else
    {
      profiles.norm<-profiles;
    }

    if(!is.null(ylim))
    {
        ylim<-as.vector(ylim);
        if(length(ylim)!=2)
        {
                stop("ylim must be of length 2");
        }
        max.y.value=max(ylim,na.rm=T);
        min.y.value=min(ylim,na.rm=T);
    }
    else
    {
        max.y.value=max(profiles.norm,na.rm=T);
        min.y.value=min(profiles.norm,na.rm=T);
    }
    #the grid limits extend 0.05 time the range beyond the max and min points 
    xmax=max(pos[,1])+0.05*(max(pos[,1])-min(pos[,1]))
    xmin=min(pos[,1])-0.05*(max(pos[,1])-min(pos[,1]))
    ymax=max(pos[,2])+0.05*(max(pos[,2])-min(pos[,2]))
    ymin=min(pos[,2])-0.05*(max(pos[,2])-min(pos[,2]))
    size.div.x<-(xmax-xmin)/n.div.x;
    size.div.y<-(ymax-ymin)/n.div.y;
    div.x<-seq(xmin,xmax,size.div.x);
    div.y<-seq(ymin,ymax,size.div.y);

    #place points in appropriate grid box
    xf<-cut(pos[,1],div.x)
    yf<-cut(pos[,2],div.y)

    xstart<-div.x+(size.div.x*0.05);
    ystart<-matrix(nrow=n.div.y,ncol=1);
    for(i in 1:(n.div.y))
    {
        ystart[i]<-(div.y[i]+div.y[i+1])/2;
    }

    grid_data<-data.frame(x=pos[,1],y=pos[,2],xbox=xf,ybox=yf,which.x=findInterval(pos[,1],div.x),which.y=findInterval(pos[,2],div.y))

    xscale<-(0.9*size.div.x)/length.profile;
    x<-matrix(nrow=(length.profile)*n.points,ncol=1)
    y<-matrix(nrow=(length.profile)*n.points,ncol=1)
    group<-matrix(nrow=(length.profile)*n.points,ncol=1)
    for(i in 1:(n.points))
    {
        xvals<-seq(0:(length.profile-1))*xscale+xstart[grid_data[["which.x"]][i]];
        
        if(clipped)
        {
                vals<-as.vector(sapply(profiles.norm[i,],function(x){clip.values(x,max.y.value,min.y.value)}));
        }
        else
        {
                vals<-profiles.norm[i,];
        }
        yvals<-(0.9*size.div.y)*((vals-min.y.value)/(max.y.value-min.y.value)-0.5)+ystart[grid_data[["which.y"]][i]];
        
        
        x[((i-1)*(length.profile)+1):(i*length.profile),1]<-xvals;
        y[((i-1)*(length.profile)+1):(i*length.profile),1]<-yvals;
        group[((i-1)*(length.profile)+1):(i*length.profile),1]<-i;
    }
    linedata<-data.frame(x=x,y=y,group=group);
    myplot<-ggplot();
    myplot+geom_path(aes_string(group="group",x="x",y="y"),data=linedata,size=0.1)+scale_x_continuous(breaks=div.x,labels=round(div.x,3))+scale_y_continuous(breaks=div.y,labels=round(div.x,3))+theme(panel.grid.minor=element_blank());
}
