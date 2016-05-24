##################################################################### 
# This function provides a visual representation of the EM Algorithm.
#
###################################################################

plot<-function(index,x,means,covariance,clusters)
{

  x<-as.matrix(x)
  x <- cbind(x[,index[1]],x[,index[2]])
  
  clusters<-as.factor(clusters)
  
  if(!(is.matrix(means)&&is.list(covariance)) )
  {
    stop("Incorrect input type for mean or covariance")
  }
  ###The next line of code is solely for the package check
  y<-V1<-V2<-NULL
  
  data<-cbind(x,clusters)
  
  
  dat<-as.data.frame(data) 
  colnames(dat)<-c("V1","V2","clusters")
  
  
  dat_ell <- data.frame() 

  
  #contour plot
  for(g in as.integer(levels(clusters))) {dat_ell <- rbind(dat_ell, 
                                                           cbind(as.data.frame(with(dat[clusters==g,], ellipse(covariance[[g]],centre=c(means[g,])) ) ),clusters=g))}
  
   
  colnames(dat_ell)<-c("x","y","clusters")
    
  plot<-0
  if(index[1]==index[2])
  {
    plot<- ggplot(data=dat, aes(x=V1, y=V2,color=factor(clusters),shape=factor(clusters))) + geom_blank()
  }
  else
  {
    plot<-ggplot(data=dat, aes(x=V1, y=V2,color=factor(clusters),shape=factor(clusters))) + geom_point() +
    geom_path(data=dat_ell, aes(x=x, y=y,color=factor(clusters)))+scale_shape_manual(name = "Clusters",values = dat$clusters)+
    scale_colour_discrete(name = "Clusters")
  
     plot<-plot +theme(axis.title.x=element_blank(),axis.title.y=element_blank())
  }
 
    return(plot)
}


cluster_plot<-function(x,means,covariance,clusters)
{
  
  
  d<-ncol(x)
  if(d>5||d<2)
  {
    stop("This function can only be used for dimension equal to 2 and less than 5")
  }
  combine<-matrix(nrow=2,ncol=d**2)
  combine[1,]<-sapply(1:d,function(x,y) rep(x,y),y=d)
  combine[2,]<-rep(c(1:d),d)
  plots<-list()
  count<-1
  for(u in 1:ncol(combine))
  {
    variables<-combine[,u]
    variables<-sort(variables)
    mean<- means[,variables]
    sigma<- lapply(1:length(covariance),function(x,y,z) rbind(y[[x]][z[[1]],z],y[[x]][z[[2]],z]),y=covariance,z=variables)
    plots[[count]]<- plot(variables,x,mean,sigma,clusters)
    count<-count+1
  }
  
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(d,d)))
  
  count<-1
  for(p in seq(1,d))
  {
    
    for(q in seq(1,d))
    {
      
      print(plots[[count]],vp= viewport(layout.pos.row=p,layout.pos.col=q))
      count=count+1;
    }
  }
               
               
               
}
