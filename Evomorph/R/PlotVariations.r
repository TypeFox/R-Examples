PlotVariations<- function(shapes,reference,path,save=F, w = 6, h = 6, d = 100, c_ref="red", c_target="black", s_ref=3, s_target=2){
  x<-y<-xmin<-ymin<-xmax<-ymax<-NULL
  if (!(is.list(shapes)))
  {stop("Shapes must be list type")}
  if (!(is.matrix(reference)))
  {stop("Reference landmarks must be matrix type")}
  if (length(shapes[[1]])!=length(reference))
  {stop("Shapes and reference must be the same length")}
  
  ref<-data.frame(x=reference[,1], y=reference[,2])
  
  #Calculate min and max values for reference
  xmin<-min(ref$x)
  xmax<-max(ref$x)
  ymin<-min(ref$y)
  ymax<-max(ref$y)
  
  newplots<-list()
  
  for(i in 1:length(shapes)){
  
  ext<-shapes[[i]]
  vari<-data.frame(x=ext[,1], y=ext[,2])
  ti<-paste("Target No.", i)
  
  newplots[[i]]<-ggplot(data=ref,aes(x=x,y=y))+ 
  geom_point(size=s_ref, color=c_ref)+coord_fixed()+
  geom_point(data=vari,aes(x=x,y=y),color=c_target,size=s_target)+
  geom_segment(data=vari,aes(x = x, y = y, xend = ref$x, yend = ref$y),show.legend = T)+
  theme_bw()+ggtitle(ti)+xlab("")+ylab("")+
  theme(axis.text=element_blank(), axis.ticks=element_blank())
  
  
  #Calculate max and min scale
  if(min(vari$x)<xmin){xmin<-min(vari$x)}
  if(max(vari$x)>xmax){xmax<-max(vari$x)}
  
  if(min(vari$y)<ymin){ymin<-min(vari$y)}
  if(max(vari$y)>ymax){ymax<-max(vari$y)}
  
  }
  
  for (i in 1:length(newplots)){newplots[[i]]<-newplots[[i]]+xlim(xmin,xmax)+ylim(ymin,ymax)}
  
  if (save==T){
  if (missing(path))
  {time=format(Sys.time(),"%a_%H_%M")
  path <- file.path(getwd(),time)
  print("Path not provided. Default path created on current working directory")
  dir.create(path)}
  for(g in 1:length(shapes)){
  if(g<10){ggsave(paste(path,"/model","000",g, ".jpeg", sep = ""),plot=newplots[[g]], width = h, height = w, units = 'in', dpi = d)}
  if(g<100 && g>=10){ggsave(paste(path,"/model","00",g, ".jpeg", sep = ""),plot=newplots[[g]], width = h, height = w, units = 'in', dpi = d)}
  if(g>=100){ggsave(paste(path,"/model","00",g, ".jpeg", sep = ""),plot=newplots[[g]], width = h, height = w, units = 'in', dpi = d)}
  }
  }
  
  return (newplots)

}