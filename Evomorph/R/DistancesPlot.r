DistancesPlot <- function(dist){
  Generations<-Procrustes_Distance<-Sample<-NULL
  if(is.list(dist)){
    df<-data.frame(x=1:length(dist[[1]]),y=dist)
    colnames(df)<-c("Generations",1:length(dist)-1)
   }
  else{
    df<-data.frame(x=1:length(dist),y=dist)
    colnames(df)<-c("Generations",0)
  }
  data<-melt(df,id="Generations")
  colnames(data)<-c("Generations","Sample","Procrustes_Distance")
  dplot=ggplot(data=data,aes(x=Generations,y=Procrustes_Distance,colour=Sample,group = Sample))+geom_point()+geom_line()+theme_bw()
  return(dplot)
}

