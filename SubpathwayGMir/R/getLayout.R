getLayout <-
function(graph){
   if(length(V(graph)$graphics_x)==0||length(V(graph)$graphics_y)==0) return (NULL)
    x_y<-c()
    graphics_x<-get.vertex.attribute(graph,"graphics_x")
    index<-which(graphics_x=="unknow")
	
    if(length(index)>1){
       temp<-as.numeric(graphics_x[which(graphics_x!="unknow")])
	   if(length(temp)<2){temp<-as.numeric(c(100,600))}
       replace_value<-seq(min(temp),max(temp),by = (max(temp)-min(temp))/(length(index)-1))
       graphics_x<-replace(graphics_x,which(graphics_x=="unknow"),replace_value)
    }else if(length(index)==1){
       temp<-as.numeric(graphics_x[which(graphics_x!="unknow")])
       graphics_x<-replace(graphics_x,which(graphics_x=="unknow"),min(temp))
    } 
    graphics_x <-as.numeric(graphics_x)
	
	graphics_y<-get.vertex.attribute(graph,"graphics_y")
    index<-which(graphics_y=="unknow")
    if(length(index)>0){
       temp<-as.numeric(graphics_y[which(graphics_y!="unknow")])
	   if(length(temp)<2){temp<-as.numeric(c(100,600))}
       graphics_y<-replace(graphics_y,which(graphics_y=="unknow"),max(temp)+100)
    } 
    graphics_y <-as.numeric(graphics_y)
	
    x_y<-as.matrix(data.frame(graphics_x=graphics_x, graphics_y=graphics_y))
    x_y[,2]<--x_y[,2]
	dimnames(x_y)<-NULL
	return (x_y)
}
