#' 
#' @name GOVenn
#' @title Venn diagram of differentially expressed genes.
#' @description The function compares lists of differentially expressed genes 
#'   and illustrates possible relations.Additionally it represents the variety 
#'   of gene expression patterns within the intersection in small pie charts 
#'   with three segements. Clockwise are shown the number of commonly up- 
#'   regulated, commonly down- regulated and contra- regulated genes.
#' @param data1 A data frame consisting of two columns: ID, logFC
#' @param data2 A data frame consisting of two columns: ID, logFC
#' @param data3 A data frame consisting of two columns: ID, logFC
#' @param title The title of the plot
#' @param label A character vector to define the legend keys
#' @param lfc.col A character vector determining the background colors of the 
#'   pie segments representing up- and down- regulated genes
#' @param circle.col A character vector to assign clockwise colors for the 
#'   circles
#' @param plot If TRUE only the venn diagram is plotted. Otherwise the function 
#'   returns a list with two items: the actual plot and a list containing the 
#'   overlap entries (default= TRUE)
#' @details The \code{plot} argument can be used to adjust the amount of 
#'   information that is returned by calling the function. If you are only 
#'   interested in the actual plot of the venn diagram, \code{plot} should be 
#'   set to TRUE. Sometimes you also want to know the elements of the 
#'   intersections. In this case \code{plot} should be set to FALSE and the 
#'   function call will return a list of two items. The first item, that can be 
#'   accessed by $plot, contains the plotting information. Additionally, a list
#'   ($table) will be returned containing the elements of the various overlaps. 
#' @import ggplot2
#' @examples
#' \dontrun{
#' #Load the included dataset
#' data(EC)
#' 
#' #Generating the circ object
#' circ<-circular_dat(EC$david, EC$genelist)
#' 
#' #Selecting terms of interest
#' l1<-subset(circ,term=='heart development',c(genes,logFC))
#' l2<-subset(circ,term=='plasma membrane',c(genes,logFC))
#' l3<-subset(circ,term=='tissue morphogenesis',c(genes,logFC))
#' 
#' GOVenn(l1,l2,l3, label=c('heart development','plasma membrane','tissue morphogenesis'))
#' }
#' @export

GOVenn<-function(data1, data2, data3, title, label, lfc.col, circle.col, plot=T){
  id <- NULL
  if (missing(label)) label<-c('List1','List2','List3')
  if (missing(lfc.col)) lfc.col<-c('firebrick1','gold','cornflowerblue')
  if (missing(circle.col)) circle.col<-c('brown1','chartreuse3','cornflowerblue')
  if (missing(title)) title<-''
  if (missing(data3)==F) {
    three<-T
    overlap<-get_overlap(data1,data2,data3)
    venn_df<-overlap$venn_df
    table<-overlap$table
  }else{
    three<-F
    overlap<-get_overlap2(data1,data2)
    venn_df<-overlap$venn_df
    table<-overlap$table
  }
  
	### calc Venn ###
  if (three){
    center<-data.frame(x=c(0.4311,0.4308,0.6380),y=c(0.6197,0.3801,0.5001),diameter=c(0.4483,0.4483,0.4483))
    outerCircle<-data.frame(x=numeric(),y=numeric(),id=numeric())
	  for (var in 1:3){
		  dat <- circleFun(c(center$x[var],center$y[var]),center$diameter[var],npoints = 100)
		  outerCircle<-rbind(outerCircle,dat)
	  }
	  outerCircle$id<-rep(c(label[1],label[2],label[3]),each=100)
	  outerCircle$id<-factor(outerCircle$id, levels=c(label[1],label[2],label[3]))
  }else{
    center<-data.frame(x=c(0.33,0.6699),y=c(0.5,0.5),diameter=c(0.6180,0.6180))
    outerCircle<-data.frame(x=numeric(),y=numeric(),id=numeric())
    for (var in 1:2){
      dat <- circleFun(c(center$x[var],center$y[var]),center$diameter[var],npoints = 100)
      outerCircle<-rbind(outerCircle,dat)
    }
    outerCircle$id<-rep(c(label[1],label[2]),each=100)
    outerCircle$id<-factor(outerCircle$id, levels=c(label[1],label[2]))
  }

	### calc single pies ### 
  if (three){
    Pie<-data.frame(x=numeric(),y=numeric(),id=numeric())
	  dat <- circleFun(c(center$x[1],max(subset(outerCircle,id==label[1])$y)-0.05),0.1,npoints = 100)
	  Pie<-rbind(Pie,dat)
	  dat <- circleFun(c(center$x[2],min(subset(outerCircle,id==label[2])$y)+0.05),0.1,npoints = 100)
	  Pie<-rbind(Pie,dat)
	  dat <- circleFun(c(max(subset(outerCircle,id==label[3])$x)-0.05,center$y[3]),0.1,npoints = 100)
	  Pie<-rbind(Pie,dat)
	  Pie$id<-rep(1:3,each=100)
	  UP<-Pie[c(1:50,100:150,200:250),]
	  Down<-Pie[c(50:100,150:200,250:300),]
  }else{
    Pie<-data.frame(x=numeric(),y=numeric(),id=numeric())
    dat <- circleFun(c(min(subset(outerCircle,id==label[1])$x)+0.05,center$y[1]),0.1,npoints = 100)
    Pie<-rbind(Pie,dat)
    dat <- circleFun(c(max(subset(outerCircle,id==label[2])$x)-0.05,center$y[2]),0.1,npoints = 100)
    Pie<-rbind(Pie,dat)
    Pie$id<-rep(1:2,each=100)
    UP<-Pie[c(1:50,100:150),]
    Down<-Pie[c(50:100,150:200),]
  }

	### calc single pie text ###
	if (three){
    x<-c();y<-c()
	  for (i in unique(Pie$id)){
		  x<-c(x,rep((min(subset(Pie,id==i)$x)+max(subset(Pie,id==i)$x))/2,2))
		  y<-c(y,(min(subset(Pie,id==i)$y)+max(subset(Pie,id==i)$y))/2+0.02)
		  y<-c(y,(min(subset(Pie,id==i)$y)+max(subset(Pie,id==i)$y))/2-0.02)
	  }
	  pieText<-data.frame(x=x,y=y,label=c(venn_df$UP[1],venn_df$DOWN[1],venn_df$UP[2],venn_df$DOWN[2],venn_df$UP[3],venn_df$DOWN[3]))
  }else{
    x<-c();y<-c()
    for (i in unique(Pie$id)){
      x<-c(x,rep((min(subset(Pie,id==i)$x)+max(subset(Pie,id==i)$x))/2,2))
      y<-c(y,(min(subset(Pie,id==i)$y)+max(subset(Pie,id==i)$y))/2+0.02)
      y<-c(y,(min(subset(Pie,id==i)$y)+max(subset(Pie,id==i)$y))/2-0.02)
    }
    pieText<-data.frame(x=x,y=y,label=c(venn_df$UP[1],venn_df$DOWN[1],venn_df$UP[2],venn_df$DOWN[2]))
  }
  
	### calc overlap pies ### 
  if (three){
    smc<-data.frame(x=c(0.6,0.59,0.31,0.5),y=c(0.66,0.34,0.5,0.5))
    PieOv<-data.frame(x=numeric(),y=numeric())
	  PieOv<-rbind(PieOv,circleFun(c(smc$x[1],smc$y[1]),0.06,npoints = 100))
	  PieOv<-rbind(PieOv,circleFun(c(smc$x[2],smc$y[2]),0.06,npoints = 100))
	  PieOv<-rbind(PieOv,circleFun(c(smc$x[3],smc$y[3]),0.06,npoints = 100))
	  PieOv<-rbind(PieOv,circleFun(c(smc$x[4],smc$y[4]),0.06,npoints = 100))
	  PieOv$id<-rep(1:4,each=100)
	  smc$id<-1:4
	  UPOv<-rbind(smc[1,],PieOv[1:33,],smc[1,],smc[2,],PieOv[100:133,],smc[2,],smc[3,],PieOv[200:233,],smc[3,],smc[4,],PieOv[300:333,],smc[4,])
	  Change<-rbind(smc[1,],PieOv[33:66,],smc[1,],smc[2,],PieOv[133:166,],smc[2,],smc[3,],PieOv[233:266,],smc[3,],smc[4,],PieOv[333:366,],smc[4,])
	  DownOv<-rbind(smc[1,],PieOv[66:100,],smc[1,],smc[2,],PieOv[166:200,],smc[2,],smc[3,],PieOv[266:300,],smc[3,],smc[4,],PieOv[366:400,],smc[4,])
  }else{
    PieOv<-data.frame(x=numeric(),y=numeric(),id=numeric())
    PieOv<-rbind(PieOv,circleFun(c(0.5,0.5),0.08,npoints = 100))
    PieOv$id<-rep(1,100)
    center<-data.frame(x=0.5, y=0.5, id=1)
    UPOv<-rbind(center[1,],PieOv[1:33,])
    Change<-rbind(center[1,],PieOv[33:66,])
    DownOv<-rbind(center[1,],PieOv[66:100,])
  }
  
  ### calc overlap pie text ###
  if (three){
    x<-c();y<-c()
	  for (i in unique(PieOv$id)){
		  x<-c(x,subset(UPOv,id==i)$x[1]+0.0115,subset(DownOv,id==i)$x[1]-0.018,subset(Change,id==i)$x[1]+0.01)
		  y<-c(y,subset(UPOv,id==i)$y[1]+0.01,subset(DownOv,id==i)$y[1],subset(Change,id==i)$y[1]-0.013)
	  }
	  small.pieT<-data.frame(x=x,y=y,label=c(venn_df$UP[5],venn_df$Change[5],venn_df$DOWN[5],venn_df$UP[6],venn_df$Change[6],venn_df$DOWN[6],venn_df$UP[4],venn_df$Change[4],venn_df$DOWN[4],venn_df$UP[7],venn_df$Change[7],venn_df$DOWN[7]))
  }else{
    x<-c(subset(UPOv,id==1)$x[1]+0.015,subset(DownOv,id==1)$x[1]-0.018,subset(Change,id==1)$x[1]+0.01)
    y<-c(subset(UPOv,id==1)$y[1]+0.015,subset(DownOv,id==1)$y[1],subset(Change,id==1)$y[1]-0.013)
    small.pieT<-data.frame(x=x,y=y,label=c(venn_df$UP[3],venn_df$Change[3],venn_df$DOWN[3]))
  }
  
	g<- ggplot()+
    geom_polygon(data=outerCircle, aes(x,y, group=id, fill=id) ,alpha=0.5,color='black')+
	  scale_fill_manual(values=circle.col)+
    guides(fill=guide_legend(title=''))+
	  geom_polygon(data=UP, aes(x,y,group=id),fill=lfc.col[1],color='white')+
	  geom_polygon(data=Down, aes(x,y,group=id),fill=lfc.col[3],color='white')+
	  geom_text(data=pieText, aes(x=x,y=y,label=label),size=5)+
    geom_polygon(data=UPOv, aes(x,y,group=id),fill=lfc.col[1],color='white')+
	  geom_polygon(data=DownOv, aes(x,y,group=id),fill=lfc.col[3],color='white')+
	  geom_polygon(data=Change, aes(x,y,group=id),fill=lfc.col[2],color='white')+
	  geom_text(data=small.pieT,aes(x=x,y=y,label=label),size=4)+
    theme_blank+
    labs(title=title)
  
  if (plot) return(g) else return(list(plot=g,table=table))
}


	
