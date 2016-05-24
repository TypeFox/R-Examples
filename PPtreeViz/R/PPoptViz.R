#' Visualize PPopt result
#' 
#' Visualize the result of projection pursuit optimization
#' @title PPopt visualization
#' @usage PPopt.Viz(PPoptOBJ)
#' @param PPoptOBJ PPoptim object. result from LDAopt, PDAopt, and PPopt
#' @references Lee, EK., Cook, D., Klinke, S., and Lumley, T.(2005) 
#' Projection Pursuit for Exploratory Supervised Classification, 
#' Journal of Computational and Graphical Statistics, 14(4):831-846.
#' @export
#' @keywords tree
#' @examples
#' data(iris)
#' PPopt.Viz(LDAopt(iris[,5],iris[,1:4],q=1))
#' PPopt.Viz(LDAopt(iris[,5],iris[,1:4],q=2))
#' @import ggplot2 grid gridExtra 

PPopt.Viz<-function(PPoptOBJ){
   proj.data<-PPoptOBJ$origdata%*%PPoptOBJ$projbest
   q<-ncol(proj.data)
   p<-ncol(PPoptOBJ$origdata)   
   vID <-1:p
   if(q==1){
      ..density..<-NULL
      origclass<-PPoptOBJ$origclass
      plot.data<-data.frame(proj.data,origclass)
      p1<-ggplot(plot.data,aes(x=proj.data,group=origclass))+
                geom_histogram(aes(y=..density..,fill=origclass))
      coef<-PPoptOBJ$projbest[,1]
      coef.data<-data.frame(vID,coef)
      bin.width<-ifelse(p>100,1,0.1)
      y.max<-max(c(abs(coef.data$coef),1/sqrt(p)))
      y.min<- -y.max
      p2<-ggplot(coef.data,aes(x=vID,y=coef))
      if(p<=10){
         p2<-p2+geom_segment(aes(yend=0,xend=vID,size=1))
      } else{
         p2<-p2+geom_segment(aes(yend=0,xend=vID))
      }       
      p2<-p2+geom_hline(yintercept=0)+
            geom_hline(yintercept=c(-1,1)*1/sqrt(p),
                       col=2,linetype="dashed")+
            ylim(y.min,y.max) +
            xlab("variable ID")+
            ggtitle("Coefficients of Best Projection")+
            theme(legend.position = "none")
      gridExtra::grid.arrange(p2,p1,nrow=1)   
   } else{
      plot.list<-list()
      list.id<-1
      for(i in 1:q){
         for(j in 1:q){
            if(i==j){
               coef<-PPoptOBJ$projbest[,i]
               coef.data<-data.frame(vID,coef)
               bin.width<-ifelse(p>100,1,0.1)
               y.max<-max(c(abs(coef.data$coef),1/sqrt(p)))

               temp.plot<-ggplot(coef.data,aes(x=vID,y=coef))
               if(p<=10){
                  temp.plot<-temp.plot+geom_segment(aes(yend=0,xend=vID,size=1))
               } else{
                  temp.plot<-temp.plot+geom_segment(aes(yend=0,xend=vID))
               }       
               temp.plot<-temp.plot+
                          geom_hline(yintercept=c(-1,1)*1/sqrt(p),
                                     col=2,linetype="dashed")+
                          geom_hline(yintercept=0)+ 
                          xlab("variable ID")+ylim(-y.max,y.max)+
                          ggtitle(paste("Coefficients of Best Projection - dim",
                                         as.character(i),sep=""))+
                          theme(legend.position = "none")
               plot.list[[list.id]]<-temp.plot
               list.id<-list.id+1
            } else{
               x<-proj.data[,j]
               y<-proj.data[,i]
               origclass<-PPoptOBJ$origclass
               plot.data<-data.frame(x,y,origclass)
               plot.list[[list.id]]<-ggplot(plot.data,aes(x=x,y=y,
                                                          color=origclass))+
                                     geom_point()+
                                     xlab(paste("dim",as.character(j)))+
                                     ylab(paste("dim",as.character(i)))                 
               list.id<-list.id+1 
            }
         }    
      }
      do.call(grid.arrange,plot.list)
   }
}