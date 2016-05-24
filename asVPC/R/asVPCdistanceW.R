#' calculate percentiles of original data using distance-related weight 
#' percentiles of simulated data with corresponding confidence interval
#' 
#' @param orig.data the original data for model fitting
#' @param sim.data the simulated data from NONMEM
#' @param n.timebin the number of bin in X axis
#' @param n.sim the number of simulation in the simulated data
#' @param n.hist the number of shifted 
#' @param q.list numeric vector of probabilities with values in [0,1]
#' @param conf.level confidence level of the interval
#' @param X.name the name of X variable in the original scatter plot
#' @param Y.name the name of Y variable in the original scatter plot
#' @param opt.DV.point option to put data point in the plot
#' @param weight.flag option to use weight in average shifted calculation 
#' @param Y.min minimum of Y range in the plot
#' @param Y.max maximum of Y range in the plot
#' @param only.med option to use only median 
#' @param plot.flag TRUE: drawing plot / FALSE: generate data for drawing plot
#' @return plot or the values to draw plot
#' @export
#' @seealso \code{\link{asVPC.binW}}
#' @references new paper...
#' @author Eun-Kyung Lee \email{lee.eunk@@gmail.com}
#' @examples
#' data(origdata)
#' data(simdata)
#' asVPC.distanceW(origdata,simdata,n.timebin=10, n.sim=100,n.hist=3)
    
asVPC.distanceW<-function(orig.data,sim.data,n.timebin,n.sim,n.hist,
                          q.list=c(0.05,0.5,0.95),
                          conf.level=0.95,
                          X.name="TIME",Y.name="DV",
                          opt.DV.point=FALSE,
                          weight.flag=FALSE,
                          Y.min=NULL,
                          Y.max=NULL,
                          only.med=FALSE,
                          plot.flag=TRUE){
   SIM.CIarea.1<-NULL
   SIM.CIarea.2<-NULL
   SIM.CIarea.3<-NULL
   DV.point<-NULL
   DV.quant<-NULL
   SIM.quant<-NULL 
   ID<-NULL;G<-NULL
   bintot.N<-n.timebin*n.hist
   time.bin<-makeCOVbin(orig.data[,X.name],N.covbin=bintot.N)
   alpha<-1-conf.level     
   Q.CI<-vector("list",3)
   orig.Q<-NULL
   bintot.N<-nrow(time.bin$COV.bin.summary)
   for(i in 1:bintot.N){
      if(i<n.hist){
         sel.id<-which(as.numeric(time.bin$COV.bin)<=i+n.hist-1)
         sel.id1<-which(as.numeric(time.bin$COV.bin)==i)
         mid.point<-median(orig.data[sel.id1,X.name])
         low.point<-time.bin$COV.bin.summary$lower.COV[i]
         upper.point<-time.bin$COV.bin.summary$upper.COV[i]
      } else if(i >(bintot.N-n.hist+1)){
         sel.id<-which(as.numeric(time.bin$COV.bin)>=i-(n.hist-1))
         sel.id1<-which(as.numeric(time.bin$COV.bin)==i)
         mid.point<-median(orig.data[sel.id1,X.name])
         low.point<-time.bin$COV.bin.summary$lower.COV[i]
         upper.point<-time.bin$COV.bin.summary$upper.COV[i]
      } else{
         sel.id<-which(as.numeric(time.bin$COV.bin)>i-n.hist & 
                         as.numeric(time.bin$COV.bin)<i+n.hist)
         sel.id1<-which(as.numeric(time.bin$COV.bin)==i)
         mid.point<-median(orig.data[sel.id1,X.name])
         low.point<-time.bin$COV.bin.summary$lower.COV[i]
         upper.point<-time.bin$COV.bin.summary$upper.COV[i]
      }   
      if(bintot.N<length(table(orig.data$TIME))){
         dist.temp<-abs(orig.data$TIME[sel.id]-mid.point)
         temp.weight<-(max(dist.temp)-dist.temp)/diff(range(dist.temp))
      } else{
         A<-as.numeric(time.bin$COV.bin[sel.id])
         temp<-abs(A-median(range(A)))
         temp.weight<-(max(temp)+1)-temp
         temp.weight<-temp.weight/max(temp.weight)
      }
      if(weight.flag){
         temp.quantile<-t(apply(sim.data[sel.id,],2,function(x) 
           Hmisc::wtd.quantile(x,weight=temp.weight,
                                                     prob=q.list,na.rm=TRUE)))
         temp.orig.q<-Hmisc::wtd.quantile(orig.data[,Y.name][sel.id],
                                   weight=temp.weight,prob=q.list,na.rm=TRUE)
      } else{
         temp.quantile<-t(apply(sim.data[sel.id,],2,function(x) 
                                           quantile(x,prob=q.list,na.rm=TRUE)))
         temp.orig.q<-quantile(orig.data[,Y.name][sel.id],
                                                        prob=q.list,na.rm=TRUE)
      } 
      orig.Q<-rbind(orig.Q,c(mid.point,temp.orig.q))
      temp<-t(apply(temp.quantile,2,function(x) 
                         quantile(x,prob=c(alpha/2,0.5,1-alpha/2),na.rm=TRUE)))
      for(j in 1:length(q.list))
         Q.CI[[j]]<-rbind(Q.CI[[j]],c(mid.point,low.point,upper.point,temp[j,]))
   } 
   keep.name<-NULL
   for(j in 1:length(q.list)){
      keep.name<-c(keep.name,paste("Q",round(q.list[j]*100),"th",sep=""))
      colnames(Q.CI[[j]])<-c("mid","Lower","upper",colnames(Q.CI[[j]])[4:6])
   }    
   names(Q.CI)<-keep.name
   colnames(orig.Q)<-c("mid","Y1","Y2","Y3")
   orig.Q<-data.frame(orig.Q)
 
   plot.data<-data.frame(orig.data,X=orig.data[,X.name],Y=orig.data[,Y.name])
   if(is.null(Y.min)) Y.min<-min(c(plot.data$Y,Q.CI[[1]][,4]),na.rm=T)
   if(is.null(Y.max)) Y.max<-max(c(plot.data$Y,Q.CI[[length(Q.CI)]][,6]),na.rm=T)

   P.temp<-ggplot(plot.data,aes(x=X,y=Y))+ylim(Y.min,Y.max)+
              labs(x=X.name,y=Y.name)+theme_bw()+
              theme(panel.grid.major=element_line(colour="white"))+
              theme(panel.grid.minor=element_line(colour="white"))

   test.LU<-Q.CI[[1]][,2:3]
   test.data.tot<-Q.CI
   X.temp<-c(test.LU[,1],test.LU[nrow(test.LU),2])
   n.temp<-nrow(test.LU)
   X<-c(test.LU[1,1],rep(test.LU[2:n.temp,1],each=2),test.LU[n.temp,2])
   X<-c(X,X[length(X):1])
   if(!only.med){
      test.data<-test.data.tot[[1]]
      Y<-c(rep(test.data[,4],each=2),rep(test.data[(n.temp:1),6],each=2))
      SIM.CIarea.1<-data.frame(X=X,Y=Y,ID=1)
      P.temp<-P.temp+geom_polygon(data= SIM.CIarea.1,
                                  aes(x=X,y=Y,group=ID,fill=ID),
                                  fill="gray80",colour="gray80")
      test.data<-test.data.tot[[3]]
      Y<-c(rep(test.data[,4],each=2),rep(test.data[(n.temp:1),6],each=2))
      SIM.CIarea.3<-data.frame(X=X,Y=Y,ID=1)
      P.temp<-P.temp+geom_polygon(data= SIM.CIarea.3,
                                  aes(x=X,y=Y,group=ID,fill=ID),
                                  fill="gray80",colour="gray80")
   }
   test.data<-test.data.tot[[2]]
   Y<-c(rep(test.data[,4],each=2),rep(test.data[(n.temp:1),6],each=2))
   SIM.CIarea.2<-data.frame(X=X,Y=Y,ID=1)
   P.temp<-P.temp+geom_polygon(data= SIM.CIarea.2,
                               aes(x=X,y=Y,group=ID,fill=ID),
                               fill="gray50",colour="gray50")  
   if(opt.DV.point==TRUE){
      P.temp<-P.temp+geom_point(,color="grey30",size=2,alpha=0.5) 
      DV.point<-data.frame(X=orig.data[,X.name],Y=orig.data[,Y.name])
   }
     
   DV.quant<-data.frame(X=rep(orig.Q$mid,length(q.list)),
                        G=factor(rep(paste("Q",round(q.list*100),"th",sep=""),
                                      each=nrow(orig.Q))),
                        Y=unlist(orig.Q[,-1]))
   P.temp<-P.temp+geom_line(data=DV.quant[DV.quant$G!="Q50th",],
                            aes(x=X,y=Y,group=G),linetype=2,
                            size=1,color="black")+
                  geom_line(data=DV.quant[DV.quant$G=="Q50th",],
                            aes(x=X,y=Y,group=G),linetype=1,
                            size=1,color="black")
   colnames(orig.Q)<-c("X.mid",paste("Q",round(q.list*100),"th",sep="")) 
   if(plot.flag){
      P.temp
   } else{
      return(list(SIM.CIarea.1=SIM.CIarea.1,SIM.CIarea.2=SIM.CIarea.2, 
                  SIM.CIarea.3=SIM.CIarea.3,DV.point=DV.point, 
                  DV.quant=DV.quant,SIM.quant=SIM.quant))
   }
}