#' Visualization tools to explore each node of PPtree
#' 
#' For the inner node, four plots are provided - the bar chart style plot with
#' projection pursuit coefficients of each variable, the histogram of the
#' projected data, the bar chart style plots with means of each variables for 
#' the left and the right group, and the image plot of correlations. 
#' @title PPtree node visualization
#' @usage PPclassNode.Viz(PPclassOBJ,node.id,Rule,legend,std,image,diff.prop)
#' @param PPclassOBJ PPregclass object
#' @param node.id node ID
#' @param Rule cutoff rule
#' @param legend flag to represent legend in the plot. Default value is TRUE
#' @param std flag to standardize data before drawing plot
#' @param image flag to draw image plot of correlation matrix
#' @param diff.prop percentage of number of variables with significant 
#'                  differences and shown in red in the bar chart style means
#' @references Lee, YD, Cook, D., Park JW, and Lee, EK(2013) 
#' PPtree: Projection Pursuit Classification Tree, 
#' Electronic Journal of Statistics, 7:1369-1386.
#' @export
#' @keywords tree
#' @examples
#' data(iris)
#' Tree.result <- PP.Tree.class(iris[,5],iris[,1:4],"LDA")
#' Tree.result
#' PPclassNode.Viz(Tree.result,1,1)
PPclassNode.Viz<-function(PPclassOBJ,node.id,Rule,
                          legend=TRUE,std=TRUE,image=FALSE,diff.prop=0.2){
   searchGroup<-function(node.id,TS,gName){
      flag<-TRUE
      sel.id<-TS[node.id,2:3]
      LR.id<-c(TRUE,FALSE)
      sel.group<-NULL
      i<-1
      while((sel.id[i]!=0)&&(i<length(sel.id))){
         if(TS[sel.id[i],2]!=0){
            sel.id<-c(sel.id,TS[sel.id[i],2:3])
            if(LR.id[i])
               LR.id<-c(LR.id,c(TRUE,TRUE))
            else
               LR.id<-c(LR.id,c(FALSE,FALSE))
         }   
         if(TS[sel.id[i+1],2]!=0){
            sel.id<-c(sel.id,TS[sel.id[i+1],2:3])
            if(LR.id[i+1])
              LR.id<-c(LR.id,c(TRUE,TRUE))
            else
              LR.id<-c(LR.id,c(FALSE,FALSE))
         }
         i<-i+2
      }
      sel.Name<-TS[sel.id[which(TS[sel.id,2]==0)],3]
      selName<-sort(gName[sel.Name])
      L.list<-sort(gName[sel.Name[LR.id[which(TS[sel.id,2]==0)]]])
      R.list<-sort(gName[sel.Name[!LR.id[which(TS[sel.id,2]==0)]]])     
      return(list(selName=selName,Llist=L.list,Rlist=R.list))
   }
   
   TS<-PPclassOBJ$Tree.Struct
   Alpha<-PPclassOBJ$projbest.node
   cut.off<-PPclassOBJ$splitCutoff.node
   origdata<-PPclassOBJ$origdata
   origclass<-PPclassOBJ$origclass
   p<-ncol(origdata)
   gName<-names(table(origclass))
   if(TS[node.id,2]!=0){
      SG.result<-searchGroup(node.id,TS,gName)
      selG<-SG.result$selName
      selL<-SG.result$Llist
      selR<-SG.result$Rlist      
      sel.id<-NULL
      LR.class<-NULL
      for(i in 1:length(selG)){
         sel.id<-c(sel.id,which(origclass==selG[i])) 
         LR.class<-c(LR.class,
                     rep(ifelse(sum(selL==selG[i])!=0,"L","R"),
                         length(which(origclass==selG[i]))))
      }

      proj.data<-c(as.matrix(origdata)%*%
                     as.matrix(Alpha[TS[node.id,4],]))[sel.id]
      proj.class<-origclass[sel.id]
      plot.data<-data.frame(proj.data=proj.data,origclass=proj.class)
      p1<- ggplot(plot.data, aes(x=proj.data,group=origclass))+
                geom_histogram( aes(fill=origclass),position="stack")+
                geom_vline(xintercept=cut.off[TS[node.id,4],Rule],
                           linetype="longdash",lwd=1,col=2)
      if(!legend) p1<-p1+theme(legend.position="none")
      vID <-1:p
      coef<-Alpha[TS[node.id,4],]
      coef.data<-data.frame(vID=vID,coef=coef)
      bin.width<-ifelse(p>100,1,0.1)
      y.max <-max(c(abs(coef.data$coef),1/sqrt(p)))
      
      p2<-ggplot(coef.data,aes(x=vID,y=coef))
      if(p<=10){
        p2<-p2+geom_segment(aes(yend=0,xend=vID,size=1))
      } else{
        p2<-p2+geom_segment(aes(yend=0,xend=vID))
      }
      p2<-p2+geom_hline(yintercept=0)+ 
          geom_hline(yintercept=c(-1,1)*1/sqrt(ncol(origdata)),
                     col=2,linetype="dashed")+
          xlab("variable ID")+ggtitle(paste("Node",node.id,sep=" "))+
          ylim(-y.max,y.max)+
          theme(legend.position="none")
      sel.data<-origdata[sel.id,]
      if(std){
         sel.data<-apply(sel.data,2,function(x) (x-mean(x))/stats::sd(x))
         ytitle<-"adjusted mean by each variable mean"
      } else{
         ytitle<-"adjusted mean by overall mean"
      }
      temp.data<-c(apply(sel.data[LR.class=="L",],2,mean),
                   apply(sel.data[LR.class!="L",],2,mean))
      if(!std) temp.data<-temp.data-mean(temp.data)
      vvID<-1;mean.data<-1;Var1<-1;Var2<-1;value<-1;   
      L.data<-temp.data[1:p];R.data<-temp.data[-(1:p)]
      LRcolor<-rep("NonSig",2*p)
      diff.LR<-abs(L.data-R.data)
      cutoff.LR<-stats::quantile(diff.LR,prob=1-diff.prop)

      LRcolor[c(diff.LR>cutoff.LR,diff.LR>cutoff.LR)]<-"Sig"
      plot.data2<-data.frame(mean.data=temp.data,
                             vvID=c(vID,vID),
                             LR=factor(c(rep("L",p),rep("R",p))),
                             LRcolor=factor(LRcolor,levels=c("Sig","NonSig")))
      y.max3 <-max(c(abs(plot.data2$mean.data)))

      p3<-ggplot(plot.data2,aes(x=vvID,y=mean.data,color=LRcolor))
      if(p<=10){
         p3<-p3+geom_segment(aes(yend=0,xend=vvID,size=1))
      } else{
         p3<-p3+geom_segment(aes(yend=0,xend=vvID))
      }
      p3<-p3+facet_grid(LR~.)+
             ylab(ytitle)+xlab("variable ID")+
             ggtitle("Mean of left and right nodes")+ylim(-y.max3,y.max3)+  
             geom_hline(yintercept=0)+
             theme(legend.position="none")
      if(image & p<=30){
         image.cor<-stats::cor(sel.data)
         colnames(image.cor)<-paste("V",1:ncol(image.cor),sep="")
         rownames(image.cor)<-paste("V",1:nrow(image.cor),sep="")
         temp.data<-data.frame(Var1=rep(colnames(image.cor),nrow(image.cor)),
                          Var2=rep(colnames(image.cor),each=nrow(image.cor)),
                          value=c(image.cor))
         p4<-ggplot(temp.data,aes(x=Var1,y=Var2,fill=value))+
         geom_tile()+
             scale_fill_gradient(low ="blue",high="yellow",limit=c(-1,1))+
             xlab("variables")+ylab("variables")+
             ggtitle("correlation matrix")+theme(aspect.ratio=1)
         gridExtra::grid.arrange(p2,p1,p3,p4,nrow=2)
      } else{
         gridExtra::grid.arrange(p2,p3,p1,nrow=1)
      }  
   } else{
      sel.id<-which(PP.classify(PPclassOBJ,origdata,Rule)$predict.class==
                                                          gName[TS[node.id,3]])
      find.i<-TS[which((TS[,2]==node.id|TS[,3]==node.id)& TS[,2]!=0),4]
      proj.data<-c(as.matrix(origdata)%*%as.matrix(Alpha[find.i,]))[sel.id]
      proj.class<-origclass[sel.id]
      plot.data<-data.frame(proj.data=proj.data,proj.class=proj.class)
      
      p1<-ggplot(plot.data,aes(x=proj.data,group=proj.class))+
           geom_histogram(aes(fill=proj.class),position="stack")+
           ggtitle(paste("Node",node.id,": ",gName[TS[node.id,3]],sep=""))  
      gridExtra::grid.arrange(p1,nrow=1)     
   }
}