#' Construct the projection pursuit classification tree 
#' 
#' Find tree structure using various projection pursuit indices of 
#' classification in each split.
#' @title Projection pursuit classification tree 
#' @usage PP.Tree.class(origclass, origdata, PPmethod="LDA",weight=TRUE,
#'                      r=1,lambda=0.1,energy=0,maxiter=50000, ...) 
#' @param origclass class information vector
#' @param origdata numeric data matrix without class information
#' @param PPmethod method for projection pursuit; "LDA", "PDA", "Lr", 
#'        "GINI", and "ENTROPY"
#' @param weight weight flag in LDA, PDA and Lr index
#' @param r r in Lr index
#' @param lambda lambda in PDA index
#' @param energy parameter for the probability to take new projection
#' @param maxiter maximum iteration number
#' @param ... arguments to be passed to methods
#' @return Tree.Struct tree structure of projection pursuit classification tree
#' @return projbest.node 1 dimensional optimal projections of each node split
#' @return splitCutoff.node cutoff values of each node split 
#' @return origclass original class 
#' @return origdata original data
#' @references Lee, YD, Cook, D., Park JW, and Lee, EK(2013) 
#' PPtree: Projection Pursuit Classification Tree, 
#' Electronic Journal of Statistics, 7:1369-1386.
#' @export
#' @keywords tree
#' @examples
#' data(iris)
#' Tree.result <- PP.Tree.class(iris[,5],iris[,1:4],"LDA")
#' Tree.result
PP.Tree.class<-function(origclass,origdata,PPmethod="LDA",weight=TRUE,r=1,
                        lambda=0.1,energy=0,maxiter=50000,...){
   TOL<-NULL
   origdata<-as.matrix(origdata)
   Find.proj<-function(origclass,origdata,PPmethod,weight,r,lambda,
                       maxiter,...){
      n<-nrow(origdata)
      p<-ncol(origdata)
      g<-table(origclass)
      g.name<-as.numeric(names(g))
      G<-length(g)
      origclass<-as.numeric(factor(origclass))
      if(PPmethod=="Lr"){
         indexbest<-Lrindex(origclass,as.matrix(origdata),weight=weight,r=r);
      } else if(PPmethod=="GINI"){
         indexbest<-0;
         for(i in 1:p){
            tempdata<-as.matrix(origdata[,i]);
            tempindex<-GINIindex1D(origclass,tempdata);  
            if(indexbest<tempindex)
               indexbest<-tempindex;
         }  
      } else if(PPmethod=="ENTROPY"){
         indexbest<-0;
         for(i in 1:p){
            tempdata<-as.matrix(origdata[,i]);
            tempindex<-ENTROPYindex1D(origclass,tempdata);  
            if(indexbest<tempindex)
               indexbest<-tempindex;
         }  
      } else{
        indexbest<-0;
      }
      energy<-ifelse(energy==0,1-indexbest,energy)
      energy.temp<-1-indexbest
      TOL<-energy.temp/1000000
  
      if(PPmethod=="LDA"){
         a<-LDAopt(as.numeric(as.factor(origclass)),origdata,weight,q=1)
      } else if(PPmethod=="PDA"){
         a<-PDAopt(as.numeric(as.factor(origclass)),origdata,weight,q=1,
                   lambda=lambda)   
      } else {    
         a<-PPopt(as.numeric(as.factor(origclass)),as.matrix(origdata),
                  weight,q=1,PPmethod=PPmethod,r=r,energy=energy,cooling=0.999,
                  TOL=TOL) 
      } 
      proj.data<-as.matrix(origdata)%*%a$projbest
      sign<-sign(a$projbest[abs(a$projbest)==max(abs(a$projbest))])
      index<-(1:p)*(abs(a$projbest)==max(abs(a$projbest)))
      index<-index[index>0]

      if(G==2){
         class<-origclass
      } else{
         m<-tapply(c(proj.data),origclass,mean)
         sd<-tapply(c(proj.data),origclass,sd)
         sd.sort<-sort.list(sd)
         m.list<-sort.list(m)
         m.sort<-sort(m)
         m.name<-as.numeric(names(m.sort))
         G<-length(m)
         dist<-0
         split<-0
         for(i in 1:(G-1)){
            if(m[m.list[i+1]]-m[m.list[i]]>dist){
               split<-i
               dist<-m[m.list[i+1]]-m[m.list[i]]
            }
         }
         class<-rep(0,n)
         for(i in 1:split) 
            class<-class+(origclass==m.name[i])
         class<-2-class
         g<-table(class)
         g.name<-as.numeric(names(g))
         G<-length(g)
         n<-nrow(origdata)
         class<-as.numeric(factor(class))
         if(PPmethod=="Lr"){
            indexbest<-Lrindex(class,as.matrix(origdata),weight=weight,r=r);
         } else if(PPmethod=="GINI"){
            indexbest<-0;
            for(i in 1:p){
               tempdata<-origdata[,i];
               tempindex<-GINIindex1D(class,as.matrix(tempdata));  
               if(indexbest<tempindex)
                  indexbest<-tempindex;
            }  
         } else if(PPmethod=="ENTROPY"){
            indexbest<-0;
            for(i in 1:p){
               tempdata<-origdata[,i];
               tempindex<-ENTROPYindex1D(class,as.matrix(tempdata));  
               if(indexbest<tempindex)
                  indexbest<-tempindex;
            }  
         } else {
           indexbest<-0
         }
         energy<-ifelse(energy==0,1-indexbest,energy)
         energy.temp<-1-indexbest
         TOL<-energy.temp/1000000
         if(PPmethod=="LDA"){
            a<-LDAopt(as.numeric(as.factor(class)),as.matrix(origdata),
                      weight,q=1)
         } else if(PPmethod=="PDA"){
            a<-PDAopt(as.numeric(as.factor(class)),as.matrix(origdata),
                      weight,q=1,lambda=lambda)   
         } else { 
            a<-PPopt(as.numeric(as.factor(class)),as.matrix(origdata),
                     PPmethod=PPmethod,r=r,q=1,energy=energy,cooling=0.999,
                     TOL=TOL)      
         }  
         if(sign!=sign(a$projbest[index])) 
            a$projbest<--a$projbest
         proj.data<-as.matrix(origdata)%*%a$projbest
      }
 
      m.LR<-tapply(proj.data,class,mean)
      temp.list<-sort.list(m.LR)
      m.LR<-m.LR[temp.list]
      sd.LR<-tapply(proj.data,class,function(x) 
                                    ifelse(length(x)>1,stats::sd(x),0))[temp.list]
      IQR.LR<-tapply(proj.data,class,function(x) 
                                    ifelse(length(x)>1,stats::IQR(x),0))[temp.list]
      median.LR<-tapply(proj.data,class,stats::median)[temp.list]
      n.LR<-table(class)[temp.list]
        
      c1<-(m.LR[1]+m.LR[2])/2
      c2<-(m.LR[1]*n.LR[2]+m.LR[2]*n.LR[1])/sum(n.LR)
      c3<-ifelse(sum(sd.LR==0)!=0,c1,
                        (m.LR[1]*sd.LR[2]+m.LR[2]*sd.LR[1])/sum(sd.LR))
      c4<-ifelse(sum(sd.LR==0)!=0,c2,(m.LR[1]*sd.LR[2]/sqrt(n.LR[2])+ 
                        m.LR[2]*sd.LR[1]/sqrt(n.LR[1]))/
                        (sd.LR[1]/sqrt(n.LR[1])+sd.LR[2]/sqrt(n.LR[2])))
      c5<-(median.LR[1]+median.LR[2])/2
      c6<-(median.LR[1]*n.LR[2]+median.LR[2]*n.LR[1])/sum(n.LR)    
      c7<-ifelse(sum(IQR.LR==0)!=0,c5,(median.LR[1]*IQR.LR[2]+
                                         median.LR[2]*IQR.LR[1])/sum(IQR.LR))      
      c8<-ifelse(sum(IQR.LR==0)!=0,c6,(median.LR[1]*(IQR.LR[2]/sqrt(n.LR[2]))+ 
                           median.LR[2]*(IQR.LR[1]/sqrt(n.LR[1])))/
                     ((IQR.LR[1]/sqrt(n.LR[1]))+(IQR.LR[2]/sqrt(n.LR[2]))))
      sel.proj<-sort(proj.data[which(proj.data>stats::quantile(proj.data,prob=0.25)& 
                                     proj.data<stats::quantile(proj.data,prob=0.75))])
      C<-c(c1,c2,c3,c4,c5,c6,c7,c8)
      Index<-a$indexbest
      Alpha<-t(a$projbest)
      IOindexR<-NULL
      IOindexL<-NULL
      sort.LR<-as.numeric(names(sort(m.LR)))
      IOindexL<-class==sort.LR[1]
      IOindexR<-class==sort.LR[2]
      list(Index=Index,Alpha=Alpha,C=C,IOindexL=IOindexL,IOindexR = IOindexR)
   }
   
   Tree.construct<-function(origclass,origdata,Tree.Struct, 
          id,rep,rep1,rep2,projbest.node,splitCutoff.node,PPmethod,r = NULL, 
          lambda=NULL,maxiter,...) {
      origclass<-as.integer(origclass)
      n<-nrow(origdata)
      g<-table(origclass)
      G<-length(g)
      if(length(Tree.Struct)==0){
         Tree.Struct<-matrix(1:(2*G-1),ncol=1)
         Tree.Struct<-cbind(Tree.Struct,0,0,0,0)
      }
      if(G==1){
         Tree.Struct[id,3]<-as.numeric(names(g))
         list(Tree.Struct=Tree.Struct,projbest.node=projbest.node, 
                splitCutoff.node=splitCutoff.node,rep=rep,rep1=rep1,rep2=rep2)
      } else {
         Tree.Struct[id,2]<-rep1
         rep1<-rep1+1
         Tree.Struct[id,3]<-rep1
         rep1<-rep1+1
         Tree.Struct[id,4]<-rep2
         rep2<-rep2+1
         a<-Find.proj(origclass,origdata,PPmethod,weight,r,lambda,
                      maxiter,...)
         splitCutoff.node<-rbind(splitCutoff.node,a$C)
         Tree.Struct[id,5]<-a$Index
         projbest.node<-rbind(projbest.node,a$Alpha)
         t.class<-origclass
         t.data<-origdata
         t.class<-t.class*a$IOindexL
         t.n<-length(t.class[t.class==0])
         t.index<-sort.list(t.class)
         t.index<-sort(t.index[-(1:t.n)])
         t.class<-t.class[t.index]
         t.data<-origdata[t.index,]
         b<-Tree.construct(t.class,t.data,Tree.Struct, 
                           Tree.Struct[id, 2],rep,rep1,rep2,projbest.node, 
                           splitCutoff.node,PPmethod,r,lambda,maxiter,...)
         Tree.Struct<-b$Tree.Struct
         projbest.node<-b$projbest.node
         splitCutoff.node<-b$splitCutoff.node
         rep<-b$rep
         rep1<-b$rep1
         rep2<-b$rep2
         t.class<-origclass
         t.data<-origdata
         t.class<-(t.class*a$IOindexR)
         t.n<-length(t.class[t.class==0])
         t.index<-sort.list(t.class)
         t.index<-sort(t.index[-(1:t.n)])
         t.class<-t.class[t.index]
         t.data<-origdata[t.index,]
         n<-nrow(t.data)
         G<-length(table(t.class))
         b<-Tree.construct(t.class,t.data,Tree.Struct, 
                Tree.Struct[id,3],rep,rep1,rep2,projbest.node, 
                splitCutoff.node,PPmethod,r,lambda,maxiter,...)
         Tree.Struct<-b$Tree.Struct
         projbest.node<-b$projbest.node
         splitCutoff.node<-b$splitCutoff.node
         rep<-b$rep
         rep1<-b$rep1
         rep2<-b$rep2
      }
      list(Tree.Struct=Tree.Struct,projbest.node=projbest.node, 
            splitCutoff.node=splitCutoff.node,rep=rep,rep1=rep1,rep2=rep2)
   }
   splitCutoff.node<-NULL
   projbest.node<-NULL
   Tree.Struct<-NULL
   id<-1
   rep1<-2
   rep2<-1
   rep<-1
   Tree.final<-Tree.construct(origclass,origdata,Tree.Struct,id,rep,rep1, 
                              rep2,projbest.node, splitCutoff.node,
                              PPmethod,r,lambda,TOL,maxiter,...)                            
   Tree.Struct<-Tree.final$Tree.Struct
   colnames(Tree.Struct)<-c("id","L.node.ID","R.F.node.ID","Coef.ID","Index")
   projbest.node<-Tree.final$projbest.node
   splitCutoff.node<-Tree.final$splitCutoff.node
   colnames(splitCutoff.node)<-paste("Rule",1:8,sep="")
   treeobj<-list(Tree.Struct=Tree.Struct,projbest.node=projbest.node, 
                 splitCutoff.node=splitCutoff.node,origclass=origclass,
                 origdata=origdata)
   class(treeobj)<-append(class(treeobj),"PPtreeclass")
   return(treeobj)
}