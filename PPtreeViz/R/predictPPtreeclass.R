#' predict projection pursuit classification tree
#' 
#' Predict class for the test set with the fitted projection pursuit 
#' classification tree and calculate prediction error.
#' @title predict PPtree
#' @param object a fitted object of class inheriting from "PP.Tree.class" 
#' @param newdata  the test dataset
#' @param Rule split rule 1: mean of two group means 
#'                        2: weighted mean of two group means 
#'                           - weight with group size
#'                        3: weighted mean of two group means 
#'                           - weight with group sd
#'                        4: weighted mean of two group means 
#'                           - weight with group se
#'                        5: mean of two group medians 
#'                        6: weighted mean of two group medians 
#'                           - weight with group size
#'                        7: weighted mean of two group median 
#'                           - weight with group IQR
#'                        8: weighted mean of two group median 
#'                           - weight with group IQR and size                                         
#' @param ... arguments to be passed to methods
#' @aliases predict
#' @references Lee, YD, Cook, D., Park JW, and Lee, EK(2013) 
#' PPtree: Projection Pursuit Classification Tree, 
#' Electronic Journal of Statistics, 7:1369-1386.
#' @export
#' @keywords tree
#' @examples
#' data(iris)
#' n <- nrow(iris)
#' tot <- c(1:n)
#' n.train <- round(n*0.9)
#' train <- sample(tot,n.train)
#' test <- tot[-train]
#' Tree.result <- PP.Tree.class(iris[train,5],iris[train,1:4],"LDA")
#' predict(Tree.result)
predict.PPtreeclass<-function(object,newdata=NULL,Rule=1,...) {
   Tree.result<-object
   if(is.null(newdata))
      newdata<-Tree.result$origdata
   test.data<-as.matrix(newdata)
   PP.Classification<-function(Tree.Struct,test.class.index,IOindex,
                               test.class,id,rep){
      if(Tree.Struct[id,4]==0){
         i.class<-test.class
         i.class[i.class>0]<-1
         i.class<-1-i.class
         test.class<-test.class+IOindex*i.class*Tree.Struct[id, 3]
         return(list(test.class=test.class,rep=rep))
      } else{  
         IOindexL<-IOindex*test.class.index[rep,]
         IOindexR<-IOindex*(1-test.class.index[rep,])
         rep<-rep+1
         a<-PP.Classification(Tree.Struct,test.class.index,IOindexL,
                              test.class,Tree.Struct[id,2],rep)
         test.class<-a$test.class
         rep<-a$rep;
         a<-PP.Classification(Tree.Struct,test.class.index,IOindexR,
                              test.class,Tree.Struct[id,3],rep)
         test.class<-a$test.class
         rep<-a$rep
      }
      list(test.class=test.class,rep=rep)
   }
  
   PP.Class.index<-function(class.temp,test.class.index,test.data,
                               Tree.Struct,Alpha.Keep,C.Keep,id,Rule){
      class.temp<-as.integer(class.temp)
      if(Tree.Struct[id,2]==0){
         return(list(test.class.index=test.class.index,class.temp=class.temp))
      } else{
         t.class<-class.temp 
         t.n<-length(t.class[t.class==0])
         t.index<-sort.list(t.class)
         if(t.n)
            t.index<-sort(t.index[-(1:t.n)])
         t.data<-test.data[t.index,]
         id.proj<-Tree.Struct[id,4]
            
         proj.test<-as.matrix(test.data)%*%as.matrix(Alpha.Keep[id.proj,])
         proj.test<-as.double(proj.test)
         class.temp<-t(proj.test<C.Keep[id.proj,Rule]) 
         test.class.index<-rbind(test.class.index,class.temp)
         a<-PP.Class.index(class.temp,test.class.index,test.data,
                           Tree.Struct,Alpha.Keep,C.Keep,
                           Tree.Struct[id,2],Rule)
         test.class.index<-a$test.class.index
         a<-PP.Class.index(1-class.temp,test.class.index,test.data,
                           Tree.Struct,Alpha.Keep,C.Keep,
                           Tree.Struct[id,3],Rule)
         test.class.index<-a$test.class.index;
      }
      list(test.class.index=test.class.index,class.temp=class.temp)
   }
    
   n<-nrow(test.data)
   class.temp<-rep(1,n)
   test.class.index<-NULL
   temp<-PP.Class.index(class.temp,test.class.index,test.data,
                        Tree.result$Tree.Struct,Tree.result$projbest.node,
                        Tree.result$splitCutoff.node,1,Rule)
   test.class<-rep(0,n)
   IOindex<-rep(1,n)
   temp<-PP.Classification(Tree.result$Tree.Struct,temp$test.class.index,
                           IOindex,test.class,1,1)
   class.name<-names(table(Tree.result$origclass))
   predict.class<-factor(class.name[temp$test.class])
   return(predict.class)
}

