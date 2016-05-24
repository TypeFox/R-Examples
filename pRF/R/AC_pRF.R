#Appease the gods of CRAN

m<-FDR<-value<-variable<-`.`<-Feature.id<-Importance<-NULL
#' @import multtest reshape2 ggplot2 permute dplyr
#' @importFrom randomForest randomForest importance
 

#PermutationRandomForest

#' @export
pRF<-function(response,predictors,n.perms,alpha=0.05,mtry=NULL,type=c("classification","regression"),ntree=500,seed=12345,...){
  
  #Set seed
  
  set.seed(seed)
  
  #depending on class check response variable class
  
  message("Checking if response variable is valid")
  
  if(type=="classification"){response<-as.factor(response)}
  if(type=="regression"){response<-as.numeric(as.character(response))}
  if(!type %in% c("classification","regression")){stop ("type parameter required to proceed, rerun function with type specified
                                                        types must be either regression or classification")}
  
  #Match up dimensions and transpose if required

  if(!nrow(predictors)==length(response)){ pid<-rownames(predictors)
  } else {
    
    pid<-colnames(predictors)
    
  }
  
  #Autotranspose
  
  if(!nrow(predictors)==length(response)){predictors=t(predictors)}
  
  #Fit original random forest
  
  message("fitting original random forest")
  RF.data<-randomForest(y=response,x=(predictors),mtry=mtry,ntree=ntree,importance=TRUE,...)
  
  if(type=="classification"){ 
    
    Observed<- importance( RF.data, type=2)  } else {
      
      Observed<-importance( RF.data, type=1)
      
    }
  
  
  #Generate permuted random forest  
  
  permutations.set<-shuffleSet(n=length(response),nset=n.perms)
  
  
  
  #Begin populating a frame of importances
  
  message("populating permuted importance table")
  
  
  list.vec<-list()
  
  for( i in 1:nrow(permutations.set)){
    
    if(type=="classification"){
      
      list.vec[[i]]<- randomForest(y=factor(response[permutations.set[i,]]),x=(predictors),mtry=mtry,ntree=ntree,importance=TRUE)%>%
        importance(.,type=2)%>%data.frame(.)%>%.$MeanDecreaseGini;print(i)
    } else {
      
      list.vec[[i]]<- randomForest(y=response[permutations.set[i,]],x=(predictors),mtry=mtry,ntree=ntree,importance=TRUE)%>%
        importance(.,type=1)%>%data.frame(.)%>%.$X.IncMSE;print(i)
      
    }
    
  }
  
  perms.imp<-do.call(cbind,list.vec)
  
  #Put everything together into a list
  
  results.list<-list()
  results.list$perms<-perms.imp
  results.list$obs<-Observed
  results.list$Model<-RF.data
  
  #Estimate significance
  
  #Here, I will use the exact p value given by Smyth's paper on permutation p values
  #This accounts for the biased nature of just seeing how many test stats as extreme are observed
  
  #That paper places the exact p value at b+1/m+1    
  #Calculate b, and m , and then estimate permutation p value
  
  message("calculating permutation p.values")
  
  
  b<-as.numeric()
  for(i in 1:nrow(Observed)) {
    
    b[[i]]<-length(which(perms.imp[i,] > Observed[i,]))
  }
  
  #Put together results table
  
  Res.table<-data.frame(b=b+1,m=ncol(perms.imp)+1)
  Res.table<-transform(Res.table,p.value=b/m)
  
  message("adjusting for FDR using two-step BH")
  
  FDR<-mt.rawp2adjp(rawp=Res.table$p.value,proc="TSBH",alpha=alpha)
 
  fdr.vec<-FDR$adjp[,2]
  index<-FDR$index
  Res.table$FDR<-fdr.vec[order(index)]
  
  Res.table$Feature.id<-pid
  
  results.list$Res.table<-Res.table
  return(results.list)
  }





#Plots for significant features, relative to null distribution

#' @export
sigplot<-function(pRF.list,threshold=0.05){
    
  #Create dataframe
  
  colnames(pRF.list$perms)<-make.unique(rep("perms",ncol(pRF.list$perms)))
  df<-cbind(pRF.list$Res.table,pRF.list$obs,pRF.list$perms)
  
  df<-melt(df,id.vars=c(1:6),measure.vars=c(7:ncol(df)))
  colnames(df)[[6]]<-"Importance"
  df<-filter(df,FDR < threshold)
  if(nrow(df)==0){stop ("No significant features below threshold")}
  #Locate points for observed value
  
  count.tab<-df%>%mutate(value=round(value))%>%count(Feature.id,value)
  
  qplot(data=df,geom='freqpoly',binwidth=1,x=value,facets=~Feature.id,alpha=I(.6),colour=I("dodgerblue3"))+
    geom_point(data=df,aes(x=Importance,y=max(count.tab$n/2)),size=I(3),colour=I("red"))+
    theme(panel.grid.minor=element_blank())
  
  
}
