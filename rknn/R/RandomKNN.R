################################################################################
# RandomKNN program                                                            #
# File:   RandomKNN.R                                                          #
# Author: Shengqiao Li                                                         #
# Date:   June 24, 2008 (initial)                                              #
# Dependency: class, Hmisc                                                     #
# Change Log:                                                                  #
#           December 11, 2008 - Add set.return.seed                            #
#           January 31, 2009 - Add randomKNN.cv                                #
#           May, 2013 -- Parallelization                                       #  
################################################################################

#imputation method in R
#impute::impute.knn

confusion<- function(obs, pred)
{
  #when print out, "classified as" actually points to pred.
  table(obs, pred, dnn=list("classified as->", ""));

}

confusion2acc<- function(ct)
{
  sum(diag(ct))/sum(ct)
}

cv.coef<- function(x)
#coefficient of variation
{
  sd(x)/mean(x)
}



rknn<- function(data, newdata, y, k=1, r=500, 
              mtry=trunc(sqrt(ncol(data))),
              cluster=NULL, seed=NULL
)
{
  #
  #Note: kknn::kknn is not good as class::knn. Low accuracy for Golub data
  #      klaR::sknn is even worse
  # 
   
   knns<- function(data, newdata, y, k, r, mtry)
   {
     p<-ncol(newdata);
     n<-nrow(newdata);
     
     selected<- matrix(integer(), nrow=r, ncol=mtry);
     pred.all<- matrix(nrow=n, ncol=r);
     
     for(j in 1:r){
  
        fset<- sample(p, mtry);
        aknn<- knn(train=data[, fset], test=newdata[, fset], cl=y, k=k);
       
        selected[j,]<- fset;
        pred.all[,j]<- as.integer(aknn);
    }
  
    return(list(selected=selected, pred.all= pred.all));   
   }
  
   res<- list(call=match.call());
   
#disable 
#   res$Random.seed<- set.return.seed(Random.seed, seed);
   
   p<- ncol(newdata);
   
   n<- nrow(newdata);

   res$k<- k;
   res$r<- r;
   res$mtry<- mtry;
   res$n<- n;
   res$p<- p;

   if(!is.factor(y)) y<- as.factor(y);

   if(!is.null(cluster)){
      getLoadedDLLs();
      
      clusterSetRNGStream(cluster, seed);      
      clusterExport(cluster, c('knns', 'data', 'newdata', 'y', 'k', 'mtry'), envir=environment());
      
      cluster.result<- clusterApply(cluster, splitR(r, length(cluster)), function(r) knns(data=data, newdata=newdata, y=y, k=k, r=r, mtry=mtry));
      selected<- matrix(unlist(lapply(cluster.result, function(x)t(x$selected))), ncol=mtry, byrow=TRUE);
      pred.all<- matrix(unlist(lapply(cluster.result, function(x)x$pred.all)), nrow=n, byrow=FALSE);
      
   }
   else{
      set.seed(seed);

      result<- knns(data=data, newdata=newdata, y=y, k=k, r=r, mtry=mtry);
      selected<- result$selected;
      pred.all<- result$pred.all;
  }

  pred<- character(n);
  for(i in 1:n) pred[i]<- names(which.max(table(pred.all[i,])));


  res$pred<- factor(pred, levels = seq_along(levels(y)), labels = levels(y));
  
   #features<- table(selected);

  #names(features)<- colnames(data)[as.integer(names(features))];
  res$features<- if(is.null(colnames(data))){1:p} else colnames(data);
  	
  res$features.used<- selected;
  
  class(res)<- "rknn";

  return(res);
}


##cross-validation

rknn.cv<- function(data, y, k=1, r=500, 
          mtry=trunc(sqrt(ncol(data))),
          cluster=NULL, seed=NULL
          )
{   

  knns.cv<- function(data, y, k, r, mtry)
  {
     p<-ncol(data);
     n<-nrow(data);
     
     selected<- matrix(integer(), nrow=r, ncol=mtry);
     pred.all<- matrix(nrow=n, ncol=r);
     
     for(j in 1:r){
  
        fset<- sample(p, mtry);
        aknn<- knn.cv(train=data[, fset], cl=y, k=k);
       
        selected[j,]<- fset;
        pred.all[,j]<- as.integer(aknn);
    }
  
    return(list(selected=selected, pred.all= pred.all));   
  }          
  
        
   res<- list(call=match.call());


   p<- ncol(data);
   n<- nrow(data);

   res$k<- k;
   res$r<- r;
   res$mtry<- mtry;
   res$n<- n;
   res$p<- p;

   if(!is.factor(y)) y<- as.factor(y);

  if(!is.null(cluster)){
      clusterSetRNGStream(cluster, seed);      
      clusterExport(cluster, c('knns.cv', 'data', 'y', 'k', 'mtry'), envir=environment());
      
      cluster.result<- clusterApply(cluster, splitR(r, length(cluster)), function(r) knns.cv(data=data, y=y, k=k, r=r, mtry=mtry));
      selected<- matrix(unlist(lapply(cluster.result, function(x)t(x$selected))), ncol=mtry, byrow=TRUE);
      pred.all<- matrix(unlist(lapply(cluster.result, function(x)x$pred.all)), nrow=n, byrow=FALSE);
      
   }
   else{
      set.seed(seed);

      result<- knns.cv(data=data, y=y, k=k, r=r, mtry=mtry);
      selected<- result$selected;
      pred.all<- result$pred.all;
  }
          
  pred<- character(n);
  for(i in 1:n) pred[i]<- names(which.max(table(pred.all[i,])));


  res$pred<- factor(pred, levels = seq_along(levels(y)), labels = levels(y));

   #features<- table(selected);

  #names(features)<- colnames(data)[as.integer(names(features))];
  res$features<- if(is.null(colnames(data))){1:p} else colnames(data);

  res$features.used<- selected;

  class(res)<- "rknn";

  return(res);
}

varUsed<- function (x, by.KNN = FALSE, count = TRUE)
{
    if (!inherits(x, "rknn"))
        stop(deparse(substitute(x)), "is not a rknn object")

    if (is.null(x$features.used))
        stop(deparse(substitute(x)), "does not contain variables")
    p <- x$p
    if (count) {
        if (by.KNN) {
            v <- apply(x$features.used, 1, function(x) {
                xx <- integer(p)
                y <- table(x)
                xx[as.integer(names(y))] <- y
                xx
            })
            v<- t(v)
        }
        else {
            v <- integer(p)
            vv <- table(x$features.used)
            v[as.integer(names(vv))] <- vv
        }
   }
   else {
        v <- t(apply(x$features.used, 1, function(x) sort(x)));
        if (!by.KNN)
            v <- sort(unique(as.vector(v)));
   }

   v;

}

varNotUsed<- function(x)
{
   if (!inherits(x, "rknn"))
        stop(deparse(substitute(x)), "is not a rknn object")

	index<- setdiff(1:x$p, varUsed(x, by.KNN=FALSE, count=FALSE))

	if(length(index)==0){ NULL} else x$features[index];

}
