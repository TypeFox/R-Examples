################################################################################
# RandomKNN Regression Functions                                               #
# File:   RandomKNN_regression.R                                               #
# Author: Shengqiao Li                                                         #
# Date:   June 24, 2008 (initial)                                              #
# 2013-5-14 parallelization                                                    #
################################################################################
pressresid<- function(obj)
{ 
  #internal function
  #PRESS residuals of a linear model
  
  h<- influence(obj)$hat;
  resid(obj)/(1-h)

}
predicted<- function(obj)
{
  #internal function
  #hat Y(i,-i) of a linear model
  h<- influence(obj)$hat;
  fitted(obj) - resid(obj)*h/(1-h);

}
PRESS<- function(obj)
{
  #PRESS staitic
  if (inherits(obj, "knnReg"))  resi<- obj$residuals 
  else  if (inherits(obj, "lm"))  resi<- pressresid(obj)
  else   stop(deparse(substitute(obj)), " is not a linear model or knn regression object.");
 
 return(sum(resi^2));
}

rsqp<- function(obj)
{
  #R-square predicted of a linear model
  y<- obj$model[,1];
  SST<- sum((y-mean(y))^2);
  1 -PRESS(obj)/SST;
}

#knn.reg<- function (train, test = NULL, y, k = 3, check.FNN = TRUE)
#  #FNN::knn.reg is 40~400 times faster than dprep for p=400 and p=5
#  #50 times faster than kknn::kknn.
#{
#    if (check.FNN && require(FNN, quietly = T))
#        return(FNN::knn.reg(train, test, y, k))
#    else if (require(dprep, quietly = T)) {
#    #using knn imputation method to do knn regression (mean of neighbor response)
#        if (is.null(test)) {
#            notest <- TRUE
#            test <- train
#        }
#        else notest <- FALSE
#        if (is.vector(train)) { #univariate
#            train <- cbind(train)
#            test <- cbind(test)
#        }
#        else if (is.vector(test))  test <- rbind(test)
#        n <- ifelse(is.null(dim(test)), length(test), nrow(test))   #number of samples to be predict
#        pred <- rep(NA, n)
#        data <- rbind(cbind(test, pred), cbind(train, y))
#        
#        res <- ec.knnimp(data, k = k) #imputation function
#        
#        pred <- res[1:n, ncol(res)] #prediction in the first n rows of the last column
#        residuals <- if (notest)   pred - y
#        else NULL #otherwise don't know observed valaues
#        
#        R2 <- if (notest) {
#            1 - sum(residuals^2)/sum((y - mean(y))^2)
#        }
#        else NULL
#
#        res <- list(call = match.call(), k = k, n = n, pred = pred, residuals = residuals,  R2_pred = R2)
#        class(res)<- "knnReg";
#        
#        return(res)
#    }
#    else stop("Either FNN or dprep package is required!\n")
#}
#
#knn.reg.cv<- function(data, y, k=3)
#{
#   if(require(FNN, quietly=T)) return(knn.reg(train=data, y=y, k=k, algorithm="VR"));
#
#   if(is.vector(data)) data<- cbind(data);
#   n<- nrow(data);
#   pred<- numeric(n);
#   for(i in 1:n) pred[i]<- knn.reg(train=data[-i,], test=data[i,], y=y[-i], k=k, check.FNN=FALSE)$pred
#
#   residuals<- y-pred;
#   PRESS<- sum(residuals^2);
#   R2Pred <- 1- PRESS/sum((y-mean(y))^2);
#
#   res<- list(call=match.call(), n=n, y=y, pred=pred, residuals=residuals, PRESS=PRESS, R2Pred=R2Pred);
#   class(res)<- "knnRegCV";
#   return(res);
#
#}
#

rknnReg<- function(data, newdata, y, k=1, r=500,
                   mtry=trunc(sqrt(ncol(data))), 
                   cluster=NULL,
                   seed=NULL
)
{

  knns<- function(data, newdata, y, k, r, mtry)
  {
     p<-ncol(newdata);
     n<-nrow(newdata);
     
     selected<- matrix(integer(), nrow=r, ncol=mtry);
     pred.all<- matrix(nrow=n, ncol=r);
     
     for(j in 1:r){
  
        fset<- sample(p, mtry);
        aknn<- knn.reg(train=data[, fset], test=newdata[, fset], y=y, k=k);
       
        selected[j,]<- fset;
        pred.all[,j]<- aknn$pred;
    }
  
    return(list(selected=selected, pred.all= pred.all));   
  }
   
   res<- list(call=match.call());


   res$k<- k;
   res$r<- r;
   res$mtry<- mtry;

   n<- nrow(newdata);
   p<- ncol(newdata);
   res$p<- p;

   if(!is.null(cluster)){
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
      
  pred<- rowMeans(pred.all);

  res$pred<- pred;
  #res$residuals<- pred-y;

  #features<- table(selected);

  #names(features)<- colnames(data)[as.integer(names(features))];
  res$features<- if(is.null(colnames(data))){1:p} else colnames(data);

  res$features.used<- selected;

  class(res)<- "rknn";

  return(res);
}

rknnRegSupport<- function(data, y, k=k, r=500,
                     mtry=trunc(sqrt(ncol(data))), 
                     fixed.partition=FALSE, 
                     cluster=NULL,
                     seed=NULL
)
{    
   	knns<- function(data, y, k, r, mtry, fixed.partition=FALSE)
	 {	 
	   selected<- matrix(0, nrow=r, ncol=mtry);
     SS<- numeric(r); #Sum of square errors for each randomKNN
     n<- nrow(data);
     p<-  ncol(data);
     
     pred<- vector("list", length=n);

    if(fixed.partition){
   		  dset<-  sample(n, n/2);   #bipartition
        tset<- setdiff(1:n, dset);
	   }
   
  	 for(i in 1:r){
  	   	if(!fixed.partition){
     		   dset<-  sample(n, n/2);   #bipartition
           tset<- setdiff(1:n, dset);
      	}
        fset<- sample(p, mtry);
        selected[i,]<- fset;
        aknn<- knn.reg(train=data[dset, fset], test=data[tset, fset], y=y[dset], k=k)$pred;
  
        SS[i]<- sum((aknn-y[tset])^2);

        for(j in 1:length(tset)){
          pred[[tset[j]]]<- c(pred[[tset[j]]], aknn[j])
        }
        for(j in 1:length(tset)){
          pred[[tset[j]]]<- c(pred[[tset[j]]], aknn[j])
        }
      
     }
    return(list(selected=selected, SS=SS, pred=pred)); 
	 }
   
   res<- list(call=match.call());

   res$k<- k;
   res$r<- r;
   res$mtry<- mtry;

   n<- nrow(data);
   p<- ncol(data);

   SSt<- sum((y-mean(y))^2);

  if(is.null(cluster)){
      set.seed(seed);
      
     result<- knns(data=data, y=y, k=k, r=r, mtry=mtry, fixed.partition=fixed.partition);
     selected<- result$selected;
     SS<- result$SS;
     pred<- result$pred;
   } 
   else{
      clusterSetRNGStream(cluster, seed);      
      clusterExport(cluster, c('knns', 'data', 'y', 'k', 'mtry'), envir=environment());
      
      cluster.result<- clusterApply(cluster, splitR(r, length(cluster)), 
        function(r) knns(data=data, y=y, k=k, r=r, mtry=mtry, fixed.partition=fixed.partition));

      SS<- unlist(lapply(cluster.result, function(x)x$SS));    
      
      selected<- matrix(unlist(lapply(cluster.result, function(x)t(x$selected))), ncol=mtry, byrow=TRUE);
      
      pred<- vector("list", length=n);
      for(i in 1:n){
        tmp<- vector("list", length=length(cluster));    
        for(j in 1:length(cluster))
        {
            tmp[[j]]<- cluster.result[[j]]$pred[[i]]      
        }      
        pred[[i]]<- unlist(tmp);
      }  
   }

	if(fixed.partition){
    tset<- which(!sapply(pred, is.null));
    
 		for(j in 1:length(tset))  pred[[tset[j]]]<- mean(pred[[tset[j]]]);

  		res$pred<- unlist(pred);

  		res$SS<- sum((res$pred-y[tset])^2);
      res$accuracy<- 1- res$SS/SSt; #R-square
  		res$confusion<- rbind("Observed"=y[tset], "Predicted"=res$pred);
	}
  else{
  		res$pred<- sapply(pred, mean);

  		res$SS<- sum((res$pred-y)^2);
    	res$accuracy<- 1- res$SS/SSt; #R-square
  		res$confusion<- rbind("Observed"=y, "Predicted"=res$pred);

	}

  selected<- as.vector(selected);

  features<- unique(selected);

  mreal<- length(features);
  support<- numeric(mreal);

  for(j in 1:mreal){
      SS.ind<- (which(selected==features[j])-1)%%r+1;
      support[j]<- 1-mean(SS[SS.ind])/SSt;   #R-square
  }

  names(support)<- colnames(data)[features];

  res$support<-  sort(support, decreasing=T);

  res$meanacc<-  1-mean(SS)/SSt;

  class(res)<- "rknnSupport";

  return(res);
}
################################################################################
