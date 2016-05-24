################################################################################
# Random KNN Backward Elimination Functions                                    #
# File:   backward_selection.R                                                 #
# Author: Shengqiao Li                                                         #
# Date:   June 24, 2008 (initial)                                              #
#         2013-08-05 parallelization                                           #
################################################################################

rknnSupport<- function(data, y, k=1, r=500,
                      mtry=trunc(sqrt(ncol(data))),
                      fixed.partition=FALSE,
                      cluster=NULL,
                      seed=NULL
)
{
	 knns<- function(data, y, k, r, mtry, index, ns, fixed.partition=FALSE)
	 {	 
	   selected<- matrix(0, nrow=r, ncol=mtry);
	   acc<- numeric(r); #accuracy for each randomKNN

     pred<- vector("list", length=nrow(data));
     if(fixed.partition){
     		dset<- unlist(mapply(sample, index, ns/2));   #bipartition each class
        tset<- setdiff(1:n, dset);
  	 }
  	 for(i in 1:r){
  	   	if(!fixed.partition){
          	dset<- unlist(mapply(sample, index, ns/2));   #bipartition each class
          	tset<- setdiff(1:n, dset);
      	}
          fset<- sample(p, mtry);
          selected[i,]<- fset;
          aknn<- knn(train=data[dset, fset], test=data[tset, fset], cl=y[dset], k=k);
  
          #convert to character for comparison in case  that there are fewer factor leves from aknn
          acc[i]<- mean(as.character(aknn)==as.character(y[tset]))
  
          for(j in 1:length(tset)){
            pred[[tset[j]]]<- c(pred[[tset[j]]], aknn[j])
          }
     }
    return(list(selected=selected, acc=acc, pred=pred)); 
	 }
	 
   res<- list(call=match.call());
   
   res$k<- k;
   res$r<- r;
   res$mtry<- mtry;

   n<- nrow(data);
   p<- ncol(data);
   res$n<- n;
   res$p<- p;

   if(!is.factor(y)) y<- as.factor(y);
   nclass<- nlevels(y); #number of classes
   index<- apply(sapply(y, "==", levels(y)), 1, which); #indice for each class

   if(is.matrix(index)) index<- as.list(data.frame(index));  #For equal sample size, index is a matrix

   ns<-  sapply(index, length) ;  #class size

  if(is.null(cluster))
  {
     set.seed(seed);
      
     result<- knns(data=data, y=y, k=k, r=r, mtry=mtry, index=index, ns=ns, fixed.partition=fixed.partition);
     selected<- result$selected;
     acc<- result$acc;
     pred<- result$pred;
   }
  else{
      clusterSetRNGStream(cluster, seed);      
      clusterExport(cluster, c('knns', 'data', 'y', 'k', 'mtry', 'index', 'ns'), envir=environment());
      
      cluster.result<- clusterApply(cluster, splitR(r, length(cluster)), 
        function(r) knns(data=data, y=y, k=k, r=r, mtry=mtry, index=index, ns=ns, fixed.partition=fixed.partition));

      acc<- unlist(lapply(cluster.result, function(x)x$acc));    
      
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
	    
 		  for(j in 1:length(tset))  pred[[tset[j]]]<- names(which.max(table(pred[[tset[j]]])));

      res$pred<- factor(unlist(pred), levels = seq_along(levels(y)), labels = levels(y));

  		res$accuracy<- mean(res$pred==y[tset]);

  		res$confusion<- confusion(obs=y[tset], pred=res$pred);

	}
  else {
  		for(i in 1:n) pred[[i]]<- names(which.max(table(pred[[i]]))); #majority voting

  		res$pred<- factor(unlist(pred), levels = seq_along(levels(y)), labels = levels(y));

  		res$accuracy<- mean(res$pred==y);

  		res$confusion<- table("classified as->"=y, res$pred);  #when print out, "classified as" actually points to pred.
	}

  selected<- as.vector(selected);

  features<- unique(selected);

  mreal<- length(features);
  support<- numeric(mreal);

  for(j in 1:mreal){
      acc.ind<- (which(selected==features[j])-1)%%r+1;
      support[j]<- mean(acc[acc.ind]);
  }

  names(support)<- colnames(data)[features];

  res$support<-  sort(support, decreasing=T);

  res$meanacc<-  mean(acc); #mean accuracy of r randomKNN

  class(res)<- "rknnSupport";

  return(res);
}

rknnBeg<-function(data, y, k=1, r=500,
				mtry=trunc(sqrt(ncol(data))),
				fixed.partition=FALSE,
                pk=0.5,             # elimination proportion
                stopat=4,           # minimum number of variables
                cluster=NULL,
                seed=NULL           # an integer seed
)
#backward elimination in geometric fashion
{

    res <- list(call = match.call())
    
    if (is.numeric(y))  #regression
        rknnSupport <- rknnRegSupport
    p <- ncol(data)
    ni <- trunc(log(stopat/p)/log(pk))
    varlist <- list()
    acc <- numeric(ni + 1)
    rknn <- rknnSupport(data, y, k, r, fixed.partition=fixed.partition, cluster=cluster, seed=seed)
    varlist[[1]] <- rknn$support
    acc[[1]] <- rknn$meanacc
    for (i in 2:(ni + 1)) {
        newp <- round(p * pk)
        vars <- names(rknn$support)[1:newp]
        p <- newp
        rknn <- rknnSupport(data[, vars], y, k, r, fixed.partition=fixed.partition, cluster=cluster, seed=seed)
        varlist[[i]] <- rknn$support
        acc[i] <- rknn$meanacc
    }
    res$vars <- varlist
    res$p <- sapply(varlist, length);
    res$mean_support<- sapply(varlist, mean);
    res$mean_accuracy <- acc

    class(res) <- c("rknnBeg", "rknnBE")
    return(res)
}
rknnBel<- function(data, y, k=1, r=500,
			mtry=trunc(sqrt(ncol(data))),
			fixed.partition=FALSE,
            d=1,  #drop d features each time
            stopat=4,
            cluster=NULL,
            seed=NULL
)
## backward elimination linearly
{

    res <- list(call = match.call())
    
    if (is.numeric(y)) 
        rknnSupport <- rknnRegSupport
    p <- ncol(data)
    ni <- trunc((p - stopat)/d)
    varlist <- list()
    acc <- numeric(ni + 1)
    rknn <- rknnSupport(data, y, k, r, fixed.partition=fixed.partition, cluster=cluster, seed=seed)
    varlist[[1]] <- rknn$support
    acc[1] <- rknn$meanacc
    for (i in 2:(ni + 1)) {
        newp <- p - d
        vars <- names(rknn$support)[1:newp]
        p <- newp
        rknn <- rknnSupport(data[, vars], y, k, r, fixed.partition=fixed.partition, cluster=cluster, seed=seed)
        varlist[[i]] <- rknn$support
        acc[i] <- rknn$meanacc
    }
    res$vars <- varlist
    res$p <- sapply(varlist, length);
    res$mean_support<- sapply(varlist, mean);
    res$mean_accuracy <- acc
    class(res) <- c("rknnBel", "rknnBE")
    return(res)
}

prebestset<- function(x, criterion=c("mean_accuracy", "mean_support"))
{
  #prebest subset of varibles
  #x begKNN or belKNN object.
  criterion<- match.arg(criterion);
  
  names(x$vars[[max(which.max(x[[criterion]])-1, 1)]])
}
bestset<- function(x, criterion=c("mean_accuracy", "mean_support"))
{
  #best subset of varibles
  #x begKNN or belKNN object.
  criterion<- match.arg(criterion);
  names(x$vars[[which.max(x[[criterion]])]])
}
################################################################################