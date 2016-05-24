# classify.r: Functions for Record Linkage with machine learning algorithms.

trainSupv <- function(rpairs,method,use.pred=FALSE,omit.possible=TRUE,
                convert.na=TRUE, include.data=FALSE, ...)
{
  # catch erronous input
  if (!("RecLinkData" %in% class(rpairs) ||
    "RecLinkResult" %in% class(rpairs)))
    stop("Wrong class for rpairs!")
  # error if only matches or non-matches are present
#  if(!isTRUE(use.pred) && length(unique(rpairs$pairs$is_match))==1)
#    stop("Training examples have the same match status!")
#  if(isTRUE(use.pred) && length(unique(rpairs$prediction))==1)
#    stop("Training examples have the same match status!")


  # in case of bumping, refer to backend function (this
  # method differs in terms of data handling)
  if (method=="bumping") return(.bumping(rpairs, use.pred=use.pred,
    omit.possible=omit.possible, convert.na=convert.na,
    include.data=include.data, ...))

	pairs=rpairs$pairs[,-c(1:2)]
	if (convert.na)
		pairs[is.na(pairs)]=0
	if (use.pred)
    {
    	if (is.null(rpairs$prediction))
    		stop ("No prediction vector in rpairs! Call with use.pred=FALSE.")
		pairs$is_match=rpairs$prediction
	} else
	{
		pairs$is_match=factor(rep("P",length(pairs$is_match)),
		  levels=c("N","P","L"))
		pairs$is_match[rpairs$pairs$is_match==1]="L"		
		pairs$is_match[rpairs$pairs$is_match==0]="N"		
	}
	# delete possible links if desired
	if (omit.possible)
		pairs=pairs[pairs$is_match!="P",,drop=FALSE]
  # drop unused levels
  pairs$is_match <- factor(pairs$is_match)

  # now check if a usable training set remains, i.e. there are
  # at least two distinct examples with different match outcome
  if (nrow(unique(pairs[,-ncol(pairs)])) < 2)
    stop("Not enough distinct training examples!")
  if (length(unique(pairs$is_match))==1)
    stop("All training examples have the same match status!")
  

	model=switch(as.character(method[1]),
		svm=svm(is_match ~ .,data=pairs,type="C-classification",...),
		rpart=rpart(is_match ~ .,data=pairs,method="class",...),
		ada=ada(is_match ~ .,data=pairs,...),
		bagging=bagging(is_match ~ .,data=pairs,method="class",...),
		nnet=nnet(is_match ~ .,data=pairs,size=ncol(pairs)*2, ...),
		stop("Illegal method"))
  ret=list()
 	if (isTRUE(include.data))
	   ret$train=rpairs
  ret$attrNames <- colnames(rpairs$pairs)[-c(1,2,ncol(rpairs$pairs))]
  ret$model=model
  ret$method=as.character(method)[1]
	class(ret)="RecLinkClassif"
    return(ret)

}

# classifySupv: moved to classifySupv-methods.r


classifyUnsup <- function(rpairs, method,...)
{
  if (!("RecLinkData" %in% class(rpairs) ||
    "RecLinkResult" %in% class(rpairs)))
    stop("Wrong class for rpairs!")

  if (nrow(rpairs$pairs) < 2)
    stop("Not enough record pairs!")
  
	if (method=="kmeans" || method=="bclust")
	{
		x=as.matrix(rpairs$pairs[,-c(1,2,ncol(rpairs$pairs))])
		x[is.na(x)]=0
		clust=switch(method,
			kmeans=kmeans(x,centers=2,...),
			bclust=bclust(x,...))
		y=clust$cluster
		# mark links and non-links. The cluster farther from 0 is
		# interpreted as link cluster
		link=ifelse(sum(clust$centers[1,])>sum(clust$centers[2,]),1,2)
		rpairs$prediction=rep("N",length(y))
		rpairs$prediction[y==link]="L"
    # refactor to ensure uniform order of levels
    rpairs$prediction=factor(rpairs$prediction,levels=c("N","P","L"))
		class(rpairs)="RecLinkResult"
		return(rpairs)	
	}
	stop("Illegal method!")
}


# backend function for bumping
.bumping <- function(rpairs, n.bootstrap=25, cp=0.01, ...)
{
  # Determine best cp value
#  cpfit=trainSupv(rpairs,method="rpart",cp=0,minsplit=2,xval=10)
#  cp=cpfit$model$cptable[which.min(cpfit$model$cptable[,"xerror"]),"CP"]
#
  # List of individual trees and error measures
  classif_list=list()
  performance_list=list()
  n_data <- nrow(rpairs$pairs)
  for (i in 1:n.bootstrap)
  {
    boot=rpairs[sample(n_data, replace=TRUE)]
    classif_list[[i]]=trainSupv(boot,"rpart",cp=cp,minsplit=2,xval=0, ...)
    performance_list[[i]]=errorMeasures(classifySupv(classif_list[[i]],rpairs))$accuracy
  }
  # add whole data for comparison
  classif_list[[n.bootstrap+1]]=trainSupv(rpairs,"rpart",cp=cp,minsplit=2,xval=0)
  performance_list[[n.bootstrap+1]]=errorMeasures(classifySupv(
      classif_list[[n.bootstrap+1]],rpairs))$accuracy
  return(classif_list[[which.max(performance_list)]])
}
