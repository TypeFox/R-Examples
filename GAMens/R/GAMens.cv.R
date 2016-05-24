GAMens.cv <- 
function (formula, data, cv, rsm_size=2, autoform=FALSE, iter=10, df=4, bagging=TRUE, rsm=TRUE, fusion="avgagg") 
{
	setdiff.data.frame <- function(A,B) A[ !duplicated( rbind(B,A) )[ -seq_len(nrow(B))] , ]
	n <- dim(data)[1]
	if(cv <2 | cv>n) {stop("The number of cross-validations should be larger than 2 and smaller than the number of observations")}
	
	bootstrap <- sample(1:n,replace=FALSE)
	obs_per_fold <- floor(n/cv)
	pred_per_fold <- as.data.frame(array(NA, c(n,cv)))
	pred <- as.data.frame(array(NA, c(n,1)))
	class_per_fold <- as.data.frame(array(NA, c(n,cv)))
	class <- as.data.frame(array(NA, c(n,1)))

    	for (c in 1:cv) {
        	if (c < cv) {fold <- data[bootstrap[(((c-1)*obs_per_fold)+(1:obs_per_fold))],]
		} else {fold <- data[bootstrap[(((c-1)*obs_per_fold):nrow(data))],]}
        	if (c < cv) {fold_ids <- bootstrap[(((c-1)*obs_per_fold)+(1:obs_per_fold))]
		} else {fold_ids <- bootstrap[(((c-1)*obs_per_fold):nrow(data))]}
		traindata_folds <- setdiff.data.frame(data,fold)
		GAMens_object <- GAMens(formula,traindata_folds, rsm_size, autoform=autoform,iter=iter,df=df,bagging=bagging,rsm=rsm,fusion=fusion)
		results <- predict.GAMens(GAMens_object,fold)
		pred_per_fold[fold_ids,c] <- results[[1]]
		names(pred_per_fold)[[c]] <- paste("fold",c,sep="")
		pred[fold_ids,1] <- results[[1]]
		class_per_fold[fold_ids,c] <- results[[2]]
		names(class_per_fold)[[c]] <- paste("fold",c,sep="")
		class[fold_ids,1] <- results[[2]]
    	}

	conf <- table(as.matrix(class), data[,as.character(formula[[2]])], dnn=c("Predicted Class", "Observed Class"))
	output<- list(foldpred=pred_per_fold,pred=pred,foldclass=class_per_fold,class=class,conf=conf)
}