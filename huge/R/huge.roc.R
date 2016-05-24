#-----------------------------------------------------------------------#
# Package: High-dimensional Undirected pathraph Estimation              #
# huge(): Draw ROC Curve for a solution path                            #
#         The ground truth is required                                  #
# Authors: Tuo Zhao and Han Liu                                         #
# Emails: <tzhao5@jhu.edu> and <hanliu@cs.jhu.edu>                      #
# Date: Jul 15th 2011                                                   #
# Version: 1.1.0                                                        #
#-----------------------------------------------------------------------#

huge.roc = function(path, theta, verbose = TRUE){
	gcinfo(verbose = FALSE)	
	ROC = list()
	
	theta = as.matrix(theta)
	d = ncol(theta)
	pos.total = sum(theta!=0)
	neg.total = d*(d-1) - pos.total
	
	if(verbose) cat("Computing F1 scores, false positive rates and true positive rates....")
	ROC$tp = rep(0,length(path))
   	ROC$fp = rep(0,length(path))
   	ROC$F1 = rep(0,length(path))
   	for (r in 1:length(path)){
   		tmp = as.matrix(path[[r]]) 
   		tp.all = (theta!=0)*(tmp!=0)
   		diag(tp.all) = 0
		ROC$tp[r] <- sum(tp.all!=0)/pos.total
		fp.all = (theta==0)*(tmp!=0)
		diag(fp.all) = 0
		ROC$fp[r] <- sum(fp.all!=0)/neg.total
		
		fn = 1 - ROC$tp[r]
		precision = ROC$tp[r]/(ROC$tp[r]+ROC$fp[r])
		recall = ROC$tp[r]/(ROC$tp[r]+fn)
		ROC$F1[r] = 2*precision*recall/(precision+recall)
		if(is.na(ROC$F1[r]))	ROC$F1[r] = 0
	}
	if(verbose) cat("done.\n")
		
	rm(precision,recall,tp.all,fp.all,path,theta,fn)
   	gc()	
		
	ord.fp = order(ROC$fp)
	
	tmp1 = ROC$fp[ord.fp]
	tmp2 = ROC$tp[ord.fp]
	par(mfrow = c(1,1))	
	plot(tmp1,tmp2,type="b",main = "ROC Curve", xlab = "False Postive Rate", ylab = "True Postive Rate",ylim = c(0,1))
	ROC$AUC = sum(diff(tmp1)*(tmp2[-1]+tmp2[-length(tmp2)]))/2
	
	rm(ord.fp, tmp1, tmp2)
	gc()
	class(ROC) = "roc"
	return(ROC)
}

print.roc = function(x, ...){
	cat("True Postive Rate: from",min(x$tp),"to",max(x$tp),"\n")
	cat("False Positive Rate: from",min(x$fp),"to",max(x$fp),"\n")
	cat("Area under Curve:",x$AUC,"\n")
	cat("Maximum F1 Score:",max(x$F1),"\n")
}

plot.roc = function(x, ...){	
	ord.fp = order(x$fp)
	par(mfrow = c(1,1))
	plot(x$fp[ord.fp],x$tp[ord.fp],type="b",main = "ROC Curve", xlab = "False Postive Rate", ylab = "True Postive Rate",ylim = c(0,1))
}