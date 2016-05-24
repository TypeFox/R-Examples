# TODO: Add comment
# 
# Author: ianfellows
###############################################################################


rocplot<-function(logistic.model,diag=TRUE,pred.prob.labels=FALSE,prob.label.digits=3,AUC=TRUE){
	Spec <- Sens <- NULL
	##Adapted from lroc in epicalc
	firsttable <- table(logistic.model$fitted.values, logistic.model$y)
	colnames(firsttable) <- c("Non-diseased", "Diseased")
	firsttable1 <- cbind(as.numeric(rownames(firsttable)), firsttable)
	rownames(firsttable1) <- rep("", nrow(firsttable1))
	colnames(firsttable1)[1] <- "predicted.prob"
	secondtable <- firsttable
	for (i in 1:length(secondtable[, 1])) {
		secondtable[i, 1] <- (sum(firsttable[, 1]) - sum(firsttable[(1:i), 
									1]))/sum(firsttable[, 1])
		secondtable[i, 2] <- (sum(firsttable[, 2]) - sum(firsttable[(1:i), 
									2]))/sum(firsttable[, 2])
	}
	secondtable <- as.data.frame(rbind((c(1, 1)), secondtable))
	colnames(secondtable) <- c("Spec", "Sens")
	rownames(secondtable)[1] <- "0"
	model.des <- paste("logit (", deparse(logistic.model$formula), 
			")", sep = "")
	auc <- 0
	for (i in 1:(nrow(secondtable) - 1)) {
		auc <- auc + (secondtable[i, 1] - secondtable[(i + 1), 
							1]) * 0.5 * (secondtable[i, 2] + secondtable[(i + 
										1), 2])
	}
	p<-qplot(Spec, Sens,data=secondtable[nrow(secondtable):1,], xlab = "1-Specificity", 
			ylab = "Sensitivity",main=model.des, xlim = (c(0, 1)), ylim = (c(0, 1)))+geom_step()
	if(diag)
		p<-p+geom_segment(x=0,y=0,xend=1,yend=1,colour="red",size=1.5)
	if(pred.prob.labels)
		p<-p+geom_text(size=3,angle=-45,hjust=-.3,aes(x=Spec,y=Sens,
						label=formatC(as.numeric(rownames(secondtable)),
								digits=prob.label.digits)))
	if(AUC)
		p<-p+annotate("text",x=.8,y=0.05,label=paste("AUC=",formatC(auc,digits=4)))
	p
}
