
auc <- function(x, ...) UseMethod("auc")


auc.default <- function(x, ...) {
  stop(paste("Do not know how to handle objects of class", class(x)))
}



auc.numeric <- function(x, y, ...)
{
	sens <- c(x,0)
	spez <- c(y,1)
	fpr <- 1-spez
	fpr <- sort(fpr)
	tpr <- sort(sens)
	xdiff <- fpr[-1]-fpr[-length(fpr)]
	ydiff <- tpr[-1]-tpr[-length(tpr)]
	ans <- sum(xdiff*tpr[-length(tpr)]) + sum(0.5*xdiff*ydiff)
	ans
}



auc.Daim <- function(x, ...)
{
	auc.loob <- auc(x$roc$sensloob,x$roc$specloob)
	auc.app <- auc(x$roc$sensapp,x$roc$specapp)
	auc.samples <- sapply(x$sample.roc,function(y) auc(1-y[,1],y[,2]))
	if(class(x)[2] != "cv"){
		auc.632p <- auc(x$roc$sens632p,x$roc$spec632p)
		auc.632 <- auc(x$roc$sens632,x$roc$spec632)
		ans <- list(auc.632p=auc.632p, auc.632=auc.632, auc.loob=auc.loob, 
			auc.app=auc.app, auc.samples=auc.samples)
	}
	else{
		ans <- list(auc.loob=auc.loob, auc.app=auc.app, auc.samples=auc.samples)
	}
	ans
}

