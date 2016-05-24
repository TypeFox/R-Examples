SubLikelihoods <- function(x,subject=c(1:x$nsubj)){
	Lframe <- as.data.frame(cbind(round(x$evalpoints,4),as.data.frame(x$RCC)[,subject])) 
	colnames(Lframe) <- c("Theta",paste("Subject:",subject,sep=" "))
	return(Lframe)
}



