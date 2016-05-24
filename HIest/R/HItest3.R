HItest3 <-
function(class,MLE,thresholds=c(2,8)){
	LLC <- apply(class[,1:15],1,max)
	LLM <- MLE$logLik
	dAIC <- (2*1-2*LLC) - (2*5-2*LLM)
	c1 <- class$LLD>thresholds[1]
	c2 <- (LLM-LLC)<thresholds[2]
	data.frame(Best.class = class$Best,LL.class=LLC,LLD.class=class$LLD,LL.max=MLE$logLik,dAIC,c1,c2)
	}
