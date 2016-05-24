timedelay.lm.batch <-
function(bspline.data, expr.data, regulator.list, target.list=rownames(bspline.data), ...) {
	
	target<- as.character(intersect(target.list,rownames(bspline.data)))
	regulator<- as.character(intersect(regulator.list,rownames(bspline.data))) 
	result<-NULL  	
	c<-0   	
	for (i in 1:length(target)) {
		if ((i*100)%/%length(target)>=c) {
			cat(c,'% done.\n',sep='')
			c<-c+1
		}		
		result.i<-timedelay.lm(bspline.data, expr.data, target[i], setdiff(regulator,target[i]),...)
		if (!is.null(result.i)) {
			result<-rbind(result,data.frame(regulator=names(result.i$delay),target=rep(target[i],length(result.i$delay)),
      coef=summary(result.i$fit)$coef[,1],delay=result.i$delay,adj.r.squared=summary(result.i$fit)$adj.r.squared))
		}
	}
	
	rownames(result)<-NULL
	return(result)
}

