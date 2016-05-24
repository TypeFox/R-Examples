Q.boot.ci <-
function(p,boot,alpha=.1)
{
	Q.hat <- sapply(as.list(p),function(h) distr(h,boot$model$dist,boot$model$par.hat,'q'))
	Q.boot <- sapply(as.list(p),function(h) distr(h,boot$model$dist,boot$par.star,'q'))
	colnames(Q.boot) <- paste('p = ',p,sep='')
	Q.ci.1 <- apply(Q.boot,2,quantile,c(alpha/2,.5,1-alpha/2))
	Q.ci.1[2,] <- Q.hat
	rownames(Q.ci.1) <- paste(c('low','Estimate','up'),'-',rownames(Q.ci.1),sep='')
	rownames(Q.ci.1)[2] <- 'Estimate'
	Q.ci.2 <- rbind(2*Q.ci.1[2,]-Q.ci.1[3,],Q.ci.1[2,],2*Q.ci.1[2,]-Q.ci.1[1,])
	rownames(Q.ci.2) <- rownames(Q.ci.1)
	list(percentile=Q.ci.1,reflexion=Q.ci.2)
}