Q.conf.int <-
function(p,model,alpha=.1,ln=FALSE)
{
	Q.hat <- distr(p,model=model,type='q')
	out <- sapply(as.list(p),function(h)
		distr(c(alpha/2,1-alpha/2),'norm',as.numeric(delta.Q(h,model,ln)),type='q'))
	if(ln) out <- exp(out)
	r.names <- paste(c('low-','up-'),c(alpha/2,1-alpha/2)*100,'%',sep='')
	r.names <- c(r.names[1],'Estimate',r.names[2])
	out <- rbind(out[1,],Q.hat,out[2,])
	dimnames(out) <- list(r.names,paste('p = ',p,sep=''))
	return(out)
}