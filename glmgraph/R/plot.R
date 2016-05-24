plot.glmgraph <- function(x,...){

	penalty.factor <- which(x$penalty.factor!=0)+1
	z <- x$z +1
	nonzeros <- inds <- list()
	for(i in 1:x$nlambda2) nonzeros[[i]] <- inds[[i]] <- double(x$nlambda2)
    for(i in 1:x$nlambda2) nonzeros[[i]] <- which(apply(abs(x$betas[[i]]), 1, sum)!=0)    
    for(i in 1:x$nlambda2){
    	inds[[i]] <- intersect(penalty.factor, nonzeros[[i]])
    	if(length(z)>0) inds[[i]] <- intersect(inds[[i]],z)
    }
    
    n <- ceiling(sqrt(x$nlambda2))    
    par(mfrow=c(n,n))
    for(i in 1:x$nlambda2){
    	matplot(t(x$betas[[i]][inds[[i]],]) , type="l",xaxt="n",xlab=expression(lambda[1]), ylab="Coefficient", 
				main= substitute(lambda[2]==p, list(p=round(x$lambda2[i],digits=3))))
		axis(1,at=c(1:length(x$lambda1s[[i]])),labels=round(x$lambda1s[[i]], digits =3))
    }
}
















