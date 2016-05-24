# Description: 	Implementation of Rubin's SIR
# Usage:	Sir(data.mat)



sir <- function(data.mat,theta.vector,theta.mat,M,m,tol=1e-06,ll.func,df=0)  {
    importance.ratio <- rep(NA,M)
    rand.draw    <- rmultnorm(M,theta.vector,theta.mat,tol = 1e-04)
    if (df > 0) 
        rand.draw <- rand.draw/(sqrt(rchisq(M,df)/df))
    
    empirical.draw.vector <- apply(rand.draw,1,ll.func,data.mat)
    if (sum(is.na(empirical.draw.vector)) == 0)  {
	print("SIR: finished generating from posterior density function")
    	print(summary(empirical.draw.vector))
    }
    else {
	print(paste("SIR: found",sum(is.na(empirical.draw.vector)),
		"NA(s) in generating from posterior density function, quiting"))
	return()
    }

    if (df == 0)  {
        normal.draw.vector <- apply(rand.draw,1,normal.posterior.ll,data.mat)
    }
    else {
	theta.mat <- ((df-2)/(df))*theta.mat
        normal.draw.vector <- apply(rand.draw,1,t.posterior.ll,data.mat,df)
    }
    if (sum(is.na(normal.draw.vector)) == 0)  {
	print("SIR: finished generating from approximation distribution")
    	print(summary(normal.draw.vector))
    }
    else {
	print(paste("SIR: found",sum(is.na(normal.draw.vector)),
		"NA(s) in generating from approximation distribution, quiting"))
	return()
    }

    importance.ratio <- exp(empirical.draw.vector - normal.draw.vector)
    importance.ratio[is.finite=F] <- 0
    importance.ratio <- importance.ratio/max(importance.ratio)
    if (sum(is.na(importance.ratio)) == 0)  {
	print("SIR: finished calculating importance weights")
    	print(summary(importance.ratio))
    }
    else  {
	print(paste("SIR: found",sum(is.na(importance.ratio)),
		"NA(s) in calculating importance weights, quiting"))
	return()
    }

    accepted.mat <- rand.draw[1:2,]
    while(nrow(accepted.mat) < m+2)  {
	rand.unif <- runif(length(importance.ratio))
	accepted.loc <- seq(along=importance.ratio)[(rand.unif-tol) <= importance.ratio]
	rejected.loc <- seq(along=importance.ratio)[(rand.unif-tol) > importance.ratio]
	accepted.mat <- rbind(accepted.mat,rand.draw[accepted.loc,])
	rand.draw <- rand.draw[rejected.loc,]
	importance.ratio <- importance.ratio[rejected.loc]
	print(paste("SIR: cycle complete,",(nrow(accepted.mat)-2),"now accepted"))
    }
    accepted.mat[3:nrow(accepted.mat),]
}

