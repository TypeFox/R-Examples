### Function to do prais-winsten transformation
### Author: Konstantin Kashin
### August 1, 2013
prais.transform <- function(i,rhos,p,N.times,units,panel.vec,yX,obs.mat){
	mat <- matrix(NA, nrow=N.times, ncol=p)
	unit.i <- units[i]
	rho.i <- rhos[i]
	mat[obs.mat[i,],] <- yX[panel.vec %in% unit.i,]
	
	mat.L.mat <- rbind(c(mat[1,],rep(NA,p)),embed(mat,2))
	mat.diff <- mat.L.mat[,1:p]-rho.i*mat.L.mat[,(p+1):(2*p)]
	
	# Prais correction for start of runs
	begin.run <- which(!is.na(mat.L.mat[,1]) & is.na(mat.L.mat[,(p+1)]))
	mat.diff[begin.run,] <- mat[begin.run,]*(1-rho.i^2)^0.5
	
	# remove missing values
	keep.index <- sort(union(begin.run,which(complete.cases(mat.diff))))
	mat.diff <- mat.diff[keep.index,]
}