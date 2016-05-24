##' Estimate the unrestricted covariance matrix using the Pearson-Lawley's Multivariate Correction
##'	
##' @param data Either a data matrix containing the restricted dataset, or a variance covariance matrix of the selected covariances. 
##' The variables used for selected (the Xs) must be in the first \code{p}
##' columns, while the "incidental" variables are in the last n-p:p columns
##' @param p the number of columns used for selection
##' @param v.pp The variance-covariance matrix of the population selection variables
##' @author Dustin Fife
##' @importFrom fifer cor2cov
##' @export
##' @references Birnbaum, Z. W., Paulson, E., & Andrews, F. C. (1950). On the effect of selection performed on some coordinates of a multi-dimensional population. Psychometrika, 15(2), 191-204.
##' @examples
##' # do a simulation to demonstrate unbiasedness of it
##'	#### specify population correlation matrix
##'	cor.mat = matrix(c(
##'		1, .2, .4, .5,
##'		.2, 1, .3, .6,
##'		.4, .3, 1, .4,
##'		.5, .6, .4, 1), ncol=4)
##'	
##'		#### sample from that matrix
##'	d = data.frame(mvrnorm(100000, mu=rep(0, times=ncol(cor.mat)), Sigma=cor.mat)	)
##'	names(d) = c("X1", "X2", "Y1", "Y2")
##'	
##'		#### restrict the data according to X1 and X2
##'	d.r = d[which((d$X1)>0 & d$X2>0),]
##'	
##'	cor.mat.corrected = mv.correction(data=d.r, p=2, v.pp = cor.mat[1:2, 1:2])
##'	
##' # compare the two estimates
##'	cor.mat
##'	round(cor.mat.corrected, digits=4)
##' # use example I in Birnbaum, Paulson, and Andrews (1950)
##' ## input population correlation matrix
##' cor = matrix(c(1, .703, .527, .499, .777,
##'				.703, 1, .404, .555, .679,
##'				.527, .404, 1, .253, .578,
##'				.499, .555, .253, 1, .427,
##'				.777, .679, .578, .427, 1), nrow=5)
##' require(fifer) ### for cor2cov function
##' v = cor2cov(cor, sd=c(9.93, 9.39, 8.80, 7.19, 8.05))
##' ### change order so selection variables occur first
##' ord = c(5,1,2,3,4)
##' v = v[ord, ord]
##'
##' # input observed correlation matrix
##' v.star = matrix(c(
##' 	43, 37.65, 37.24, 25.59, 14.46,
##' 	37.65, 72.24, 48.95, 29.04, 23.90,
##' 	37.24, 48.95, 82.60, 20.10, 31.24,
##' 	25.59, 29.04, 20.10, 68.04, 4.71,
##' 	14.46, 23.90, 31.24, 4.71, 48.31), ncol=5)
##' mv.correction(v.star, p=2, v.pp = v[1:2, 1:2])
mv.correction = function(data, p, v.pp){

	### check if dim(v.pp) = p	
	if(length(v.pp)==1){
		if (p!=1){
			stop("The dimensions of v.pp do not equal p")			
		}
	} else {
		if (dim(v.pp)[1] != p | dim(v.pp)[2] != p) {
			stop("The dimensions of v.pp do not equal p")	        
	    }
    }


	### if they supply a data matrix, convert to a variance covariance matrix
	if (ncol(data) != nrow(data)){
		data = cov(data.matrix(data))
	}
	var.covar = data
	
	### get sufficient estimates
	n = ncol(var.covar)
	vstar.pp = var.covar[1:p, 1:p]
    vstar.np.p = matrix(var.covar[(p + 1):n, 1:p], nrow=length((p+1):n), ncol=length(1:p))
    vstar.np.np = matrix(var.covar[(p + 1):n, (p + 1):n], nrow=length((p+1):n))

	### output actual estimates
	correction.off.diag = vstar.np.p %*% solve(vstar.pp) %*% v.pp
	correction.diag = vstar.np.np - vstar.np.p%*%(solve(vstar.pp) - solve(vstar.pp)%*%v.pp%*%solve(vstar.pp))%*%t(vstar.np.p)
	v.prime = rbind(v.pp, correction.off.diag)
	v.prime2 = rbind(t(correction.off.diag), correction.diag)
	v.prime = cbind(v.prime, v.prime2)
	v.prime
}